// SPDX-License-Identifier: GPL-3.0-only
#include "engine.hpp"
#include "dsp_fft.hpp"
#include "formant.hpp"
#include "fractional_sinc.hpp"
#include "glottal.hpp"
#include "harmonic_mask.hpp"
#include "harmonic_phase.hpp"
#include "harmonic_synth.hpp"
#include "pitch.hpp"
#include "source_analysis.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <random>
#include <stdexcept>
#include <utility>

namespace tn {
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr std::size_t max_samples = 16'000'000;
using Clock = std::chrono::steady_clock;

double milliseconds(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

double unit(double x) { return std::clamp(x, 0.0, 1.0); }
double mix(double a, double b, double x) { return a + (b - a) * x; }
double smooth(double x) {
    x = unit(x);
    return x * x * (3 - 2 * x);
}

using detail::Sinc;

int stretch_hop(int sr) { return std::max(1, static_cast<int>(std::lround(sr * .012))); }

struct TimeMap {
    double ni, no, front_i, front_o, front_nominal;
    int source_begin = 0;  // offset of the selected samples in the analysis buffer
    double loop_a = 0, loop_b = 0, loop_out = 0;
    int passes = 1;
    bool hold_end = false;
    bool identity() const {
        // An unpitched, fully fixed render can still be padded without
        // changing the source clock. A pitched sustain holds the end model.
        return !hold_end && (ni == no || (front_i == ni && no >= ni)) && front_i == front_nominal;
    }
    bool looping() const { return passes > 1; }
    bool needs_sustain_guard(int hop) const {
        const double tail = ni - front_i, duration = no - front_o;
        // Keep the original end-support trigger and its analysis context.
        // Residual replacement is controlled separately.
        return hold_end || (tail > 0. && tail <= 4.*hop+2. && duration >= 4.*tail && duration-tail >= 2.*hop);
    }
    double front_source() const { return source_begin + front_i; }
    int last_source() const { return source_begin + static_cast<int>(ni) - 1; }
    int active_samples() const {
        return static_cast<int>(!hold_end && front_i == ni ? std::min(no, std::ceil(front_o)) : no);
    }
    double operator()(double z) const {
        const double last = last_source();
        if (z < front_o && front_nominal > 0) return std::min(last, source_begin + z * front_i / front_nominal);
        if (identity()) return std::min(last, source_begin + z);
        const double tail_i = ni - front_i, tail_o = no - front_o, u = z - front_o;
        if (tail_i <= 0 || tail_o <= 0) return last;
        double src;
        if (!looping()) {
            src = front_source() + u * tail_i / tail_o;
        } else {
            const double intro = loop_a - front_source();
            if (u < intro) return front_source() + u;
            if (u >= intro + loop_out) return std::min(last, loop_b + u - intro - loop_out);
            const double width = loop_b - loop_a;
            const double phase = std::fmod((u - intro) * passes / loop_out, 2.);
            src = loop_a + width * (phase <= 1. ? phase : 2. - phase);
        }
        return std::clamp(src, double(source_begin), last);
    }
};

// Extend the stable sustain only, never replay the attack or the release.
// The ping-pong is a cursor over forward grains, not reversed PCM.
void choose_loop(TimeMap& t, const std::vector<double>& marks, const detail::HarmonicSplit& split,
                 int sr, bool force) {
    const double extra = t.no - t.front_o - (t.ni - t.front_i);
    if (extra < (force ? .030 : .085) * sr) return;
    double best = 0;
    const auto stable = [&](std::size_t i) { return split.strength[i] > .72f && split.periodic_share[i] > .50f; };
    for (std::size_t i = 0; i < split.strength.size();) {
        if (!stable(i)) { ++i; continue; }
        const auto start = i++;
        while (i < split.strength.size() && stable(i)) ++i;
        const double lo = std::max({t.front_source() + .006 * sr, marks.front() + .012 * sr,
                                    start * split.hop + .012 * sr});
        const double hi = std::min({t.source_begin + t.ni - .010 * sr, marks.back() - .012 * sr,
                                    (i-1) * split.hop - .012 * sr});
        const double width = std::min({hi - lo, .160 * sr, extra * .49});
        if (width < (force ? .012 : .035) * sr) continue;
        const double center = (lo + hi) / 2;
        const double a = center - width / 2, b = center + width / 2;
        double prior = 0., cost = 0.;
        int count = 0;
        auto it = std::lower_bound(marks.begin(), marks.end(), a);
        for (; it != marks.end() && *it < b; ++it) {
            if (it == marks.begin()) continue;
            const double p = *it - *(it - 1);
            if (prior > 0) cost += std::abs(std::log(p / prior));
            prior = p;
            ++count;
        }
        if (count < 3 || cost / count > .13) continue;
        const double score = width / (1. + 25. * cost / count);
        if (score > best) {
            best = score;
            t.loop_a = a;
            t.loop_b = b;
            t.loop_out = width + extra;
            t.passes = 1 + 2 * static_cast<int>(extra / (2 * width));
        }
    }
}

pitch::Track rebase_source_track(const pitch::Track& track, std::size_t start, std::size_t count) {
    pitch::Track local;
    local.step = track.step;
    const std::size_t last = count-1;
    const std::size_t frames = 1+(last+local.step-1)/local.step;
    local.frames.reserve(frames);
    for (std::size_t i = 0; i < frames; ++i)
        local.frames.push_back(pitch::sample(track, static_cast<double>(start+std::min(i*local.step, last))));
    return local;
}

pitch::Track source_track(const detail::SourceAnalysisInput& in) {
    const auto& x = in.segment;
    const auto sr = in.rate, start = in.first;
    if (!in.mrq_f0.frames.empty()) {
        if (start == 0 && x.size() == in.recording.size()) return in.mrq_f0;
        return rebase_source_track(in.mrq_f0, start, x.size());
    }
    if (in.f0.empty() && in.forced_hz == 0) {
#ifdef TN_FINDS_NEBULA_F0
        // OpenMoresampler estimates on the complete recording. Rebase the
        // resulting F0 to the OTO analysis buffer in source samples.
        auto track = pitch::estimate_source(in.recording, sr);
        if (start == 0 && x.size() == in.recording.size()) return track;
        return rebase_source_track(track, start, x.size());
#else
        return pitch::estimate(x, sr);
#endif
    }
    pitch::Track tr;
    tr.step = std::max(1, static_cast<int>(std::lround(sr * .0025)));
    const std::size_t count = (x.size() - 1 + tr.step - 1) / tr.step + 1;
    tr.frames.resize(count);
    std::size_t j = 0;
    for (std::size_t i = 0; i < count; ++i) {
        double hz = in.forced_hz;
        if (hz == 0) {
            const double time = (start + static_cast<double>(i) * tr.step) / sr;
            const auto& p = in.f0;
            if (time <= p.front().time_s) hz = p.front().hz;
            else if (time >= p.back().time_s) hz = p.back().hz;
            else {
                while (j + 1 < p.size() && p[j+1].time_s < time) ++j;
                const double u = (time - p[j].time_s) / (p[j+1].time_s - p[j].time_s);
                if (p[j].hz && p[j+1].hz) hz = mix(p[j].hz, p[j+1].hz, u);
                else hz = u < .5 ? p[j].hz : p[j+1].hz;
            }
        }
        if (hz > 0) tr.frames[i] = {hz, 1.};
    }
    return tr;
}

detail::SourceAnalysis analyze_source(const detail::SourceAnalysisInput& in) {
    detail::SourceAnalysis a;
    if (in.forced_hz < 0.) return a;
    const auto& x = in.segment;
    const int sr = in.rate;
    auto anchors = source_track(in);
    const bool mrq = !in.mrq_f0.frames.empty();
    const bool supplied = !in.f0.empty() || in.forced_hz > 0.;
    const auto voiced = [](pitch::Frame f) { return f.hz > 0.; };
    if (mrq && std::none_of(anchors.frames.begin(), anchors.frames.end(), voiced)) {
        // An all-unvoiced MRQ (including an OTO slice of a voiced file) is
        // still an authoritative track, not a reason to estimate F0.
        a.clock = std::move(anchors);
        return a;
    }
#ifdef TN_FINDS_NEBULA_F0
    // Nebula and MRQ both provide the voicing and frequency of each frame.
    const bool bridge = !supplied;
#else
    const bool bridge = mrq;
#endif
    // The bridged clock only supplies phase continuity through UV spans.
    // The original track still supplies the separator's voicing evidence.
    a.clock = bridge ? pitch::bridge_clock(anchors) : pitch::candidate_clock(x, sr, anchors, supplied);
    auto support = pitch::locate(x, sr, a.clock);
    if (support.empty()) return a;
    const std::vector<pitch::Island>* located = nullptr;
    // With no UV frames, bridging leaves the track and its pulse locations unchanged.
    const bool fully_voiced = bridge && std::all_of(anchors.frames.begin(), anchors.frames.end(), voiced);
    if (fully_voiced) located = &support;
    a.split = detail::split_harmonic(x, sr, support.front().marks, a.clock, in.cleanup,
                                    in.mask_return, supplied ? nullptr : &anchors, located,
                                    mrq && !fully_voiced);
    a.pulses = std::move(support.front());
    if (in.need_rd && !a.split.harmonic.empty())
        a.rd = detail::estimate_glottal_track(a.split.harmonic, a.pulses.marks, a.split, sr);
    return a;
}

double center_pitch(const pitch::Track& tr, const TimeMap& t) {
    std::vector<double> hz;
    std::size_t prefix = static_cast<std::size_t>(t.front_source() / tr.step);
    const auto first = static_cast<std::size_t>(t.source_begin / tr.step);
    const auto end = std::min(tr.frames.size(), std::size_t((t.last_source()+tr.step-1) / tr.step + 1));
    for (std::size_t i = prefix; i < end; ++i) {
        if (tr.frames[i].hz > 0) hz.push_back(tr.frames[i].hz);
    }
    if (hz.empty()) {
        for (std::size_t i = first; i < end; ++i)
            if (tr.frames[i].hz > 0) hz.push_back(tr.frames[i].hz);
    }
    if (hz.empty()) return 0;
    auto mid = hz.begin() + hz.size() / 2;
    std::nth_element(hz.begin(), mid, hz.end());
    const double median = *mid;
    double sum = 0, w = 0;
    for (auto it = tr.frames.begin() + std::min(prefix, end); it != tr.frames.begin() + end; ++it) {
        if (it->hz <= 0) continue;
        const double d = std::log(it->hz / median);
        const double weight = std::max(.1, it->confidence) * std::exp(-d*d*20);
        sum += std::log(it->hz) * weight;
        w += weight;
    }
    return w > 0 ? std::exp(sum / w) : median;
}

// Short forward-reading WSOLA grains. The search is bounded near the requested
// map; an uncorrelated join uses approximately equal-power overlap, not a dip.
struct Warp {
    std::vector<float> dry, harmonic, residual_correction;
};

Warp stretch(const std::vector<float>& x, const detail::HarmonicSplit& split,
             const std::vector<float>& correction, const TimeMap& map, int sr, int count) {
    const auto& h = split.harmonic;
    if (map.identity()) return {x, h, correction};
    Warp wave{std::vector<float>(count), std::vector<float>(h.empty() ? 0 : count),
              std::vector<float>(correction.empty() ? 0 : count)};
    const int begin = map.source_begin, last = map.last_source();
    if (count <= 1 || begin == last) {
        std::fill_n(wave.dry.begin(), count, x[begin]);
        if (!h.empty()) std::fill_n(wave.harmonic.begin(), count, h[begin]);
        if (!correction.empty())
            std::fill_n(wave.residual_correction.begin(), count, correction[begin]);
        return wave;
    }
    const int hop = stretch_hop(sr);
    const int search = std::max(1, static_cast<int>(std::lround(sr * .004)));
    const int stride = std::max(1, sr / 8000);
    const auto residual = [&](int at) {
        const int p = std::clamp(at, begin, last);
        return double(x[p])-(h.empty() ? 0.f : h[p]);
    };
    std::vector<double> reference, candidates, cross, energy;
    reference.reserve((hop+stride-1)/stride);
    // At velocity 100, no WSOLA overlap should reach back into the fixed
    // consonant: a sharp attack cannot move just because the vowel stretches.
    const int keep = map.front_i == map.front_nominal ? static_cast<int>(map.front_o) : 0;
    std::copy_n(x.begin()+begin, keep, wave.dry.begin());
    if (!h.empty()) std::copy_n(h.begin()+begin, keep, wave.harmonic.begin());
    if (!correction.empty()) std::copy_n(correction.begin()+begin, keep, wave.residual_correction.begin());
    int a = begin+keep;
    for (int i = keep; i < count; i += hop) {
        const int span = std::min(hop, count - i);
        const int want = static_cast<int>(std::lround(map(i + span)));
        int b = want;
        double agreement = (a + span == b) ? 1. : 0.;
        if (a + span != b && span >= hop / 2) {
            const double periodic_share = split.share(map(i + span * .5));
            reference.clear();
            double aa = 0.;
            for (int j = stride/2; j < span; j += stride) {
                const double v = residual(a+j);
                reference.push_back(v);
                aa += v*v;
            }
            const auto score = [&](const double* q) {
                double ab = 0., bb = 0.;
                for (std::size_t k = 0; k < reference.size(); ++k) {
                    const double w = q[k];
                    ab += reference[k]*w;
                    bb += w*w;
                }
                return aa*bb > 1.e-28 ? ab/std::sqrt(aa*bb) : 0.;
            };
            candidates.resize(reference.size());
            for (std::size_t k = 0; k < candidates.size(); ++k)
                candidates[k] = residual(want-span+stride/2+static_cast<int>(k)*stride);
            agreement = score(candidates.data());
            if (periodic_share < .92 && want > begin+search && want <= last-search) {
                // Every trial is a view of the same strided residual. Read
                // and clamp it once; the scoring and tie-breaking stay the same.
                const int trials = 1+2*search/stride;
                candidates.resize(reference.size()+trials-1);
                for (std::size_t k = 0; k < candidates.size(); ++k)
                    candidates[k] = residual(want-search-span+stride/2+static_cast<int>(k)*stride);
                // Vectorize across independent delay trials rather than a
                // sum over time. Each trial retains its original sum order.
                cross.assign(trials, 0.);
                energy.assign(trials, 0.);
                for (std::size_t k = 0; k < reference.size(); ++k) {
                    const double v = reference[k];
                    const auto* q = candidates.data()+k;
                    for (int trial = 0; trial < trials; ++trial) {
                        const double w = q[trial];
                        cross[trial] += v*w;
                        energy[trial] += w*w;
                    }
                }
                double best = agreement;
                for (int trial = 0, d = -search; trial < trials; ++trial, d += stride) {
                    const double norm = aa*energy[trial];
                    const double corr = norm > 1.e-28 ? cross[trial]/std::sqrt(norm) : 0.;
                    const double cost = corr - (.12 + .60 * periodic_share) * std::abs(d) / search;
                    if (cost > best) { best = cost; b = want + d; agreement = corr; }
                }
            }
        }
        agreement = unit(agreement);
        for (int j = 0; j < span; ++j) {
            const int p = std::clamp(a + j, begin, last), q = std::clamp(b + j - span, begin, last);
            const double w = smooth(j / static_cast<double>(span)), v = 1. - w;
            const double gain = p == q ? 1. : 1. / std::sqrt(v*v + w*w + 2*agreement*v*w);
            wave.dry[i + j] = static_cast<float>((x[p] * v + x[q] * w) * gain);
            if (!h.empty()) wave.harmonic[i + j] = static_cast<float>((h[p] * v + h[q] * w) * gain);
            if (!correction.empty())
                wave.residual_correction[i + j] = static_cast<float>((correction[p] * v + correction[q] * w) * gain);
        }
        a = b;
    }
    return wave;
}

// The spectral windows give consonants their color but can spread a short
// unvoiced attack into nearby silence. Follow the original short-time level
// in the fixed region without importing its waveform or boosting voiced leakage.
template<class Residual>
void follow_consonant_envelope(std::vector<float>& noise, const Residual& sample,
                               const detail::HarmonicSplit& split, const TimeMap& map, int edge, int sr) {
    const int hop = std::max(1, static_cast<int>(std::lround(sr*.0015)));
    const auto level = [&](const auto& read, int at, int first, int last) {
        const int a = std::max(first, at-hop), b = std::min(last, at+hop+1);
        double sum = 0., energy = 0.;
        for (int i = a; i < b; ++i) {
            const double v = read(i);
            sum += v;
            energy += v*v;
        }
        const double mean = sum/(b-a);
        return std::sqrt(std::max(0., energy/(b-a)-mean*mean));
    };
    const auto read_noise = [&](int i) { return noise[i]; };
    const int frames = (edge+hop-1)/hop;
    std::vector<double> gains(frames+1, 1.);
    for (int j = 0; j < frames; ++j) {
        const int output = j*hop;
        const int source = static_cast<int>(std::lround(map(output)));
        const double uv = 1.-smooth((split.activity(source)-.10)/.25);
        if (uv <= 0.) continue;
        const double measured = level(sample, source, map.source_begin, map.last_source()+1);
        const double generated = level(read_noise, output, 0, edge);
        if (generated < 1.e-12) continue;
        const double ratio = measured/generated;
        // Leave stationary texture alone within a small tolerance. Boost
        // only short, measured unvoiced events, not missing high harmonics.
        const double g = ratio < .77 ? std::min(1., 1.30*ratio) : std::clamp(ratio/1.30, 1., 2.2);
        gains[j] = mix(1., g, uv*smooth((edge-output)/double(2*hop)));
    }
    for (int j = 0; j < frames; ++j) {
        const int a = j*hop, b = std::min(edge, a+hop);
        for (int i = a; i < b; ++i)
            noise[i] *= static_cast<float>(mix(gains[j], gains[j+1], smooth(double(i-a)/(b-a))));
    }
}

// Generate independent noise from a time-varying spectrum at the mapped
// source position, without stretching or repeating the recorded waveform.
std::vector<float> synthesize_residual(const std::vector<float>& x, const detail::HarmonicSplit& split,
                                      const std::vector<float>& correction, const pitch::Track& clock,
                                      const TimeMap& map, int sr, int first, int edge, int stop, int count,
                                      bool follow_fixed_envelope) {
    const auto& harmonic = split.harmonic;
    const auto residual = [&](int i) {
        return (x[i] - (harmonic.empty() ? 0.f : harmonic[i])) - (correction.empty() ? 0.f : correction[i]);
    };
    const int n = std::clamp(1 << static_cast<int>(std::lround(std::log2(sr*.046))), 256, 8192);
    const int half = n/2, radius = std::min(half, static_cast<int>(std::lround(sr*.016)));
    const double df = double(sr)/n;
    const double left = std::min(map.front_source(), double(map.last_source()));
    const double width = map.last_source()-left;
    const double fixed_width = left-map.source_begin;
    const int fixed_frames = static_cast<int>(std::ceil(fixed_width/(sr*.006)));
    const int tail_frames = std::max(1, static_cast<int>(std::ceil(width/(sr*.006))));
    const int intervals = fixed_frames+tail_frames;
    const int input_size = static_cast<int>(x.size());
    const int join_center = static_cast<int>(std::lround(left))-radius;
    const double join_weight = join_center >= map.source_begin && split.activity(join_center) > .25
        ? 1.-smooth(width/radius) : 0.;
    const detail::Fft fft(n);
    std::vector<detail::Complex> spec(n);
    std::vector<float> window(2*radius+1), smoothed(half+1);
    std::vector<double> sum(half+1);
    struct PowerFrame { int index = -1, center = -1; std::vector<float> power; };
    std::array<PowerFrame, 4> powers;
    for (int j = 0; j <= 2*radius; ++j)
        window[j] = static_cast<float>(.5-.5*std::cos(pi*j/radius));
    const auto smooth_power = [&](std::vector<float>& power, int bins) {
        sum[0] = 0.;
        for (int k = 1; k <= half; ++k) sum[k] = sum[k-1]+power[k];
        for (int k = 1; k <= half; ++k) {
            const int a = std::max(1, k-bins), b = std::min(half, k+bins);
            smoothed[k] = static_cast<float>((sum[b]-sum[a-1])/(b-a+1));
        }
        power.swap(smoothed);
    };
    // Cache frames on demand; replacement may cover any source length.
    const auto power_at = [&](int frame) -> const std::vector<float>& {
        auto& cached = powers[frame % powers.size()];
        auto& power = cached.power;
        if (cached.index == frame) return power;
        cached.index = frame;
        const double position = frame < fixed_frames
            ? map.source_begin+fixed_width*frame/fixed_frames
            : left+width*(frame-fixed_frames)/tail_frames;
        const int center = static_cast<int>(std::lround(position));
        // Use a one-sided window at a voiced cutoff if the inner part is
        // still voiced. It avoids estimating aspiration from the last,
        // least well separated pulse while keeping UV transitions local.
        const int inner = std::max(0, center-radius/2);
        const bool inner_voice = split.activity(center) > .25 && split.activity(inner) > .25;
        int analysis = inner_voice ? std::max(inner, std::min(center, input_size-1-radius/2)) : center;
        if (inner_voice && join_weight > 0.) {
            // With less than half a window of stable vowel left, trust the
            // complete window just before the join more than a truncated
            // last cycle. Do not pull a genuine UV transition back in time.
            const double blend = join_weight*smooth((center-join_center)/double(radius));
            analysis = static_cast<int>(std::lround(mix(analysis, std::min(analysis, join_center), blend)));
        }
        // Near a short tail, distinct grid points can map to the same actual
        // one-sided analysis window. Keep the same PSD in the requested slot
        // so the pair of adjacent-frame references below remains stable.
        if (cached.center == analysis) return power;
        cached.center = analysis;
        for (const auto& other : powers) {
            if (&other != &cached && other.center == analysis) {
                power = other.power;
                return power;
            }
        }
        power.assign(half+1, 0.f);
        const int start = analysis-radius;
        const int begin = std::max(0, -start), end = std::min(2*radius, input_size-start);
        const auto frame_pitch = pitch::sample(clock, analysis);
        const double hz = frame_pitch.hz;
        const bool voiced = split.activity(analysis) > .20 && hz > 0.;
        // At an actual cutoff the separator has fewer right-hand periods;
        // the residual can contain a stronger copy of the vowel. Estimate
        // the noise floor beside each partial before smoothing the envelope.
        const auto apodize = [&](int j) {
            // For a truncated half window, taper the new edge instead of
            // transforming a hard join to the padding as broadband noise.
            const double edge = std::max(1, radius/4);
            return window[j]*smooth((end-j-1)/edge)*smooth((j-begin)/edge);
        };
        double energy = 0., weight = 0., average = 0.;
        for (int j = begin; j < end; ++j) {
            const double w = apodize(j);
            energy += w*w;
            weight += w;
            average += w*residual(start+j);
        }
        if (energy <= 1.e-12) return power;
        average /= weight;
        std::fill(spec.begin(), spec.end(), detail::Complex{});
        for (int j = begin; j < end; ++j)
            spec[j] = {static_cast<float>(apodize(j)*(residual(start+j)-average)), 0.f};
        fft.run(spec);
        // Average out the harmonic comb of any leakage in the separation;
        // preserve broad breath/frication bands, not peaks at the old F0.
        const int bins = std::clamp(static_cast<int>(std::lround(std::max(170., .6*hz)/df)), 1, half-1);
        for (int k = 1; k <= half; ++k) power[k] = static_cast<float>(std::norm(spec[k])/energy);
        if (voiced) {
            // Bound the lower-band peaks by power halfway between adjacent
            // source harmonics. This affects noise estimation only, never H.
            const int limit = std::min(half, static_cast<int>(3000./df));
            for (int k = 1; k <= limit; ++k) {
                const double order = k*df/hz;
                const double left_skirt = std::floor(order-.5)+.5;
                const auto bin = [&](double at) {
                    return std::clamp(static_cast<int>(std::lround(at*hz/df)), 1, half);
                };
                const double floor = mix(power[bin(left_skirt)], power[bin(left_skirt+1.)], order-left_skirt);
                smoothed[k] = std::min(power[k], static_cast<float>(1.5*floor));
            }
            std::copy_n(smoothed.begin()+1, limit, power.begin()+1);
        }
        smooth_power(power, bins);
        smooth_power(power, bins);
        for (int k = 1; k <= half; ++k) power[k] *= static_cast<float>(smooth((k*df-45.)/75.));
        return power;
    };

    std::uint32_t seed = static_cast<std::uint32_t>(count) ^ (static_cast<std::uint32_t>(edge)*16777619u);
    const int stride = std::max(1, (map.last_source()-map.source_begin)/64);
    for (int i = map.source_begin; i <= map.last_source(); i += stride)
        seed = (seed ^ static_cast<std::uint32_t>(std::lround(x[i]*32768.)))*16777619u;
    std::mt19937 fixed_random(seed ^ 0x9e3779b9u), tail_random(seed);
    detail::RealInverseFft inverse(n);
    std::vector<detail::Complex> excitation(half+1);
    std::vector<float> noise(stop-first), block, bell(n);
    for (int j = 0; j < n; ++j) bell[j] = static_cast<float>(std::sin(pi*j/n));
    // Use the same grid, overlap and excitation past the boundary whether
    // or not the optional fixed-region noise was needed before it.
    const int origin = edge-(edge-first+half-1)/half*half;
    for (int center = origin; center-half < stop; center += half) {
        auto& random = center < edge ? fixed_random : tail_random;
        const auto draw = [&] { return static_cast<float>(random() >> 8)*(2.f/16777216.f)-1.f; };
        const double source = map(center);
        const double at = fixed_frames > 0 && source < left
            ? fixed_frames*unit((source-map.source_begin)/fixed_width)
            : fixed_frames+(width > 0. ? tail_frames*unit((source-left)/width) : 0.);
        const int i = std::min(intervals-1, static_cast<int>(at));
        const double frac = at-i;
        const auto& a = power_at(i);
        const auto& b = power_at(i+1);
        for (int k = 1; k < half; ++k) {
            // Each uniform component has variance 1/3: E[|X[k]|²] = n * PSD.
            const float gain = static_cast<float>(std::sqrt(1.5*n*mix(a[k], b[k], frac)));
            const float real = draw(), imag = draw();
            excitation[k] = detail::Complex(real, imag)*gain;
        }
        inverse.run(excitation, block, 1.f/n);
        const int start = center-half;
        // Independent blocks use equal-power overlap. A unit-sum amplitude
        // window would audibly modulate the noise at the synthesis hop.
        for (int j = std::max(0, first-start); j < std::min(n, stop-start); ++j)
            noise[start+j-first] += bell[j]*block[j];
    }
    if (follow_fixed_envelope && edge > 0)
        follow_consonant_envelope(noise, residual, split, map, edge, sr);
    return noise;
}

struct ResidualMix {
    int edge, first, stop, join_start, join_end;
    double fixed, tail;
    bool follow_fixed_envelope;

    ResidualMix(const TimeMap& map, int sr, int count, const Options& o)
        : edge(std::min(count, static_cast<int>(std::ceil(map.front_o)))),
          first(o.fixed_residual_mix > 0. && edge > 0 ? 0 : edge), stop(edge),
          join_start(edge), join_end(edge),
          fixed(edge > 0 ? o.fixed_residual_mix/100. : o.residual_mix/100.),
          tail(edge < count ? o.residual_mix/100. : fixed),
          follow_fixed_envelope(o.fixed_residual_mix > 0.) {
        const double source_tail = map.ni-map.front_i;
        const bool short_join = tail > 0. && !map.identity() && source_tail <= stretch_hop(sr);
        if (edge > 0 && edge < count) {
            if (short_join && tail > fixed) {
                // Near a cutoff, the last fixed cycles can still leak into
                // the residual. Hand off before those cycles reach the join.
                // In particular a fully synthetic tail never reads WSOLA
                // from either side of the cutoff after the OTO boundary.
                const int width = static_cast<int>(std::lround(sr*.016));
                join_start = edge-std::min(width, std::max(1, edge/2));
                first = std::min(first, join_start);
            } else if (fixed != tail || (short_join && fixed < 1.)) {
                join_end += std::min(count-edge, static_cast<int>(std::lround(sr*.006)));
            }
        }
        stop = tail > 0. ? count : join_end;
    }
    bool enabled() const { return first < stop; }
    double transition(int n) const {
        if (n < join_start) return 0.;
        if (n >= join_end) return 1.;
        return smooth(double(n-join_start)/(join_end-join_start));
    }
};

// The unwarped source spectrum sets the noise level. Keep the old WSOLA
// branch from overpowering it, including a nonzero cutoff plateau, but
// never use a repeated fragment as gain reference for the generated noise.
void mix_residual(std::vector<float>& noise, const Warp& wave, const ResidualMix& edit, int sr) {
    const auto recorded = [&](int at) {
        const float h = wave.harmonic.empty() ? 0.f : wave.harmonic[at];
        const float removed = wave.residual_correction.empty() ? 0.f : wave.residual_correction[at];
        return (wave.dry[at]-h)-removed;
    };
    const int tail_count = edit.stop-edit.edge;
    const int span = std::max(1, static_cast<int>(std::lround(.046*sr)));
    const int frames = edit.tail > 0. && edit.tail < 1. ? (tail_count+span-1)/span : 0;
    const auto end = [&](int i) { return std::min(tail_count, (i+1)*span); };
    const auto center = [&](int i) { return .5*(i*span+end(i)-1); };
    std::vector<double> gains(frames);
    for (int i = 0; i < frames; ++i) {
        double old_energy = 0., new_energy = 0.;
        for (int j = i*span; j < end(i); ++j) {
            const int at = edit.edge+j;
            const double d = recorded(at), r = noise[at-edit.first];
            old_energy += d*d;
            new_energy += r*r;
        }
        // Half a dB of headroom avoids chasing random short-window variance.
        gains[i] = std::sqrt(unit(1.122*new_energy/std::max(1.e-24, old_energy)));
    }
    const double trim = smooth(2.*edit.tail);
    const double fixed_dry = std::sqrt(1.-edit.fixed), fixed_noise = std::sqrt(edit.fixed);
    const double tail_dry = std::sqrt(1.-edit.tail), tail_noise = std::sqrt(edit.tail);
    const double fixed_angle = std::asin(fixed_noise), tail_angle = std::asin(tail_noise);
    int frame = 0;
    for (int n = edit.first; n < edit.stop; ++n) {
        const int j = n-edit.edge;
        float& generated = noise[n-edit.first];
        const double w = edit.transition(n);
        if (w == 0.) {
            if (edit.fixed < 1.) generated = static_cast<float>(fixed_dry*recorded(n) + fixed_noise*generated);
            continue;
        }
        if (w == 1. && edit.tail >= 1.) continue;
        double gain = 1.;
        if (frames > 0) {
            while (frame+1 < frames && j >= center(frame+1)) ++frame;
            gain = gains[frame];
            if (frame+1 < frames && j > center(frame)) {
                const double u = (j-center(frame))/(center(frame+1)-center(frame));
                gain = mix(gain, gains[frame+1], smooth(u));
            }
            gain = mix(1., gain, trim*w);
        }
        double sampled_weight = tail_dry, generated_weight = tail_noise;
        if (w < 1. && edit.fixed != edit.tail) {
            // Both equal-power weights have zero slope at the endpoints.
            const double angle = mix(fixed_angle, tail_angle, w);
            sampled_weight = std::cos(angle);
            generated_weight = std::sin(angle);
        }
        generated = static_cast<float>(sampled_weight*gain*recorded(n) + generated_weight*generated);
    }
}

// Hold a supported *interior* cycle model if the source cursor approaches a
// hard cutoff. Do not move the target F0 lookup, the OTO cursor or a UV exit.
// With less than one period to the last mark, either pulse window would
// otherwise read zero padding for part of the cycle throughout the sustain.
class SustainSource {
    const TimeMap& time_;
    double edge_ = 0., limit_ = 0.;
public:
    SustainSource(const TimeMap& time, const std::vector<double>& marks,
                  const detail::HarmonicSplit& split, int sr, bool guard) : time_(time) {
        if (!guard || marks.size() < 3) return;
        const double end = time.last_source();
        const double period = marks.back()-marks[marks.size()-2];
        const double limit = std::min(marks[marks.size()-2], end-period-.5*split.hop);
        const double probe = std::max(period, .006*sr);
        if (limit < time.source_begin || end-limit > 3.*probe ||
            split.activity(end) < .24 || split.activity(limit) < .26 ||
            split.activity(.5*(limit+end)) < .24) return;
        limit_ = limit;
        edge_ = .5*std::min(probe, end-limit);
    }

    double operator()(double n) const {
        const double at = time_(n);
        if (edge_ <= 0. || at <= limit_-edge_) return at;
        const double u = unit((at-limit_+edge_)/(2.*edge_));
        const double safe = limit_-edge_*(1.-u)*(1.-u);
        return mix(at, safe, smooth((n-time_.front_o+2.*edge_)/(2.*edge_)));
    }
};

struct Grain {
    double center, a, b, mix_ab, period_l, period_r;
};

Grain grain(double q, double source, const pitch::Island& is) {
    const auto& p = is.marks;
    auto end = std::upper_bound(p.begin(), p.end(), source);
    if (end == p.begin()) {
        const double period = p[1]-p[0], a = p[0]-period;
        return {q, a, p[0], smooth((source-a)/period), period, period};
    }
    if (end == p.end()) {
        const double period = p.back()-p[p.size()-2];
        return {q, p.back(), p.back()+period, smooth((source-p.back())/period), period, period};
    }
    const auto j = static_cast<std::size_t>(end - p.begin());
    const double width = p[j] - p[j-1];
    const double w = unit((source - p[j-1]) / width);
    const double prev = j >= 2 ? p[j-1] - p[j-2] : width;
    const double next = j + 1 < p.size() ? p[j+1] - p[j] : width;
    return {q, p[j-1], p[j], smooth(w), mix(prev, width, w), mix(width, next, w)};
}

float read_grain(const Sinc& s, const Grain& g, int source_begin, int source_end, double d) {
    const auto read = [&](double c) { return s.read(c + d, source_begin, source_end); };
    if (g.mix_ab <= 0) return read(g.a);
    if (g.mix_ab >= 1) return read(g.b);
    return static_cast<float>(mix(read(g.a), read(g.b), g.mix_ab));
}

class TargetPitch {
    const TimeMap& time_;
    const pitch::Track& track_;
    const Options& options_;
    double center_, fallback_, bend_step_, transpose_, mod_, maximum_;
public:
    TargetPitch(const TimeMap& time, const pitch::Island& is, const pitch::Track& tr,
                const Options& o, double center, int sr)
        : time_(time), track_(tr), options_(o), center_(center), fallback_(sr/(is.marks[1]-is.marks[0])),
          bend_step_(sr*60./(96.*o.tempo)), transpose_(std::exp2(o.pitch_st/12.)),
          mod_(o.modulation*.01), maximum_(std::min(2400., sr/3.)) {}

    double operator()(double n) const {
        const auto& o = options_;
        double src = pitch::sample(track_, time_(n)).hz;
        if (src <= 0.) src = fallback_;
        double result = o.target_hz > 0 ? o.target_hz : src;
        if (o.target_hz > 0 && mod_ != 0 && center_ > 0)
            result *= std::pow(std::clamp(src/center_, .25, 4.), mod_);
        if (!o.bends.empty()) {
            const double u = n/bend_step_;
            const auto k = static_cast<std::size_t>(u);
            if (k < o.bends.size()) {
                const int next = k+1 < o.bends.size() ? o.bends[k+1] : 0;
                result *= std::exp2(mix(o.bends[k], next, u-k)/1200.);
            }
        }
        return std::clamp(result*transpose_, 30., maximum_);
    }
};

std::vector<Grain> target_grains(int active, const TimeMap& time, const pitch::Island& is,
                                const pitch::Track& tr, const Options& o, double center, int sr,
                                bool unchanged_clock, detail::TargetClock& target) {
    std::vector<Grain> out;
    if (unchanged_clock) {
        const double first = is.marks.front() - (is.marks[1] - is.marks.front());
        out.push_back(grain(first, first, is));
        for (double p : is.marks) out.push_back(grain(p, p, is));
        const double last = is.marks.back() + (is.marks.back() - is.marks[is.marks.size()-2]);
        out.push_back(grain(last, last, is));
        return out;
    }
    const TargetPitch frequency(time, is, tr, o, center, sr);
    const int step = std::max(1, static_cast<int>(std::lround(sr * .001)));
    target.hop = step;
    constexpr double start = 0.;
    const auto seed = grain(start, time(start), is);
    // Local phase at the first output position; positive-time crossings are
    // found by integrating a linearly interpolated target F0 in Hz.
    const double anchor = std::abs(seed.a - time(start)) < std::abs(seed.b - time(start)) ? seed.a : seed.b;
    double phase = (time(start) - anchor) / seed.period_l;
    phase -= std::floor(phase);
    target.phase = phase;
    if (phase < 1e-9) out.push_back(seed);
    double f = frequency(start);
    target.f0.reserve(1+(active+step-1)/step);
    target.f0.push_back(f);
    if (phase >= 1e-9) out.push_back(grain(start - phase * sr / f, time(start), is));
    for (int n = 0; n < active; n += step) {
        const int len = std::min(step, active - n);
        const double next = frequency(n + len);
        target.f0.push_back(next);
        const double df = (next - f) / len;
        const double add = .5 * (f + next) * len / sr;
        for (double delta = 1 - phase; delta <= add; delta += 1.) {
            const double a = delta * sr;
            const double root = std::sqrt(std::max(0., f*f + 2*df*a));
            const double q = n + std::min(static_cast<double>(len), 2*a / (f + root));
            out.push_back(grain(q, time(q), is));
        }
        phase += add;
        phase -= std::floor(phase);
        f = next;
    }
    if (out.empty() || out.back().center < active)
        out.push_back(grain(active + (1. - phase) * sr / f, time(active), is));
    return out;
}

// A large upward edit can put two positive portions of the source pulse
// in the same output cycle, even when the source H had no DC. Remove that
// local cycle mean with continuous interpolation, never an IIR across UV.
void remove_cycle_drift(std::vector<float>& wet, const std::vector<Grain>& grains, double ratio,
                        const detail::HarmonicSplit& split, const TimeMap& time) {
    struct Mean { double at, value; };
    std::vector<Mean> mean;
    mean.reserve(grains.size()-1);
    bool needed = false;
    const double last = static_cast<double>(wet.size());
    for (std::size_t j = 0; j+1 < grains.size(); ++j) {
        const auto& a = grains[j];
        const auto& b = grains[j+1];
        const double period = b.center-a.center;
        const double gain = smooth((1. - 2.*ratio*period/(a.period_r+b.period_l)-.13)/.15);
        const double start = std::max(0., a.center), end = std::min(last, b.center);
        if (end <= start || period <= 0.) continue;
        double integral = 0.;
        if (gain > 0.) {
            for (int n = static_cast<int>(start); n < static_cast<int>(std::ceil(end)); ++n)
                integral += wet[n] * (std::min(n+1., end) - std::max(double(n), start));
            needed = true;
        }
        mean.push_back({.5*(start+end), gain*integral/std::max(1., end-start)});
    }
    if (!needed) return;
    std::size_t j = 0;
    for (std::size_t n = 0; n < wet.size(); ++n) {
        while (j+1 < mean.size() && n >= mean[j+1].at) ++j;
        const auto& a = mean[j];
        const double dc = j+1 == mean.size() ? a.value :
            mix(a.value, mean[j+1].value, smooth((n-a.at)/(mean[j+1].at-a.at)));
        if (dc != 0.) wet[n] -= static_cast<float>(dc * smooth(split.activity(time(static_cast<double>(n))) / .20));
    }
}

void pad_tail(std::vector<float>& y, int count, int sr) {
    const int active = static_cast<int>(y.size());
    if (active >= count) return;
    const int edge = std::min(active, static_cast<int>(std::lround(.003 * sr)));
    for (int j = 0; j < edge; ++j) y[active - 1 - j] *= static_cast<float>(smooth(j / static_cast<double>(edge)));
    y.resize(count);
}

pitch::Track output_track(const pitch::Track& source, const TimeMap& time, const std::vector<Grain>& grains,
                          int sr, int length) {
    pitch::Track track;
    track.step = std::max(1, static_cast<int>(std::lround(.005*sr)));
    const int count = 1 + (length + track.step-1)/track.step;
    track.frames.reserve(count);
    std::size_t i = 0;
    for (int j = 0; j < count; ++j) {
        const double n = std::min(j*track.step, length);
        if (grains.size() >= 2) {
            while (i+2 < grains.size() && n >= grains[i+1].center) ++i;
            track.frames.push_back({sr/(grains[i+1].center-grains[i].center), 1.});
        } else {
            track.frames.push_back(pitch::sample(source, time(n)));
        }
    }
    return track;
}

void gain_stage(std::vector<float>& y, const Options& o, double reference_peak = -1.) {
    double gain = o.volume / 100.;
    if (o.peak_normalize) {
        double peak = reference_peak;
        if (peak < 0.) {
            peak = 0.;
            for (float v : y) peak = std::max(peak, double(std::abs(v)));
        }
        if (peak > 0) gain *= .5 / peak;
    }
    if (gain != 1) for (float& v : y) v = static_cast<float>(v * gain);
}

void validate(const std::vector<float>& input, int sr, const Options& o) {
    if (sr < 8000 || sr > 192000 || input.empty() || input.size() > max_samples)
        throw std::invalid_argument("expected 1..16000000 samples at 8..192 kHz");
    for (double v : {o.offset_ms, o.cutoff_ms, o.length_ms, o.consonant_ms, o.velocity,
                     o.target_hz, o.pitch_st, o.modulation, o.tempo, o.formant_st, o.forced_hz,
                     o.consonant_boost, o.opening, o.dynamics, o.volume, o.mask_return, o.residual_mix, o.fixed_residual_mix,
                     o.tension, o.breathiness, o.mouth, o.resonance, o.formant_emphasis})
        if (!std::isfinite(v)) throw std::invalid_argument("non-finite option");
    if (o.offset_ms < 0 || o.consonant_ms < 0 || (o.length_ms < 0 && o.length_ms != -1) ||
        o.velocity < 0 || o.velocity > 400 ||
        o.tempo < 1 || o.tempo > 1000 || o.target_hz < 0 || o.target_hz > 2400 ||
        (o.target_hz > 0 && o.target_hz < 30) || std::abs(o.pitch_st) > 48 ||
        std::abs(o.formant_st) > 24 || o.modulation < -200 || o.modulation > 200 ||
        (o.forced_hz != -1 && o.forced_hz != 0 && (o.forced_hz < 30 || o.forced_hz > 2400)) ||
        o.consonant_boost < 0 || o.consonant_boost > 100 ||
        std::abs(o.opening) > 100 || o.dynamics < 0 || o.dynamics > 100 ||
        o.volume < 0 || o.volume > 1000 || o.bends.size() > 2'000'000 || (!o.f0.empty() && o.forced_hz != 0))
        throw std::invalid_argument("option outside supported range");
    if (o.mask_return < 0. || o.mask_return > 1.)
        throw std::invalid_argument("mask return must be 0..1");
    if (o.residual_mix < 0. || o.residual_mix > 100. || o.fixed_residual_mix < 0. || o.fixed_residual_mix > 100.)
        throw std::invalid_argument("residual mixes must be 0..100");
    for (double v : {o.tension, o.breathiness, o.mouth, o.resonance, o.formant_emphasis})
        if (std::abs(v) > 100.) throw std::invalid_argument("voice effect must be -100..100");
    for (float v : input) if (!std::isfinite(v) || std::abs(v) > 64.f)
        throw std::invalid_argument("expected finite audio with magnitude <= 64");
    for (int v : o.bends) if (v < -2048 || v > 2047) throw std::invalid_argument("invalid pitch bend");
    if (!o.mrq_f0.frames.empty()) {
        if (!o.f0.empty() || o.forced_hz != 0. || o.mrq_f0.step < 1 || o.mrq_f0.step > sr ||
            o.mrq_f0.frames.size() > 1'000'000)
            throw std::invalid_argument("invalid MRQ F0 track or conflicting source F0 options");
        for (const auto f : o.mrq_f0.frames)
            if (!(f.hz == 0. || (f.hz >= 30. && f.hz <= std::min(2400., sr / 3.))) ||
                f.confidence != (f.hz == 0. ? 0. : 1.))
                throw std::invalid_argument("invalid MRQ F0 frame");
    }
    double last = -1;
    for (const auto& p : o.f0) {
        if (!std::isfinite(p.time_s) || p.time_s < 0 || p.time_s <= last || !std::isfinite(p.hz) ||
            (p.hz != 0 && (p.hz < 30 || p.hz > std::min(2400., sr / 3.))))
            throw std::invalid_argument("invalid F0: use increasing times and 0 or 30..2400 Hz");
        last = p.time_s;
    }
}
} // namespace

static std::vector<float> render_impl(const std::vector<float>& input, int sr, const Options& o,
                                      Report* report, NoteModel* model) {
    validate(input, sr, o);
    Report result;
    result.input_samples = input.size();
#ifdef TN_FINDS_NEBULA_F0
    result.f0_backend = "nebula";
#endif
    const double ms = sr / 1000.;
    const double ni = static_cast<double>(input.size());
    const double start = std::round(std::min(ni, o.offset_ms * ms));
    const double end = std::round(std::clamp(o.cutoff_ms < 0 ? start - o.cutoff_ms * ms
                                                           : ni - o.cutoff_ms * ms, 0., ni));
    if (start >= end) throw std::invalid_argument("offset and cutoff select no input samples");
    const double nout = o.length_ms == -1 ? end - start : std::round(o.length_ms * ms);
    if (nout < 1 || nout > max_samples) throw std::invalid_argument("output length is outside 1..16000000 samples");
    const int count = static_cast<int>(nout);
    const int first = static_cast<int>(start), last = static_cast<int>(end);
    result.output_samples = count;
    const double front = std::min(o.consonant_ms * ms, end - start);
    const double fixed_out = front * std::exp2(1. - o.velocity / 100.);
    TimeMap time{end - start, nout, front, std::min(nout, fixed_out), fixed_out};
    const bool bend = std::any_of(o.bends.begin(), o.bends.end(), [](int v) { return v != 0; });
    const bool untouched_pitch = o.target_hz == 0 && o.pitch_st == 0 && !bend;
    const int hop = stretch_hop(sr);
    time.hold_end = !untouched_pitch && time.ni-front <= 1. && nout - time.front_o >= 2.*hop;
    const bool sustain_guard = time.needs_sustain_guard(hop);
    const int active = time.active_samples();
    const ResidualMix residual_mix(time, sr, active, o);
    const bool replace_residual = residual_mix.enabled();
    // A cut can leave fewer than two periods at the offset. Use preceding
    // context for analysis, but never let the next syllable supply the tail.
    const int context = sustain_guard ? static_cast<int>(std::lround(.050*sr)) : 0;
    const int analysis_first = std::max(0, first-context);
    time.source_begin = first-analysis_first;
    std::vector<float> x(input.begin()+analysis_first, input.begin()+last);
    const bool color_residual = o.mouth != 0. || o.resonance != 0. || o.formant_emphasis != 0.;
    const bool color_harmonic = color_residual || o.tension != 0.;
    const bool no_voice_effects = o.formant_st == 0 && o.opening == 0 && o.dynamics == 0 &&
                                  !color_harmonic && o.breathiness == 0.;
    const bool unchanged_clock = !model && time.identity() && untouched_pitch;
    if (unchanged_clock && no_voice_effects && o.consonant_boost == 0 && !replace_residual) {
        const auto ts = Clock::now();
        pad_tail(x, count, sr);
        gain_stage(x, o);
        result.bypass = true;
        result.synthesis_ms = milliseconds(ts);
        if (report) *report = result;
        return x;
    }

    auto ts = Clock::now();
    const bool spectral_formant = !o.grain_formant && o.formant_st != 0;
    const double formant_ratio = std::exp2(o.formant_st / 12.);
    const double ratio = spectral_formant ? 1. : formant_ratio;
    const bool resynthesize = !unchanged_clock || ratio != 1.;
    const bool source_filter = resynthesize && ratio == 1. && !o.phase_lock;
    if (!o.mrq_f0.frames.empty()) result.f0_backend = "mrq";
    else if (!o.f0.empty()) result.f0_backend = "csv";
    else if (o.forced_hz != 0.) result.f0_backend = "forced";
    const detail::SourceAnalysisInput analysis_input{input, x, o.f0, o.mrq_f0, sr, analysis_first, o.forced_hz,
                                                     o.mask_return, resynthesize, o.tension != 0. || source_filter};
#ifdef TN_FINDS_ANALYSIS_CACHE_ID
    const auto analysis = detail::cached_source_analysis(analysis_input, analyze_source, result.analysis_cache_hit);
#else
    const auto analysis = analyze_source(analysis_input);
#endif
    const auto& clock = analysis.clock;
    const auto& split = analysis.split;
    const auto& source_rd = analysis.rd;
    const auto& h = split.harmonic;
    if (!h.empty()) {
        // Bridging supplies a clock, even through breaths. Those interpolated
        // pitches must not contribute to the representative source note.
        auto measured = clock;
        for (std::size_t i = 0; i < measured.frames.size(); ++i) {
            auto& f = measured.frames[i];
            const auto pos = static_cast<double>(i * measured.step);
            if (split.activity(pos) < .22) f.hz = 0.;
            f.confidence *= split.activity(pos);
        }
        result.source_pitch_hz = center_pitch(measured, time);
        if ((o.loop || o.force_loop) && !time.identity())
            choose_loop(time, analysis.pulses.marks, split, sr, o.force_loop);
    }
    result.islands = split.regions();
    result.loop_used = time.looping();
    result.analysis_ms = milliseconds(ts);
    ts = Clock::now();
    if (model) {
        model->sample_rate = sr;
        model->active = active;
    }
    std::vector<float> correction;
    if (spectral_formant && !h.empty()) {
        // Shift the residual after low-band leakage cleanup. Carry the edit
        // through the same WSOLA windows as the unedited recording, while
        // keeping the alignment search on the original residual.
        const auto& cleanup = split.cleanup;
        const auto original = [&](std::size_t i) {
            return (x[i] - h[i]) - (cleanup.empty() ? 0.f : cleanup[i]);
        };
        correction.resize(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) correction[i] = original(i);
        constexpr double breath_smoothing_hz = 220.;
        detail::shift_formant_envelope(correction, sr, formant_ratio, clock, breath_smoothing_hz);
        for (std::size_t i = 0; i < x.size(); ++i) {
            const float voiced = static_cast<float>(smooth(split.activity(static_cast<double>(i)) / .30));
            const float removed = cleanup.empty() ? 0.f : cleanup[i];
            correction[i] = removed - voiced * (correction[i] - original(i));
        }
    }
    const auto& dry_correction = correction.empty() ? split.cleanup : correction;
    auto wave = stretch(x, split, dry_correction, time, sr, active);
    std::vector<float> mixed;
    if (replace_residual) {
        mixed = synthesize_residual(x, split, dry_correction, clock, time, sr, residual_mix.first,
                                    residual_mix.edge, residual_mix.stop, active,
                                    residual_mix.follow_fixed_envelope);
        mix_residual(mixed, wave, residual_mix, sr);
    }
    // Anchor the existing normalization to Tns0 so a different noise
    // crest factor never changes the gain on the harmonics or attack.
    const bool anchor_peak = replace_residual && o.peak_normalize;
    double recorded_peak = 0.;
    const auto measure = [&](float v) {
        if (anchor_peak) recorded_peak = std::max(recorded_peak, double(std::abs(v)));
    };
    auto& y = wave.dry;
    const double boost = 5. * std::pow(o.consonant_boost / 100., 2.);
    if (!h.empty()) {
        std::vector<float> shaped;
        if (spectral_formant) {
            shaped = h;
            detail::shift_formant_envelope(shaped, sr, formant_ratio, clock);
            // A formant-EQ window can overlap a fully dry consonant. Limit
            // just the EQ delta by the measured support (not a UV label).
            for (std::size_t i = 0; i < shaped.size(); ++i) {
                const double a = smooth(split.activity(static_cast<double>(i)) / .30);
                shaped[i] = h[i] + static_cast<float>(a) * (shaped[i] - h[i]);
            }
        }
        const auto& source = spectral_formant ? shaped : h;
        const auto& is = analysis.pulses;
        const SustainSource sustain_source(time, is.marks, split, sr, sustain_guard);
        detail::TargetClock target;
        const auto grains = target_grains(active, time, is, clock, o, result.source_pitch_hz, sr,
                                          unchanged_clock, target);
        std::vector<float> periodic_out;
        if (resynthesize && grains.size() >= 2) {
            if (!model) periodic_out.resize(active);
            if (ratio == 1.) {
                std::vector<detail::PitchEpoch> epochs;
                epochs.reserve(grains.size());
                for (std::size_t j = 0; j < grains.size(); ++j) {
                    const auto left = j == 0 ? j : j-1;
                    const auto right = std::min(j+1, grains.size()-1);
                    const double at = grains[j].center;
                    const double hz = sr*double(right-left)/(grains[right].center-grains[left].center);
                    epochs.push_back({at, sustain_source(std::clamp(at, 0., double(active))), hz});
                }
                const detail::PhaseReference reference{h, source_rd, formant_ratio};
                if (model) {
                    detail::prepare_harmonics(source, is.marks, epochs, split, sr, active,
                                              source_filter ? &reference : nullptr, *model);
                    model->pitch_hop = target.hop;
                    model->f0.reserve(target.f0.size());
                    for (double hz : target.f0)
                        model->f0.push_back(static_cast<float>(hz));
                } else {
                    detail::synthesize_harmonics(source, is.marks, epochs, split, sr, periodic_out, target,
                                                 source_filter ? &reference : nullptr);
                }
            } else {
                // Explicit --grain-formant mode warps each pulse in sample time.
                const Sinc sinc(source, ratio);
                const int source_end = static_cast<int>(source.size());
                for (std::size_t j = 0; j + 1 < grains.size(); ++j) {
                    const auto& l = grains[j];
                    const auto& r = grains[j+1];
                    const double period = r.center-l.center;
                    if (period <= 0.) continue;
                    const int lo = std::max(0, static_cast<int>(std::ceil(l.center)));
                    const int hi = std::min(active, static_cast<int>(std::ceil(r.center)));
                    const double extent_l = std::min(1., l.period_r / (ratio * period));
                    const double extent_r = std::min(1., r.period_l / (ratio * period));
                    for (int n = lo; n < hi; ++n) {
                        const double u = unit((n-l.center) / period), blend = smooth(u);
                        const double v = (1.-blend) * smooth((extent_l-u) / std::min(.16, .6*extent_l));
                        const double w = blend * smooth((extent_r-(1.-u)) / std::min(.16, .6*extent_r));
                        float wet = 0.f;
                        if (v > 0.) wet += static_cast<float>(v) * read_grain(sinc, l, 0, source_end, (n-l.center)*ratio);
                        if (w > 0.) wet += static_cast<float>(w) * read_grain(sinc, r, 0, source_end, (n-r.center)*ratio);
                        periodic_out[n] = wet / static_cast<float>(std::max(.32, v+w));
                    }
                }
            }
            if (o.phase_lock && grains.size() >= 5) {
                std::vector<double> dst, at;
                dst.reserve(grains.size()); at.reserve(grains.size());
                for (const auto& g : grains) {
                    dst.push_back(g.center);
                    at.push_back(sustain_source(g.center));
                }
                detail::lock_harmonic_phase(source, is.marks, dst, at, sr, 0, periodic_out);
            }
            // The spectral tables already have an empty DC bin.
            if (ratio != 1.) remove_cycle_drift(periodic_out, grains, ratio, split, time);
        }

        // The same WSOLA windows carry the recording, harmonic estimate,
        // and residual corrections. Replace the warped harmonic estimate
        // with the new periodic waveform and mix the corrected residual.
        const float shelf_k = static_cast<float>(1. - std::exp(-2 * pi * 1700. / sr));
        const double shelf_gain = std::exp2(o.opening / 125.) - 1.;
        const double mb = o.breathiness * .01;
        const double whisper = std::max(0., mb);
        const double breath_gain = 1. + mb;
        const bool edit_voice = color_harmonic || breath_gain != 1. || whisper > 0.;
        const double level_hz = o.target_hz > 0. ? o.target_hz * std::exp2(o.pitch_st / 12.) : result.source_pitch_hz;
        const bool effects = edit_voice || shelf_gain != 0. || boost > 0. || o.dynamics > 0.;
        const auto& wet_wave = !periodic_out.empty() ? periodic_out : spectral_formant ? source : wave.harmonic;
        pitch::Track edited_clock;
        if (!unchanged_clock && (o.tension != 0. || o.mouth != 0. || o.formant_emphasis != 0.))
            edited_clock = output_track(clock, time, grains, sr, active);
        const auto& effect_clock = unchanged_clock ? clock : edited_clock;
        detail::GlottalTrack mapped_rd;
        if (!model && !unchanged_clock && !source_rd.rd.empty()) {
            mapped_rd.hop = effect_clock.step;
            mapped_rd.rd.resize(effect_clock.frames.size());
            for (std::size_t i = 0; i < mapped_rd.rd.size(); ++i)
                mapped_rd.rd[i] = source_rd.sample(sustain_source(std::min(double(i)*mapped_rd.hop, double(active))));
        }
        std::vector<float> colored_wet, colored_dry, colored_mixed;
        if (model) {
            for (std::size_t j = 0; j < model->frames.size(); ++j) {
                const double pos = sustain_source(std::min(double(active), double(j*model->hop)));
                detail::color_harmonics(model->frames[j], sr, o, source_rd.sample(pos),
                                        smooth(split.activity(pos)/.30), level_hz);
            }
        } else if (color_harmonic) {
            colored_wet = wet_wave;
            detail::shape_voice_envelope(colored_wet, sr, effect_clock, o, true,
                                         unchanged_clock ? &source_rd : &mapped_rd);
        }
        if (color_residual) {
            colored_dry.resize(active);
            for (int n = 0; n < active; ++n) {
                const float removed = wave.residual_correction.empty() ? 0.f : wave.residual_correction[n];
                colored_dry[n] = (y[n]-wave.harmonic[n])-removed;
            }
            if (replace_residual) {
                colored_mixed = colored_dry;
                std::copy(mixed.begin(), mixed.end(), colored_mixed.begin()+residual_mix.first);
                detail::shape_voice_envelope(colored_mixed, sr, effect_clock, o, false);
            }
            detail::shape_voice_envelope(colored_dry, sr, effect_clock, o, false);
        }
        const int color_fade = std::min(static_cast<int>(std::lround(sr*.006)),
                                        residual_mix.join_end-residual_mix.join_start);
        float lp = 0;
        std::size_t g = 0;
        for (int n = 0; n < active; ++n) {
            const float base = wave.harmonic[n];
            const float residue = y[n]-base;
            const float removed = wave.residual_correction.empty() ? 0.f : wave.residual_correction[n];
            const bool use_mixed = n >= residual_mix.first && n < residual_mix.stop;
            const float alternate = use_mixed ? mixed[n-residual_mix.first] : 0.f;
            const auto color = [&](float& dry, float& edited, float amount) {
                const float dry_delta = colored_dry[n]-dry;
                dry += amount*dry_delta;
                if (!use_mixed) return;
                // Spectral EQ windows can see the splice before its first
                // sample. Bring their correction in and out continuously too.
                double weight = 1.;
                if (color_fade > 0) {
                    if (residual_mix.first > 0)
                        weight *= smooth(double(n-residual_mix.first)/color_fade);
                    if (residual_mix.stop < active)
                        weight *= smooth(double(residual_mix.stop-n)/color_fade);
                }
                edited += amount*static_cast<float>(mix(dry_delta, colored_mixed[n]-edited, weight));
            };
            if (model) {
                float dry = residue-removed, edited = alternate;
                if (color_residual || breath_gain != 1. || boost > 0.) {
                    const double pos = time(n);
                    const double voice_mix = smooth(split.activity(pos)/.30);
                    if (color_residual) color(dry, edited, static_cast<float>(voice_mix));
                    const float gain = static_cast<float>((1.+(breath_gain-1.)*voice_mix)*(1.+(1.-split.share(pos))*boost));
                    dry *= gain;
                    edited *= gain;
                }
                measure(dry);
                y[n] = use_mixed ? edited : dry;
                continue;
            }
            const float wet = wet_wave[n];
            if (!effects) {
                const float reference = y[n] + ((wet-base)-removed);
                measure(reference);
                y[n] = use_mixed ? wet + alternate : reference;
                continue;
            }
            const double pos = time(n);
            const float activity = static_cast<float>(split.activity(pos));
            float dry = residue - removed, edited = alternate;
            float voiced = wet;
            const double voice_mix = edit_voice ? smooth(activity / .30) : 0.;
            if (color_harmonic && voice_mix > 0.) {
                const float gain = static_cast<float>(voice_mix);
                voiced += gain*(colored_wet[n]-voiced);
                if (color_residual) color(dry, edited, gain);
            }
            if (shelf_gain != 0.) {
                lp += shelf_k * (wet-lp);
                voiced += static_cast<float>(shelf_gain) * (wet-lp);
            }
            // Mb scales the voiced part after the residual source is chosen.
            if (breath_gain != 1.) {
                const float gain = static_cast<float>(1. + (breath_gain - 1.) * voice_mix);
                dry *= gain;
                edited *= gain;
            }
            if (o.dynamics != 0. && level_hz > 0. && grains.size() >= 2) {
                while (g+2 < grains.size() && n >= grains[g+1].center) ++g;
                const double hz = sr/(grains[g+1].center-grains[g].center);
                voiced *= static_cast<float>(1. + o.dynamics * .0012 *
                    std::clamp(std::log2(hz / level_hz) * 4., -1., 1.));
            }
            if (whisper > 0.) voiced *= static_cast<float>(1. - whisper * voice_mix);
            const float gain = static_cast<float>(1. + (1. - split.share(pos)) * boost);
            measure(voiced + dry*gain);
            y[n] = voiced + (use_mixed ? edited : dry)*gain;
        }
    } else if (replace_residual) {
        const float gain = static_cast<float>(1. + boost);
        for (int n = 0; n < active; ++n) {
            measure(y[n]*gain);
            const bool use_mixed = n >= residual_mix.first && n < residual_mix.stop;
            y[n] = (use_mixed ? mixed[n-residual_mix.first] : y[n])*gain;
        }
    } else if (boost > 0.) {
        for (float& v : y) v *= static_cast<float>(1. + boost);
    }
    pad_tail(y, count, sr);
    if (model) {
        double gain = o.volume*.01;
        if (o.peak_normalize) {
            double dry_peak = recorded_peak, harmonic_peak = 0.;
            if (!anchor_peak)
                for (float v : y) dry_peak = std::max(dry_peak, double(std::abs(v)));
            for (const auto& frame : model->frames) {
                double bound = 0.;
                for (const auto coefficient : frame.harmonics) bound += 2.*std::abs(coefficient);
                harmonic_peak = std::max(harmonic_peak, bound);
            }
            const double peak = std::hypot(dry_peak, harmonic_peak);
            if (peak > 0.) gain *= .5/peak;
        }
        if (gain != 1.) {
            for (float& v : y) v = static_cast<float>(v*gain);
            for (auto& frame : model->frames)
                for (auto& coefficient : frame.harmonics) coefficient *= static_cast<float>(gain);
        }
        model->residual = std::move(y);
        return {};
    }
    gain_stage(y, o, anchor_peak ? recorded_peak : -1.);
    result.synthesis_ms = milliseconds(ts);
    if (report) *report = result;
    return std::move(y);
}
std::vector<float> render(const std::vector<float>& input, int sr, const Options& o, Report* report) {
    return render_impl(input, sr, o, report, nullptr);
}

NoteModel prepare_note(const std::vector<float>& input, int sr, const Options& o) {
    if (o.target_hz <= 0. || o.grain_formant || o.phase_lock)
        throw std::invalid_argument("note preparation requires classic harmonic synthesis");
    NoteModel model;
    render_impl(input, sr, o, nullptr, &model);
    return model;
}
} // namespace tn
