// SPDX-License-Identifier: GPL-3.0-only
#include "harmonic_mask.hpp"
#include "harmonic_ridges.hpp"
#include "dsp_fft.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace tn::detail {
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr int bands = 12;

float ease(double x) {
    x = std::clamp(x, 0., 1.);
    return static_cast<float>(x*x*(3. - 2.*x));
}

double sample_curve(const std::vector<float>& x, int hop, double sample) {
    if (x.empty()) return 0.;
    const double p = std::clamp(sample / hop, 0., double(x.size() - 1));
    const auto k = static_cast<std::size_t>(p);
    const auto w = ease(p - k);
    return x[k] + w * (x[std::min(k+1, x.size()-1)] - x[k]);
}

float cubic(const std::vector<float>& x, double p) {
    const int n = static_cast<int>(x.size());
    if (p < 0. || p > n - 1.) return 0.f;
    const int k = static_cast<int>(p);
    const float t = static_cast<float>(p-k), b = x[k];
    if (k == n - 1) return b;
    const float c = x[k+1];
    if (k == 0 || k + 2 >= n) return b + t * (c-b);
    const float a = x[k-1], d = x[k+2];
    return b + .5f * t * (c-a + t*(2*a-5*b+4*c-d + t*(3*(b-c)+d-a)));
}

// Neighbor directions are two complex, phase-aligned predictions of the
// same source-time window. Do not derive the mask from the magnitude of the
// input alone, or from the periodic average (whose center contains noise).
// Pack x+i*x_left and x_right in two FFTs; only the averaged periodic
// prediction goes through the inverse transform.
void conjugate_pair(const std::vector<Complex>& from, int k, int n, Complex& re, Complex& im) {
    const Complex a = from[k], b = std::conj(from[n-k]);
    re = .5f * (a + b);
    const Complex difference = a-b;
    im = Complex(.5f*difference.imag(), -.5f*difference.real());
}

void pack_pair(std::vector<Complex>& to, int k, int n, Complex a, Complex b) {
    to[k] = Complex(a.real()-b.imag(), a.imag()+b.real());
    to[n-k] = Complex(a.real()+b.imag(), b.real()-a.imag());
}

struct Point {
    double left = 0., right = 0.;
    float wl = 0.f, wr = 0.f;
};

// A single phase mapping is shared by every frequency range. Crossfade the
// availability of a neighbor at the real file edges, never at a UV label.
class Alignment {
    std::vector<std::array<double, 3>> spans_;
    int len_;
    int i_ = 1;
public:
    Alignment(const std::vector<double>& p, int len) : len_(len) {
        spans_.reserve(p.size()+4);
        const int last = static_cast<int>(p.size())-1;
        const double l = p[1]-p[0], r = p.back()-p[p.size()-2];
        for (int i = -2; i <= last+2; ++i) {
            const double a = i < 0 ? p.front()+i*l : i > last ? p.back()+(i-last)*r : p[i];
            const double span = i < 0 ? l : i >= last ? r : p[i+1]-a;
            spans_.push_back({a, span, 1./span});
        }
    }
    Point at(int n) {
        while (i_+2 < static_cast<int>(spans_.size()) && n >= spans_[i_+1][0]) ++i_;
        const auto& a = spans_[i_];
        const double u = (n-a[0]) * a[2];
        const auto& l = spans_[i_-1];
        const auto& r = spans_[i_+1];
        const double left = l[0]+u*l[1], right = r[0]+u*r[1];
        // Half a period of continuous support for a missing neighbor. At a
        // file end, the other direction can still provide acoustic evidence.
        const auto support = [&](double pos) {
            return ease(2. * std::min(pos, len_ - 1. - pos) * a[2]);
        };
        return {left, right, support(left), support(right)};
    }
};

struct BinEvidence {
    double power, left_power, right_power, left_cross, right_cross;
    double agreement, control, pair_cross;

    static float dot(Complex a, Complex b) {
        return a.real()*b.real() + a.imag()*b.imag();
    }

    BinEvidence(Complex a, Complex l, Complex r, float wl, float wr, bool paired)
        : power(std::norm(a)), left_power(std::norm(l)), right_power(std::norm(r)),
          left_cross(dot(a,l)), right_cross(dot(a,r)), agreement(left_cross*wl + right_cross*wr),
          control(wl*wl+wr*wr), pair_cross(paired ? dot(l,r) : 0.f) {}
};

struct Evidence {
    double power = 0., neighbor = 0., agree = 0.;
    double control = 0., paired = 0., agreement_lr = 0.;
    double left_power = 0., right_power = 0., left_cross = 0., right_cross = 0.;
    void add(double weight, const BinEvidence& bin, bool both) {
        const double sum = bin.left_power + bin.right_power;
        power += weight * bin.power;
        neighbor += weight * sum;
        agree += weight * bin.agreement;
        left_power += weight * bin.left_power;
        right_power += weight * bin.right_power;
        left_cross += weight * bin.left_cross;
        right_cross += weight * bin.right_cross;
        control += weight * bin.control * bin.power;
        if (both) {
            paired += weight * sum;
            agreement_lr += 2. * weight * bin.pair_cross;
        }
    }
    float mask(int band, double floor) const {
        if (power <= floor || neighbor <= floor || control <= floor) return 0.f;
        // Repeatability includes both normalized correlation and amplitude
        // agreement. At a gradual onset/release retain the available side;
        // do not let a very different second neighbor veto a true tail.
        const double repeat = 2. * agree / (control + neighbor);
        const double both = paired > floor ? agreement_lr / paired : repeat;
        const double verified = .64 * repeat + .36 * std::min(repeat, both + .13);
        // The shared low-order check needs stronger evidence than an
        // individual band, especially when the available clip is short.
        const double low = band < 0 ? .78 : band < 3 ? .67 : band < 5 ? .68 : .62;
        const double high = band < 0 ? .94 : band < 3 ? .92 : band < 5 ? .93 : .95;
        return ease((verified-low)/(high-low)) * ease((power-floor)/(9.*floor));
    }

    // During attacks and releases the amplitude-sensitive score can reject
    // a repeatable waveform. Check the two directions after normalizing
    // their levels. Both must contain measurable, aligned sound.
    float transition(double floor, float wl, float wr) const {
        if (std::min(wl, wr) < .8f || std::min({power, left_power, right_power}) <= floor)
            return 0.f;
        const double l = std::sqrt(left_power/power)/wl;
        const double r = std::sqrt(right_power/power)/wr;
        const double rho_l = left_cross/std::sqrt(power*left_power);
        const double rho_r = right_cross/std::sqrt(power*right_power);
        const double motion = std::max(std::abs(l-1.), std::abs(r-1.));
        // A transition may reach only one side of this window. The closer
        // neighbor should then be an especially good match to the center.
        const double match = std::max(rho_l, rho_r);
        const double agree_both = std::min(rho_l, rho_r);
        const double weaker_level = std::min(l,r);
        return ease((motion-.16)/.36) * ease((match-.78)/.17) *
               ease((agree_both-.56)/.24) * ease((weaker_level-.09)/.14) *
               ease((power-floor)/(9.*floor));
    }

    // A neighboring voiced region can lend evidence in one direction at
    // a boundary. Require a good local match and a genuinely weaker or
    // less coherent opposite side before enabling that fallback.
    float one_sided(double floor, float wl, float wr, bool left) const {
        const float available = left ? wl : wr;
        const double ref_power = left ? left_power : right_power;
        const double other_power = left ? right_power : left_power;
        const double ref_cross = left ? left_cross : right_cross;
        const double other_cross = left ? right_cross : left_cross;
        if (available < .8f || std::min(power, ref_power) <= floor) return 0.f;
        const double rho = ref_cross/std::sqrt(power*ref_power);
        const double other = other_power > floor ? other_cross/std::sqrt(power*other_power) : 0.;
        const double level = std::sqrt(ref_power/power)/available;
        const double weak = std::sqrt(other_power/ref_power);
        const float asymmetry = std::max(1.f-ease((weak-.28)/.53), ease((rho-other-.08)/.23));
        return ease((rho-.65)/.26) * ease((level-.74)/.42) * asymmetry *
               ease((power-floor)/(9.*floor));
    }

    float neighbor_gain(bool left) const {
        const double energy = left ? left_power : right_power;
        const double cross = left ? left_cross : right_cross;
        return energy > 0. ? static_cast<float>(std::clamp(cross/energy, 0., 3.)) : 0.f;
    }

    std::array<float, 2> neighbor_gains() const {
        // Give a phase-matched neighbor more weight and exclude negative
        // phase matches from the normalization.
        const double l = std::max(0., left_cross), r = std::max(0., right_cross);
        const double total = std::max(1e-30, l*l + r*r);
        return {static_cast<float>(l*l/total)*neighbor_gain(true),
                static_cast<float>(r*r/total)*neighbor_gain(false)};
    }
};

// Keep the pre-mask average and its one-cycle fades at voiced island edges.
std::vector<float> cycle_average(const std::vector<float>& x, const std::vector<pitch::Island>& islands) {
    constexpr std::array<float, 5> weights{.12f, .22f, .32f, .22f, .12f};
    std::vector<float> average(x.size());
    for (const auto& island : islands) {
        const auto& marks = island.marks;
        const int end = static_cast<int>(marks.size()) - 1;
        for (int cycle = 0; cycle < end; ++cycle) {
            const int first = std::max(0, cycle-2), last = std::min(end-1, cycle+2);
            const int begin = static_cast<int>(std::ceil(marks[cycle]));
            const int stop = std::min(static_cast<int>(x.size()), static_cast<int>(std::ceil(marks[cycle+1])));
            const double inv = 1./(marks[cycle+1]-marks[cycle]);
            const bool interior = first > 0 && last < end-1;
            for (int i = begin; i < stop; ++i) {
                const double phase = (i-marks[cycle])*inv;
                const float up = interior ? 1.f : ease(phase);
                const float down = interior ? 1.f : ease(1.-phase);
                const float edge = (cycle == 0 ? up : 1.f)*(cycle == end-1 ? down : 1.f);
                float sum = 0.f;
                double mass = 0.;
                for (int j = first; j <= last; ++j) {
                    const float w = weights[j-cycle+2]*(j == 0 ? up : 1.f)*(j == end-1 ? down : 1.f);
                    const float sample = j == cycle ? x[i] : cubic(x, marks[j]+phase*(marks[j+1]-marks[j]));
                    sum += w*sample;
                    mass += w;
                }
                if (mass > 1e-9) average[i] = static_cast<float>(sum/mass*edge);
            }
        }
    }
    return average;
}

float upper_band(double hz) { return ease((hz-1350.)/750.); }

} // namespace

double HarmonicSplit::activity(double at) const { return sample_curve(strength, hop, at); }
double HarmonicSplit::share(double at) const { return sample_curve(periodic_share, hop, at); }

std::size_t HarmonicSplit::regions() const {
    std::size_t count = 0;
    bool prior = false;
    for (const auto x : strength) {
        const bool active = x > .025f;
        count += active && !prior;
        prior = active;
    }
    return count;
}

HarmonicSplit split_harmonic(const std::vector<float>& x, int sr,
                             const std::vector<double>& marks, const pitch::Track& track,
                             bool cleanup, double amount, const pitch::Track* automatic_track) {
    HarmonicSplit out;
    if (marks.size() < 2 || x.size() < 8) return out;
    // The trial clock spans the whole clip, including frames with no reliable
    // pitch estimate. 32 ms windows use overlapping short-time evidence,
    // while the aligned waveform itself retains the low note's full period.
    const int n = std::clamp(1 << static_cast<int>(std::lround(std::log2(sr * .032))), 256, 8192);
    const int half = n/2, hop = n/4, len = static_cast<int>(x.size());
    out.hop = hop;
    const auto count = static_cast<std::size_t>((len-1)/hop+2);
    out.strength.resize(count);
    out.periodic_share.resize(count);
    out.harmonic.resize(x.size());
    if (cleanup) out.cleanup.resize(x.size());
    const Fft fft(n);
    std::vector<Complex> pair(n), next(n), edge_spectrum(n);
    std::vector<float> win(n);
    for (int j = 0; j < n; ++j) win[j] = static_cast<float>(std::sin(pi*j/n));
    // Closer spacing around the harmonics-to-frication crossover keeps the
    // top of a vowel from extending the mask through the adjacent air band.
    constexpr std::array<double, bands> centers{90., 275., 600., 1175., 2050., 2900., 3700., 4600., 5300., 6600., 9000., 15000.};
    struct Weight { int a, b; float w, dc, low, upper; };
    std::vector<Weight> weights(half);
    for (int k = 1; k < half; ++k) {
        const double hz = k * (double(sr)/n);
        const auto it = std::upper_bound(centers.begin(), centers.end(), hz);
        const int b = std::min(bands-1, static_cast<int>(it-centers.begin())), a = std::max(0, b-1);
        const float v = b == a ? 0.f : ease(std::log(hz/centers[a])/std::log(centers[b]/centers[a]));
        weights[k] = {a, b, v, ease((hz-12.)/18.), 1.f-ease((hz-650.)/1450.), upper_band(hz)};
    }
    const double peak = std::abs(*std::max_element(x.begin(), x.end(),
        [](float a, float b) { return std::abs(a) < std::abs(b); }));
    if (peak < 1e-8) { out.harmonic.clear(); out.cleanup.clear(); return out; }
    const auto measured = automatic_track ? pitch::Track{} : pitch::estimate(x, sr);
    const auto& seed = automatic_track ? *automatic_track : measured;
    auto voiced_track = seed;
    if (!automatic_track) {
        // Supplied F0 controls the phase, while measured voicing prevents it
        // from turning an isolated resonance or a breath into a voiced island.
        for (std::size_t i = 0; i < voiced_track.frames.size(); ++i) {
            auto& f = voiced_track.frames[i];
            if (f.hz > 0.) f.hz = pitch::sample(track, static_cast<double>(i)*seed.step).hz;
        }
    }
    const auto voiced = pitch::locate(x, sr, voiced_track);
    auto average = voiced.empty() ? std::vector<float>{} : cycle_average(x, voiced);
    if (!average.empty()) return_breath(x, sr, marks, voiced, average, amount);
    std::vector<float> high_band(average.size());
    const double floor = 1e-8 * peak*peak * n*n + 1e-17;
    const double df = double(sr)/n;
    auto island = voiced.begin();
    Alignment aligned(marks, len);
    // Four overlapping windows share the same phase-aligned samples. Keep
    // only the current window in a ring and interpolate each point once.
    struct AlignedSamples { float left, right, wl, wr; };
    std::vector<AlignedSamples> neighbors(n);
    int filled = 0;
    std::vector<Complex> current(half), prediction(half), previous(half);
    bool any = false;
    for (int center = -hop; center-half < len; center += hop) {
        const int start = center-half;
        while (island != voiced.end() && start >= island->marks.back()) ++island;
        const bool has_average = island != voiced.end() && start+n > island->marks.front();
        const auto clock = pitch::sample(track, center);
        const double hz = clock.hz;
        const double anchor_stop = 2.3*hz + 2.*df;
        const double anchor_taper = hz + df;
        for (; filled < std::min(len, start+n); ++filled) {
            const auto p = aligned.at(filled);
            neighbors[filled & (n-1)] = {cubic(x, p.left), cubic(x, p.right), p.wl, p.wr};
        }
        double el = 0., er = 0., wav = 0.;
        for (int j = 0; j < n; ++j) {
            const int q = start+j;
            if (q < 0 || q >= len) { pair[j] = next[j] = Complex(); continue; }
            const auto& p = neighbors[q & (n-1)];
            const float w = win[j], a = x[q]*w;
            pair[j] = Complex(a, w*p.wl*p.left);
            next[j] = Complex(w*p.wr*p.right, has_average ? w*average[q] : 0.f);
            el += w*w*p.wl*p.wl; er += w*w*p.wr*p.wr; wav += w*w;
        }
        if (el+er <= 1e-12 || wav <= 1e-12) continue;
        const float left = static_cast<float>(std::sqrt(el/wav)), right = static_cast<float>(std::sqrt(er/wav));
        fft.run(pair);
        fft.run(next);
        std::array<Evidence, bands> scores{};
        Evidence anchor;
        const bool paired = left > .6f && right > .6f;
        for (int k = 1; k < half; ++k) {
            conjugate_pair(pair, k, n, current[k], previous[k]);
            const auto l = previous[k];
            Complex r = next[k], high{};
            if (has_average) conjugate_pair(next, k, n, r, high);
            next[k] = r;
            next[n-k] = high;
            const auto& w = weights[k];
            const BinEvidence bin(current[k], l, r, left, right, paired);
            scores[w.a].add((1.-w.w)*w.dc, bin, paired);
            if (w.w > 0) scores[w.b].add(w.w*w.dc, bin, paired);
            const double low_support = 1. - ease((df*k-anchor_stop)/anchor_taper);
            if (low_support > 0.) anchor.add(low_support*w.dc, bin, paired);
            // Include the center in the estimate only after verifying the
            // independent neighbors. It preserves rapid harmonic envelope
            // motion without transferring a band of raw noise wholesale.
            prediction[k] = (.46f * current[k] + .27f * (l+r)) / (.46f+.27f*(left+right));
        }
        std::array<float, bands> mask{};
        float maximum = 0.f;
        double power = 0., coherent = 0.;
        for (const auto& band : scores) power += band.power;
        // A resonant patch of noise can repeat at the trial lag in isolation.
        // Require compatible low orders before treating it as part of a voice.
        // This includes H2/H3, so a weak or missing H1 need not veto the split.
        const double anchor_floor = std::max(floor, power*9e-7);
        const float supported = anchor.mask(-1, anchor_floor);
        // The phase clock may continue across pure frication. The correction
        // needs a measured or supplied period as well as the acoustic tests.
        const auto period_support = [](double q) { return ease((q-.50)/.24); };
        const float changing = period_support(clock.confidence)*anchor.transition(anchor_floor, left, right);
        const double probe = sr/std::max(30., hz) + .5*hop;
        const auto same_side = [&](double at) {
            const auto f = pitch::sample(track, at);
            const double mismatch = std::abs(f.hz/std::max(30., hz)-1.);
            return ease((f.confidence-.55)/.27) * (1.f-ease((mismatch-.048)/.075));
        };
        const double one_side_floor = std::max(anchor_floor, power*7.5e-5);
        const bool interior = start >= 0 && start+n <= len;
        const float on_left = interior ? same_side(center-probe)*anchor.one_sided(one_side_floor, left, right, true) : 0.f;
        const float on_right = interior ? same_side(center+probe)*anchor.one_sided(one_side_floor, left, right, false) : 0.f;
        std::array<float, bands> recovery{}, left_gain{}, right_gain{}, dry_gain{};
        float recovered = 0.f, borrowed_left = 0.f, borrowed_right = 0.f;
        for (int b = 0; b < bands; ++b) {
            const float local = scores[b].mask(b, floor) * supported;
            mask[b] = local * ease((local-.025)/.075);
            maximum = std::max(maximum, mask[b]);
            coherent += scores[b].power * mask[b];
            if (b >= 4) continue;
            const auto& evidence = scores[b];
            const float two = changing > 0.f ? std::min(changing, evidence.transition(floor, left, right)) : 0.f;
            const float one_l = on_left > 0.f ? std::min(on_left, evidence.one_sided(floor, left, right, true)) : 0.f;
            const float one_r = on_right > 0.f ? std::min(on_right, evidence.one_sided(floor, left, right, false)) : 0.f;
            const float one = std::max(one_l, one_r);
            recovery[b] = std::max(0.f, .98f*std::max(two, one)-mask[b]);
            if (recovery[b] <= 0.f) continue;
            if (one > two) {
                // Keep the current envelope when only one prediction is reliable.
                const bool use_left = one_l >= one_r;
                (use_left ? left_gain[b] : right_gain[b]) = .38f*recovery[b]*evidence.neighbor_gain(use_left);
                dry_gain[b] = .62f*recovery[b];
                auto& borrowed = use_left ? borrowed_left : borrowed_right;
                borrowed = std::max(borrowed, one);
            } else {
                const auto gain = evidence.neighbor_gains();
                left_gain[b] = .65f*recovery[b]*gain[0];
                right_gain[b] = .65f*recovery[b]*gain[1];
                dry_gain[b] = .35f*recovery[b];
            }
            recovered = std::max(recovered, recovery[b]);
        }
        if (center >= 0 && static_cast<std::size_t>(center/hop) < count) {
            // Keep activity and sustain selection tied to measured periodicity.
            // Neither depends on the setting of the upper-band mask.
            out.strength[center/hop] = maximum;
            out.periodic_share[center/hop] = static_cast<float>(coherent/std::max(floor, power));
        }
        const bool extract_low = maximum > 0.f;
        if (!extract_low && recovered <= 0.f && !has_average) continue;
        any |= extract_low || recovered > 0.f;
        const bool inverse_edge = recovered > 0.f || (extract_low && has_average && cleanup);
        pair[0] = Complex{};
        pair[half] = Complex(0.f, has_average ? next[half].imag() : 0.f);
        const double edge_full = std::clamp(3.15*hz, 420., 1100.);
        const double edge_width = std::clamp(1.2*hz, 140., 450.);
        const float borrow_l = .83f*ease(borrowed_left/.4), borrow_r = .83f*ease(borrowed_right/.4);
        const double margin = std::min<double>(hop, 1.2*sr/std::max(30., hz));
        const auto edge_support = [&](int at) {
            const auto fade = [&](double pos) { return ease((pitch::sample(track, pos).confidence-.30)/.42); };
            float support = fade(at);
            if (borrow_l > support) support = std::max(support, borrow_l*fade(at-margin));
            if (borrow_r > support) support = std::max(support, borrow_r*fade(at+margin));
            return support;
        };
        const float edge_center = recovered > 0.f ? edge_support(center) : 0.f;
        for (int k = 1; k < half; ++k) {
            const auto& w = weights[k];
            const float m = (mask[w.a] + w.w*(mask[w.b]-mask[w.a])) * w.dc;
            const auto h = m * prediction[k];
            const auto low = (1.f-w.upper)*h;
            Complex delta{};
            if (recovered > 0.f && df*k < edge_full+edge_width) {
                const float dry = dry_gain[w.a] + w.w*(dry_gain[w.b]-dry_gain[w.a]);
                const float l = left_gain[w.a] + w.w*(left_gain[w.b]-left_gain[w.a]);
                const float r = right_gain[w.a] + w.w*(right_gain[w.b]-right_gain[w.a]);
                const float taper = 1.f - ease((df*k-edge_full)/edge_width);
                delta = (w.dc*taper)*(dry*current[k]+l*previous[k]+r*next[k]);
            }
            const float suppression = cleanup ? .50f*w.low*m*m : 0.f;
            const auto residue = suppression*(current[k]-h-edge_center*delta);
            // Reuse the packed transform to keep the high band available
            // for the final sample-wise fade. In dry frames it carries cleanup.
            const auto second = has_average ? w.upper*next[n-k] : residue;
            pack_pair(pair, k, n, low, second);
            if (inverse_edge) {
                const auto clean = has_average ? residue : Complex{};
                pack_pair(edge_spectrum, k, n, delta, clean);
            }
        }
        fft.run(pair, true);
        if (inverse_edge) {
            edge_spectrum[0] = edge_spectrum[half] = Complex();
            fft.run(edge_spectrum, true);
        }
        const float gain = .5f/n;
        for (int j = std::max(0,-start); j < std::min(n,len-start); ++j) {
            const int at = start+j;
            const float w = gain * win[j];
            if (extract_low) out.harmonic[at] += w*pair[j].real();
            if (has_average) high_band[at] += w*pair[j].imag();
            if (extract_low && cleanup)
                out.cleanup[at] += w*(has_average ? edge_spectrum[j].imag() : pair[j].imag());
            if (recovered > 0.f) {
                // A reliable inner cycle can support a weaker boundary.
                out.harmonic[at] += w*edge_support(at)*edge_spectrum[j].real();
            }
        }
    }
    if (!any) {
        out.harmonic.clear();
        out.cleanup.clear();
    } else if (!high_band.empty()) {
        // The high band still needs real low-order voicing evidence.
        // Fade its filter tails before they reach unvoiced frication.
        for (const auto& region : voiced) {
            const auto& p = region.marks;
            const int first = static_cast<int>(std::ceil(p.front()));
            const int last = std::min(len, static_cast<int>(std::ceil(p.back())));
            for (int i = first; i < last; ++i) {
                const float edge = std::min(ease(2.*(i-p.front())/(p[1]-p.front())),
                                            ease(2.*(p.back()-i)/(p.back()-p[p.size()-2])));
                out.harmonic[i] += edge*ease((out.activity(i)-.025)/.175)*high_band[i];
            }
        }
    }
    return out;
}
} // namespace tn::detail
