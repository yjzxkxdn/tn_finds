// SPDX-License-Identifier: GPL-3.0-only
#include "phrase.hpp"
#include "dsp_fft.hpp"
#include "fractional_sinc.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>

namespace tn::wavtool {
namespace {
using detail::Complex;
constexpr double pi = 3.14159265358979323846;

double smooth(double t) {
    t = std::clamp(t, 0., 1.);
    return t*t*(3.-2.*t);
}

struct Position {
    std::size_t a, b;
    float u;
};

Position at_grid(double at, int hop, int active, std::size_t count) {
    const auto a = std::min(static_cast<std::size_t>(at/hop), count-2);
    const double left = double(a)*hop, right = std::min(left+hop, double(active));
    return {a, a+1, static_cast<float>(std::clamp((at-left)/(right-left), 0., 1.))};
}

double source_at(const Clip& clip, int at, int sr) {
    return std::round(clip.offset*clip.note.sample_rate/1000.) +
           double(at-clip.left)*clip.note.sample_rate/sr;
}

struct Spectrum {
    double hz = 0.;
    std::vector<Complex> orders;
};

Complex order(const std::vector<Complex>& values, std::size_t k) {
    return k < values.size() ? values[k] : Complex{};
}

double angle(Complex coefficient) {
    return std::atan2(double(coefficient.imag()), double(coefficient.real()));
}

Complex polar(double magnitude, double phase) {
    return {static_cast<float>(magnitude*std::cos(phase)), static_cast<float>(magnitude*std::sin(phase))};
}

// Every partial has its own amplitude envelope and a phase relative to H1.
// Straight-line complex interpolation would put opposed partials through zero.
Complex interpolate_order(Complex a, Complex b, double level_mix, double phase_mix) {
    constexpr float negligible = 1.e-12f;
    const double x = std::sqrt(std::norm(a)), y = std::sqrt(std::norm(b));
    if (x <= negligible) return b*static_cast<float>(level_mix);
    if (y <= negligible) return a*static_cast<float>(1.-level_mix);
    if (a.imag() == 0.f && b.imag() == 0.f && a.real() > 0.f && b.real() > 0.f)
        return {static_cast<float>(x+level_mix*(y-x)), 0.f};
    const double phase = angle(a);
    return polar(x+level_mix*(y-x), phase+phase_mix*std::remainder(angle(b)-phase, 2*pi));
}

Complex add_order(Complex a, Complex b, double& arc) {
    const double x = std::sqrt(std::norm(a)), y = std::sqrt(std::norm(b));
    if (x <= 1.e-12) return b;
    if (y <= 1.e-12) return a;
    const double phase = angle(a), difference = angle(b)-phase;
    // Retain the same rotation through +/-pi as long as this pair plays.
    arc = std::isfinite(arc) ? arc+std::remainder(difference-arc, 2*pi)
                             : std::remainder(difference, 2*pi);
    return polar(x+y, phase+arc*y/(x+y));
}

class ClipSweep {
    struct Playing {
        std::size_t clip, prior;
        std::vector<double> phase_arcs;
    };
    const std::vector<Clip>& clips_;
    const std::vector<int>& ends_;
    std::vector<std::size_t> starts_;
    std::vector<Playing> playing_;
    std::size_t next_ = 0;
    int sr_;

public:
    ClipSweep(const std::vector<Clip>& clips, const std::vector<int>& ends, int sr)
        : clips_(clips), ends_(ends), sr_(sr) {
        for (std::size_t i = 0; i < clips.size(); ++i)
            if (ends[i] > std::max(clips[i].left, 0) && clips[i].note.frames.size() >= 2)
                starts_.push_back(i);
        std::stable_sort(starts_.begin(), starts_.end(), [&](auto a, auto b) {
            return clips[a].left < clips[b].left;
        });
    }

    Spectrum at(int n) {
        const auto no_clip = clips_.size();
        const double no_arc = std::numeric_limits<double>::quiet_NaN();
        while (next_ < starts_.size() && clips_[starts_[next_]].left <= n)
            playing_.push_back({starts_[next_++], no_clip, {}});
        playing_.erase(std::remove_if(playing_.begin(), playing_.end(),
                       [&](const Playing& p) { return ends_[p.clip] <= n; }), playing_.end());
        Spectrum result;
        double support = 0., carried_hz = 0.;
        int prior_end = n;
        std::size_t prior = no_clip;
        for (auto& p : playing_) {
            const auto& clip = clips_[p.clip];
            const auto& note = clip.note;
            const double local = source_at(clip, n, sr_);
            if (local < 0. || local >= note.active) continue;
            const auto f = at_grid(local, note.pitch_hop, note.active, note.f0.size());
            const double hz = note.f0[f.a] + f.u*(note.f0[f.b]-note.f0[f.a]);
            carried_hz = hz;
            const auto grid = at_grid(local, note.hop, note.active, note.frames.size());
            const auto& a = note.frames[grid.a];
            const auto& b = note.frames[grid.b];
            const double mix = smooth(grid.u);
            const double voice = a.support + mix*(b.support-a.support);
            if (voice <= 1.e-6) {
                p.phase_arcs.clear();
                continue;
            }
            double gain = clip.gain(n);
            if (note.active < static_cast<int>(note.residual.size()))
                gain *= smooth((note.active-1.-local)/(.003*note.sample_rate));
            if (result.hz == 0.) {
                result.hz = hz;
                support = voice;
            } else {
                const double span = prior_end-clip.left;
                double u = span > 0. ? std::clamp((n-clip.left)/span, 0., 1.) : 1.;
                const double x = (1.-u)*support, y = u*voice;
                if (x+y > 0.) u = y/(x+y);
                result.hz += u*(hz-result.hz);
                support = std::max(support, voice);
            }
            const auto count = std::max(a.harmonics.size(), b.harmonics.size());
            if (result.orders.size() < count) result.orders.resize(count);
            if (p.prior != prior) {
                p.phase_arcs.clear();
                p.prior = prior;
            }
            if (prior != no_clip && p.phase_arcs.size() < count) p.phase_arcs.resize(count, no_arc);
            const auto level = static_cast<float>(gain);
            for (std::size_t k = 0; k < count; ++k) {
                const auto value = level*interpolate_order(order(a.harmonics,k), order(b.harmonics,k), mix, grid.u);
                if (prior == no_clip) result.orders[k] = value;
                else result.orders[k] = add_order(result.orders[k], value, p.phase_arcs[k]);
            }
            prior = p.clip;
            prior_end = std::max(prior_end, ends_[p.clip]);
        }
        // If confidence drops inside a model clip, keep advancing the known
        // pitch clock through the residual without adding any harmonics.
        if (result.hz == 0.) result.hz = carried_hz;
        return result;
    }
};

void bandlimit(Spectrum& frame, double max_hz, int sr) {
    if (max_hz <= 0.) { frame.orders.clear(); return; }
    frame.orders.resize(std::min(frame.orders.size(), static_cast<std::size_t>(.475*sr/max_hz)));
    for (std::size_t k = 0; k < frame.orders.size(); ++k) {
        const auto edge = static_cast<float>(1.-smooth(((k+1)*max_hz/sr-.425)/.05));
        frame.orders[k] *= edge;
    }
}

int polar_steps(const Spectrum& a, const Spectrum& b, int length) {
    float peak = 0.f;
    for (const auto* frame : {&a, &b})
        for (const auto coefficient : frame->orders)
            peak = std::max(peak, std::norm(coefficient));
    double turn = 0.;
    const double min_power = std::max(1.e-24, 6.25e-6*peak);
    for (std::size_t k = 0; k < std::min(a.orders.size(), b.orders.size()); ++k) {
        if (std::min(std::norm(a.orders[k]), std::norm(b.orders[k])) <= min_power) continue;
        const double arc = std::remainder(angle(b.orders[k])-angle(a.orders[k]), 2*pi);
        turn = std::max(turn, std::abs(arc));
    }
    // Evaluate extra polar knots for rapid partials. Between the resulting
    // tables a chord of at most pi/8 loses at most 0.17 dB at equal levels.
    return std::clamp(static_cast<int>(std::ceil(turn/(pi/8.))), 1, length);
}

class Tables {
    std::array<std::optional<detail::RealInverseFft>, 14> plans_;
    std::vector<Complex> spectrum_;
public:
    void build(const std::vector<Complex>& orders, std::vector<float>& table) {
        if (orders.empty()) { table.clear(); return; }
        const int count = static_cast<int>(orders.size());
        int size = 32, bits = 5;
        while (size < 8*(count+1)) { size *= 2; ++bits; }
        auto& fft = plans_[bits];
        if (!fft) fft.emplace(size);
        spectrum_.assign(size/2+1, {});
        std::copy(orders.begin(), orders.end(), spectrum_.begin()+1);
        fft->run(spectrum_, table);
    }
};

float cycle_at(const std::vector<float>& table, double phase) {
    if (table.empty()) return 0.f;
    const int mask = static_cast<int>(table.size())-1;
    const double at = phase*table.size();
    const int k = static_cast<int>(at);
    const auto t = static_cast<float>(at-k);
    const float a = table[(k-1)&mask], b = table[k&mask], c = table[(k+1)&mask], d = table[(k+2)&mask];
    return b + .5f*t*(c-a + t*(2*a-5*b+4*c-d + t*(3*(b-c)+d-a)));
}
} // namespace

double Clip::gain(int at) const {
    const double t = (at-left)/frames_per_ms;
    std::size_t k = 0;
    while (k+2 < positions.size() && t >= positions[k+1]) ++k;
    const double width = positions[k+1]-positions[k];
    const double mix = width <= 0. ? 0. : std::clamp((t-positions[k])/width, 0., 1.);
    return gains[k]+mix*(gains[k+1]-gains[k]);
}

std::vector<float> synthesize_phrase(const std::vector<Clip>& clips, int size, int sr) {
    std::vector<float> result(size);
    std::vector<int> ends(clips.size());
    int end = size;
    for (std::size_t i = clips.size(); i-- > 0;) {
        end = std::min(end, clips[i].end);
        ends[i] = end;
    }
    bool harmonics = false;
    for (std::size_t i = 0; i < clips.size(); ++i) {
        const auto& clip = clips[i];
        const auto& note = clip.note;
        const int first = std::max(0, clip.left);
        if (ends[i] <= first || note.residual.empty()) continue;
        const int available = static_cast<int>(note.residual.size());
        std::optional<detail::Sinc> sinc;
        if (note.sample_rate != sr) sinc.emplace(note.residual, double(note.sample_rate)/sr);
        harmonics = harmonics || !note.frames.empty();
        const double offset = std::round(clip.offset*note.sample_rate/1000.);
        if (!sinc) {
            // At equal sample rates the original rounded offset and source
            // cursor are integers. Clip to the recorded interval once.
            const double origin = clip.left-offset;
            const double lo = std::max(double(first), origin);
            const double hi = std::min(double(ends[i]), origin+available);
            if (lo >= hi) continue;
            int at = static_cast<int>(lo-origin);
            const int stop = static_cast<int>(hi);
            for (int n = static_cast<int>(lo); n < stop; ++n, ++at)
                result[n] += static_cast<float>(clip.gain(n)*note.residual[at]);
        } else {
            for (int n = first; n < ends[i]; ++n) {
                const double at = offset+double(n-clip.left)*note.sample_rate/sr;
                if (at < 0. || at >= available) continue;
                result[n] += static_cast<float>(clip.gain(n)*sinc->read(at, 0, available));
            }
        }
    }
    if (!harmonics || size == 0) return result;

    // Insert clip boundaries into the synthesis grid. An interval outside
    // all model clips contributes no harmonic audio, including real rests.
    std::vector<int> boundaries{0, size};
    std::vector<int> coverage(size+1);
    for (std::size_t i = 0; i < clips.size(); ++i) {
        const auto& clip = clips[i];
        const auto& note = clip.note;
        if (note.frames.empty()) continue;
        const double offset = std::round(clip.offset*note.sample_rate/1000.);
        const double ratio = double(sr)/note.sample_rate;
        const auto limit = [&](double local) {
            return static_cast<int>(std::ceil(std::clamp(clip.left+(local-offset)*ratio, 0., double(size))));
        };
        const int first = std::max({0, clip.left, limit(0.)});
        const int last = std::min(ends[i], limit(note.active));
        if (first >= last) continue;
        boundaries.push_back(first);
        boundaries.push_back(last);
        ++coverage[first];
        --coverage[last];
    }
    const int hop = std::max(1, static_cast<int>(std::lround(.0029*sr)));
    for (int n = hop; n < size; n += hop) boundaries.push_back(n);
    std::sort(boundaries.begin(), boundaries.end());
    boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());
    ClipSweep sweep(clips, ends, sr);
    auto a = sweep.at(0), b = sweep.at(boundaries[1]);
    Tables tables;
    std::vector<float> left, right, inner_left, inner_right;
    std::vector<Complex> intermediate;
    bandlimit(a, std::max(a.hz,b.hz), sr);
    tables.build(a.orders, left);
    double phase = 0.;
    int covered = 0;
    for (std::size_t i = 1; i < boundaries.size(); ++i) {
        const int first = boundaries[i-1], last = boundaries[i], len = last-first;
        covered += coverage[first];
        auto next = i+1 < boundaries.size() ? sweep.at(boundaries[i+1]) : Spectrum{};
        // Use the same bandlimited spectrum and table on both sides of each knot.
        bandlimit(b, std::max({a.hz, b.hz, next.hz}), sr);
        tables.build(b.orders, right);
        if (covered > 0 && (a.hz > 0. || b.hz > 0.)) {
            const double f0 = a.hz > 0. ? a.hz : b.hz;
            const double f1 = b.hz > 0. ? b.hz : f0;
            const double step = (f1-f0)/len;
            const int steps = polar_steps(a,b,len);
            for (int part = 1, from = 0; part <= steps; ++part) {
                const int to = len*part/steps;
                const auto& ta = part == 1 ? left : inner_left;
                const auto* tb = &right;
                if (part < steps) {
                    const double u = double(to)/len;
                    intermediate.resize(std::max(a.orders.size(), b.orders.size()));
                    for (std::size_t k = 0; k < intermediate.size(); ++k)
                        intermediate[k] = interpolate_order(order(a.orders,k), order(b.orders,k), u, u);
                    tables.build(intermediate, inner_right);
                    tb = &inner_right;
                }
                for (int j = from; j < to; ++j) {
                    const auto u = static_cast<float>(j-from)/(to-from);
                    result[first+j] += (1.f-u)*cycle_at(ta, phase) + u*cycle_at(*tb, phase);
                    phase += (f0+(j+.5)*step)/sr;
                    phase -= std::floor(phase);
                }
                std::swap(inner_left, inner_right);
                from = to;
            }
        }
        a = std::move(b);
        b = std::move(next);
        std::swap(left,right);
    }
    return result;
}
} // namespace tn::wavtool
