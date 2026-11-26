// SPDX-License-Identifier: GPL-3.0-only
#include "harmonic_ridges.hpp"
#include "dsp_fft.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace tn::detail {
namespace {
constexpr double pi = 3.14159265358979323846;

float smooth(double u) {
    u = std::clamp(u, 0., 1.);
    return static_cast<float>(u*u*(3.-2.*u));
}

float interpolate(const std::vector<float>& x, double p) {
    if (p < 0. || p >= static_cast<double>(x.size())) return 0.f;
    const auto i = static_cast<std::size_t>(p);
    const float u = static_cast<float>(p-i);
    const float b = x[i], c = i+1 < x.size() ? x[i+1] : 0.f;
    if (i == 0 || i+2 >= x.size()) return b+u*(c-b);
    const float a = x[i-1], d = x[i+2];
    return b+.5f*u*(c-a+u*(2*a-5*b+4*c-d+u*(3*(b-c)+d-a)));
}

// A few dozen microseconds of lag error do not matter much to a low F0
// estimate, but move a high partial by an appreciable fraction of a cycle.
// For ridge analysis only, run a fine fractional lock on a zero phase
// low/mid-band copy. Never replace the marks used by the actual renderer.
std::vector<double> refine_ridge_marks(const std::vector<float>& x, int sr,
                                       const std::vector<double>& marks) {
    if (std::adjacent_find(marks.begin(), marks.end(), [=](double a, double b) {
        return b-a > sr/155.;
    }) == marks.end()) return {};
    const int dec = std::max(1, sr/7600);
    const int radius = 5*dec;
    const double cut = std::min(1850./sr, .41/dec);
    std::vector<float> fir(2*radius+1);
    double sum = 0.;
    for (int j = -radius; j <= radius; ++j) {
        const double u = pi*j/radius;
        sum += fir[j+radius] = static_cast<float>((.42+.5*std::cos(u)+.08*std::cos(2.*u)) *
                         (j ? std::sin(2.*pi*cut*j)/(pi*j) : 2.*cut));
    }
    for (auto& v : fir) v = static_cast<float>(v/sum);
    const int length = static_cast<int>(x.size());
    std::vector<float> low(static_cast<std::size_t>((length-1)/dec+1));
    for (int i = 0; i < static_cast<int>(low.size()); ++i) {
        const int p = i*dec, begin = std::max(-radius,-p), end = std::min(radius,length-1-p);
        float v = 0.;
        for (int j = begin; j <= end; ++j) v += fir[j+radius]*x[p+j];
        low[i] = v;
    }
    constexpr int count = 81;
    std::array<float, count> window{}, saved{};
    float total = 0.;
    for (int j = 0; j < count; ++j)
        total += window[j] = static_cast<float>(.63-.37*std::cos(2.*pi*j/(count-1)));
    std::vector<float> delta(marks.size()), support(marks.size());
    for (std::size_t i = 1; i < marks.size(); ++i) {
        const double a = marks[i-1]/dec, b = marks[i]/dec, span = b-a;
        const float pitch_weight = smooth((span*dec/sr-1./155.)/(1./120.-1./155.));
        const double prev = i >= 2 ? .5*(span+(marks[i-1]-marks[i-2])/dec) : span;
        const double next = i+1 < marks.size() ? .5*(span+(marks[i+1]-marks[i])/dec) : span;
        const double range = std::min(.015*span,.00021*sr/dec);
        const double a_start = a-.41*prev, b_start = b-.41*next;
        if (pitch_weight <= 0.f || a_start < 1. || b_start-range < 1. ||
            b+.41*next+range >= static_cast<double>(low.size())-2.) continue;
        const double da = .82*prev/(count-1), db = .82*next/(count-1);
        double mass = 0., energy = 0.;
        for (int j = 0; j < count; ++j) {
            const double v = saved[j] = interpolate(low,a_start+j*da);
            mass += window[j]*v;
            energy += window[j]*v*v;
        }
        energy -= mass*mass/total;
        if (energy < 1e-11) continue;
        const auto agree = [&](double offset) {
            double mb = 0., eb = 0., dot = 0.;
            for (int j = 0; j < count; ++j) {
                const double v = interpolate(low,b_start+j*db+offset);
                mb += window[j]*v;
                eb += window[j]*v*v;
                dot += window[j]*saved[j]*v;
            }
            eb -= mb*mb/total;
            return eb > 1e-11 ? (dot-mass*mb/total)/std::sqrt(energy*eb) : 0.;
        };
        const double center = agree(0.);
        if (center < .94) continue;
        const double l = agree(-range), r = agree(range);
        const double curvature = 2.*center-l-r;
        if (curvature <= 5e-5) continue;
        delta[i] = static_cast<float>(dec*range*std::clamp(.5*(r-l)/curvature,-1.,1.));
        support[i] = pitch_weight*smooth((center-.94)/.055);
    }
    // Bounded leaky integration in both directions: correct the differences
    // between adjacent marks without shifting the whole note or accumulating
    // drift when the source has poor acoustic evidence.
    auto out = marks;
    const float max_shift = .00028f*sr;
    float shift = 0.f;
    for (std::size_t i = 1; i < marks.size(); ++i) {
        shift = std::clamp(support[i]*(.82f*shift+delta[i]),-max_shift,max_shift);
        out[i] += .5f*shift;
    }
    shift = 0.f;
    for (std::size_t i = marks.size()-1; i > 0; --i) {
        shift = std::clamp(support[i]*(.82f*shift-delta[i]),-max_shift,max_shift);
        out[i-1] += .5f*shift;
    }
    return out;
}

// Linear within each measured cycle. The same position mapping is used on
// analysis and synthesis, including a small extension at each real file edge.
class PhaseMap {
    const std::vector<double>& marks_;
    int last_;
public:
    explicit PhaseMap(const std::vector<double>& m) : marks_(m), last_(static_cast<int>(m.size())-1) {}
    double position(int i) const {
        if (i < 0) return marks_[0]+i*(marks_[1]-marks_[0]);
        if (i > last_) return marks_.back()+(i-last_)*(marks_[last_]-marks_[last_-1]);
        return marks_[i];
    }
    int cycle(double pos) const {
        if (pos < marks_.front()) return static_cast<int>(std::floor((pos-marks_[0])/(marks_[1]-marks_[0])));
        if (pos >= marks_.back())
            return last_+static_cast<int>(std::floor((pos-marks_.back())/(marks_[last_]-marks_[last_-1])));
        return static_cast<int>(std::upper_bound(marks_.begin(), marks_.end(), pos)-marks_.begin())-1;
    }
    double phase(double pos) const {
        const int i = cycle(pos);
        const double a = position(i);
        return i+(pos-a)/(position(i+1)-a);
    }
    double time(double phase) const {
        const int i = static_cast<int>(std::floor(phase));
        const double a = position(i);
        return a+(phase-i)*(position(i+1)-a);
    }
};

struct Work {
    Fft fft;
    RealInverseFft inverse;
    std::vector<float> window, power, keep, removed;
    std::vector<Complex> frames, average;
    explicit Work(int n) : fft(n), inverse(n), window(n), power(n/2), keep(n/2), removed(n),
                           frames(n), average(n/2) {
        for (int i = 0; i < n; ++i) window[i] = static_cast<float>(std::sin(pi*i/n));
    }
};

struct Peak {
    float bin = 0.f, weight = 0.f;
};

Peak partial(const Work& w, int order, int cycles, float floor) {
    const int expected = order*cycles;
    const int spread = cycles/3;
    int b = expected;
    for (int i = expected-spread; i <= expected+spread; ++i)
        if (w.power[i] > w.power[b]) b = i;
    const float l = w.power[b-1], m = w.power[b], r = w.power[b+1];
    if (!(m > floor && m >= l && m >= r)) return {};
    const float curvature = 2*m-l-r;
    const float bin = b+(curvature > 1e-20f ? .5f*(r-l)/curvature : 0.f);
    const float left = .5f*(w.power[expected-cycles/2] + w.power[expected-cycles/2+1]);
    const float right = .5f*(w.power[expected+cycles/2] + w.power[expected+cycles/2-1]);
    // Compare the peak to both nearby between-harmonic bands. A broadband
    // consonant is unlikely to supply several peaks at these cycle phases.
    const float background = std::max(floor, .5f*(left+right));
    const float contrast = (m-background)/(m+background);
    const float near = 1.f-smooth((std::abs(bin-expected)/cycles-.19)/.19);
    const float weight = smooth((contrast-.58)/.31)*near;
    return {bin, weight};
}

} // namespace

void return_breath(const std::vector<float>& input, int sr, const std::vector<double>& marks,
                   const std::vector<pitch::Island>& voiced, std::vector<float>& average, double amount) {
    if (amount <= 0. || voiced.empty()) return;
    const int len = static_cast<int>(input.size());
    const int hop = std::max(1, static_cast<int>(std::lround(.008 * sr)));
    const auto refined = refine_ridge_marks(input, sr, marks);
    const PhaseMap phase(refined.empty() ? marks : refined);
    std::vector<float> removed(len), normalizer(len);
    struct Scratch {
        std::array<std::unique_ptr<Work>, 15> plans{};
        std::vector<Peak> peaks;
        std::vector<std::array<double, 2>> comb;
    };
    struct Block {
        const Work* work = nullptr;
        int start = 0, stop = 0;
        double left = 0., phase_step = 0.;
    };
    const auto analyze = [&](int center, Scratch& scratch) -> Block {
        const double middle = phase.phase(center);
        const double period = (phase.time(middle+2.)-phase.time(middle-2.))*.25;
        // Shorter bass windows limit phase drift without merging the ridges.
        const double local_hz = sr/period;
        const double min_cycles = 8.-2.*smooth((148.-local_hz)/28.);
        const int cycles = static_cast<int>(std::lround(std::clamp(.050*local_hz,min_cycles,16.)));
        const double left = middle-.5*cycles, right = left+cycles;
        const double begin = phase.time(left), end = phase.time(right);
        const int start = std::max(0, static_cast<int>(std::ceil(begin)));
        const int stop = std::min(len, static_cast<int>(std::ceil(end)));
        if (stop <= start) return {};
        const double f0 = sr*cycles/(end-begin);
        int exponent = 8;
        while (exponent < 14 && (1<<exponent) < 1.10*(end-begin)) ++exponent;
        const int n = 1<<exponent, half = n/2;
        auto& plan = scratch.plans[exponent];
        if (!plan) plan = std::make_unique<Work>(n);
        auto& w = *plan;
        const double phase_step = double(cycles)/n;
        // With a high supplied F0 the fundamental itself can enter this band.
        const int first_order = std::max(1, static_cast<int>(std::floor(1050./f0)));
        const int last_order = std::min((half-1-cycles)/cycles, static_cast<int>(std::min(21000.,sr*.49)/f0));
        int cycle = static_cast<int>(std::floor(left));
        double a = phase.position(cycle), b = phase.position(cycle+1);
        for (int j = 0; j < n; ++j) {
            const double p = left+j*phase_step;
            if (p >= cycle+1) {
                a = b;
                b = phase.position(++cycle+1);
            }
            const double at = a+(p-cycle)*(b-a);
            w.frames[j] = w.window[j]*Complex(interpolate(input, at), interpolate(average, at));
        }
        w.fft.run(w.frames);
        float frame_max = 0.f;
        for (int k = 1; k < half; ++k) {
            const auto a0 = w.frames[k], a1 = std::conj(w.frames[n-k]);
            w.power[k] = .25f*std::norm(a0+a1);
            w.average[k] = Complex(0.f,-.5f)*(a0-a1);
            frame_max = std::max(frame_max,w.power[k]);
        }
        const float floor = std::max(1e-16f, frame_max*3e-8f);
        auto& peaks = scratch.peaks;
        auto& comb = scratch.comb;
        peaks.resize(static_cast<std::size_t>(last_order+1));
        comb.resize(peaks.size());
        const int core = std::max(1, cycles/5);
        for (int order = first_order; order <= last_order; ++order) {
            peaks[order] = partial(w, order, cycles, floor);
            double on = 0., between = 0.;
            for (int k = order*cycles-core; k <= order*cycles+core; ++k) {
                on += w.power[k];
                between += .5*(w.power[k-cycles/2]+w.power[k+(cycles+1)/2]);
            }
            comb[order] = {on, between};
        }
        std::fill(w.keep.begin(), w.keep.end(), 0.f);
        const float bass_weight = 1.f+.9f*smooth((148.-f0)/34.);
        for (int order = first_order; order <= last_order; ++order) {
            const auto peak = peaks[order];
            if (peak.weight <= 0.f) continue;
            // Random narrow peaks also occur in breath. Look for excess
            // energy near several neighboring multiples of the same F0.
            double on = 0., between = 0.;
            for (int q = std::max(first_order,order-8); q <= std::min(last_order,order+8); ++q) {
                on += comb[q][0];
                between += comb[q][1];
            }
            const double contrast = (on-between)/std::max<double>(floor, on+between);
            const float keep = std::min(1.f, 1.8f*bass_weight*peak.weight)*smooth((contrast-.14)/.28);
            const float width = std::min(4.f, .46f*cycles);
            const int first = std::max(1, static_cast<int>(std::ceil(peak.bin-width)));
            const int last = std::min(half-1, static_cast<int>(std::floor(peak.bin+width)));
            for (int k = first; k <= last; ++k) {
                const float shape = 1.f-smooth((std::abs(k-peak.bin)-.42f*width)/(.58f*width));
                w.keep[k] = std::max(w.keep[k], keep*shape);
            }
        }
        // The real mask preserves phase and only subtracts material from the average.
        const double bin_hz = f0/cycles;
        w.frames[0] = w.frames[half] = Complex{};
        for (int k = 1; k < half; ++k) {
            const float high = smooth((k*bin_hz-900.)/350.);
            w.frames[k] = (high*(1.f-w.keep[k]))*w.average[k];
        }
        w.inverse.run(w.frames, w.removed, 1.f/n);
        return {&w, start, stop, left, phase_step};
    };
    const auto overlap = [&](const Block& block) {
        if (!block.work) return;
        const auto& w = *block.work;
        const auto left = block.left, phase_step = block.phase_step;
        int cycle = phase.cycle(block.start);
        for (int i = block.start; i < block.stop;) {
            const double p = phase.position(cycle), q = phase.position(cycle+1);
            const double scale = 1./((q-p)*phase_step);
            const int limit = std::min(block.stop, static_cast<int>(std::ceil(q)));
            for (; i < limit; ++i) {
                const double pos = (cycle-left)/phase_step+(i-p)*scale;
                const float win = interpolate(w.window,pos);
                normalizer[i] += win*win;
                removed[i] += win*interpolate(w.removed,pos);
            }
            ++cycle;
        }
    };
    std::vector<int> centers;
    auto voice = voiced.begin();
    for (int center = 0; center < len; center += hop) {
        while (voice != voiced.end() && center >= voice->marks.back()) ++voice;
        if (voice == voiced.end()) break;
        if (center >= voice->marks.front()) centers.push_back(center);
    }
#ifdef _OPENMP
    // Classic hosts can run many resamplers at once. Bound the work per
    // process and avoid launching a team for ordinary short phonemes.
    const int workers = std::min({4, omp_get_max_threads(), static_cast<int>(centers.size()/32)});
    if (workers > 1 && !omp_in_parallel()) {
#pragma omp parallel num_threads(workers)
        {
            Scratch scratch;
#pragma omp for ordered schedule(static, 1)
            for (int i = 0; i < static_cast<int>(centers.size()); ++i) {
                const auto block = analyze(centers[i], scratch);
                // Preserve the sequence of float additions at each sample.
#pragma omp ordered
                overlap(block);
            }
        }
    } else
#endif
    {
        Scratch scratch;
        for (const int center : centers) overlap(analyze(center, scratch));
    }
    const auto gain = static_cast<float>(amount);
    for (const auto& island : voiced) {
        const auto& p = island.marks;
        const int start = static_cast<int>(std::ceil(p.front()));
        const int stop = std::min(len, static_cast<int>(std::ceil(p.back())));
        for (int i = start; i < stop; ++i) {
            const float edge = std::min(smooth((i-p.front())/(p[1]-p.front())),
                                        smooth((p.back()-i)/(p.back()-p[p.size()-2])));
            average[i] -= gain*edge*removed[i]/std::max(1.f,normalizer[i]);
        }
    }
}
} // namespace tn::detail
