// SPDX-License-Identifier: GPL-3.0-only
#include "pitch.hpp"
#include "dsp_fft.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace tn::pitch {
namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double fmin = 55., fmax = 1300.;
constexpr int hypotheses = 9;  // silence plus up to eight local lag maxima

bool valid_rate(int sr) { return sr >= 8000 && sr <= 192000; }

Frame clean(Frame f) {
    if (!std::isfinite(f.hz) || f.hz <= 0.) f.hz = 0.;
    f.confidence = std::isfinite(f.confidence) ? std::clamp(f.confidence, 0., 1.) : 0.;
    return f;
}

struct Reduced {
    int dec;
    double rate;
    double peak = 0.;
    std::vector<float> y;

    // Symmetric Blackman-windowed sinc: no phase delay or aliasing by simple
    // subsampling. Approximately 10.5 kHz (usually 11--12 kHz) gives
    // high F0 candidates and their strongest harmonics room below Nyquist.
    Reduced(const std::vector<float>& x, int sr, int target_rate = 10500,
            double cutoff_hz = std::numeric_limits<double>::infinity(), bool round_rate = false)
        : dec(std::max(1, (sr + (round_rate ? target_rate/2 : 0)) / target_rate)), rate(double(sr) / dec) {
        const int m = 6 * dec;
        const double cutoff = std::min(.38 / dec, cutoff_hz / sr);
        std::vector<float> taps(2 * m + 1);
        double dc = 0.;
        for (int j = -m; j <= m; ++j) {
            const double p = pi * j / m;
            const double w = .42 + .5 * std::cos(p) + .08 * std::cos(2. * p);
            dc += taps[j + m] = static_cast<float>(w * (j ? std::sin(2. * pi * cutoff * j) / (pi * j) : 2. * cutoff));
        }
        for (float& t : taps) t = static_cast<float>(t / dc);

        double input_peak = 0.;
        for (float v : x) if (std::isfinite(v)) input_peak = std::max(input_peak, std::abs(double(v)));
        const double scale = 1. / std::max(1., input_peak);
        const auto n = static_cast<int>(x.size());
        std::vector<float> padded(x.size() + 2 * m, 0.f);
        for (int i = 0; i < n; ++i) if (std::isfinite(x[i])) padded[i+m] = static_cast<float>(x[i] * scale);
        std::fill(padded.begin(), padded.begin() + m, padded[m]);
        std::fill(padded.end() - m, padded.end(), padded[m+n-1]);
        y.resize(1 + (n - 1) / dec);
        for (std::size_t i = 0; i < y.size(); ++i) {
            const float* in = padded.data() + i * dec;
            float a = 0.f, b = 0.f, c = 0.f, d = 0.f;
            int j = 0;
            for (; j + 3 < 2 * m + 1; j += 4) {
                a += in[j] * taps[j]; b += in[j+1] * taps[j+1];
                c += in[j+2] * taps[j+2]; d += in[j+3] * taps[j+3];
            }
            float u = (a + b) + (c + d);
            for (; j < 2 * m + 1; ++j) u += in[j] * taps[j];
            y[i] = u;
            peak = std::max(peak, std::abs(double(u)));
        }
    }
};

struct Sums {
    const std::vector<float>& y;
    std::vector<double> mass, energy;

    explicit Sums(const std::vector<float>& input) : y(input), mass(y.size()+1), energy(y.size()+1) {
        for (std::size_t i = 0; i < y.size(); ++i) {
            const double v = y[i];
            mass[i+1] = mass[i] + v;
            energy[i+1] = energy[i] + v*v;
        }
    }

    double power(int start, int len) const {
        if (len <= 0) return 0.;
        const double m = mass[start+len] - mass[start];
        return std::max(0., (energy[start+len] - energy[start] - m*m/len) / len);
    }

    float nccf(int center, int lag, int halfwin, double floor) const {
        const int a = lag / 2, b = lag - a;
        const int start = std::max(center - halfwin, a);
        const int stop = std::min(center + halfwin + 1, static_cast<int>(y.size()) - b);
        const int n = stop - start;
        if (n < std::max(12, lag)) return 0.f;
        const int x0 = start - a, x1 = start + b;
        const double m0 = (mass[stop-a] - mass[x0]) / n, m1 = (mass[stop+b] - mass[x1]) / n;
        const double e0 = energy[stop-a] - energy[x0] - n*m0*m0;
        const double e1 = energy[stop+b] - energy[x1] - n*m1*m1;
        if (e0 <= floor*n || e1 <= floor*n) return 0.f;
        const float* p = y.data() + x0;
        const float* q = y.data() + x1;
        float s0 = 0.f, s1 = 0.f, s2 = 0.f, s3 = 0.f;
        int i = 0;
        for (; i + 3 < n; i += 4) {
            s0 += p[i] * q[i]; s1 += p[i+1] * q[i+1];
            s2 += p[i+2] * q[i+2]; s3 += p[i+3] * q[i+3];
        }
        double dot = (s0 + s1) + (s2 + s3);
        for (; i < n; ++i) dot += double(p[i]) * q[i];
        const double norm = std::sqrt(e0*e1);
        const double corr = std::clamp((dot - n*m0*m1) / norm, -1., 1.);
        // Energy imbalance reduces the appeal of a repeated pattern during
        // attacks, but does not exclude a strong crescendo or decrescendo.
        return static_cast<float>(corr * (.74 + .26 * 2. * norm / (e0 + e1)));
    }
};

struct Hypothesis {
    double hz = 0., strength = 0., reward = 0., log_hz = 0.;
};

struct LatticeFrame {
    std::array<Hypothesis, hypotheses> h{};
    std::array<unsigned char, hypotheses> prev{};
    int count = 1;  // index zero denotes unvoiced
};

// A local cosine fit avoids nearly a quarter-sample bias for the shortest
// lags. If the lobe is not locally cosine-shaped, use the parabola instead.
std::pair<double, double> refine(double l, double m, double r) {
    const double quad = .5 * (r - l) / (2.*m - l - r);
    const double curvature = (l + r) / (2.*m);
    if (curvature > -.85 && curvature < .9999) {
        const double w = std::acos(curvature);
        const double b = (r - l) / (2. * std::sin(w));
        const double off = std::atan2(b, m) / w;
        if (std::abs(off) < .63) return {off, std::min(1., std::hypot(m, b))};
    }
    return {quad, std::clamp(m + .25 * (r - l) * quad, 0., 1.)};
}

void add(LatticeFrame& f, Hypothesis v) {
    if (f.count < hypotheses) f.h[f.count++] = v;
    else {
        auto worst = std::min_element(f.h.begin() + 1, f.h.end(), [](auto a, auto b) {return a.reward < b.reward;});
        if (v.reward > worst->reward) *worst = v;
    }
}

// Autocorrelation of real data using two half-length complex transforms.
// The real and imaginary parts hold the even and odd time samples.
class RealAutocorrelation {
    using Complex = detail::Complex;
    detail::Fft fft_;
    std::vector<Complex> packed_, power_, rotation_;

public:
    explicit RealAutocorrelation(int size)
        : fft_(size/2), packed_(size/2), power_(size/2), rotation_(size/2) {
        for (int k = 0; k < size/2; ++k) {
            const double angle = -2.*pi*k/size;
            rotation_[k] = {static_cast<float>(std::cos(angle)), static_cast<float>(std::sin(angle))};
        }
    }

    void run(const std::vector<float>& signal) {
        std::fill(packed_.begin(), packed_.end(), Complex{});
        for (std::size_t i = 0; i < signal.size(); ++i) {
            if (i & 1) packed_[i/2].imag(signal[i]);
            else packed_[i/2].real(signal[i]);
        }
        fft_.run(packed_);
        const auto mask = packed_.size()-1;
        for (std::size_t k = 0; k < packed_.size(); ++k) {
            const auto a = packed_[k], b = std::conj(packed_[(packed_.size()-k)&mask]);
            const auto even = .5f*(a+b), difference = .5f*(a-b);
            const Complex odd(difference.imag(), -difference.real());
            const auto twiddle = rotation_[k];
            const float wr = twiddle.real(), wi = twiddle.imag();
            const Complex rotated(wr*odd.real()-wi*odd.imag(), wr*odd.imag()+wi*odd.real());
            const float common = std::norm(even)+std::norm(odd);
            const float cross = 2.f*(even.real()*rotated.real()+even.imag()*rotated.imag());
            // The sum and difference of the two real power bins are the
            // transformed even and odd autocorrelation, respectively.
            power_[k] = {common+wi*cross, wr*cross};
        }
        fft_.run(power_, true);
    }

    float lag(int n) const { return n & 1 ? power_[n/2].imag() : power_[n/2].real(); }
};

constexpr double octave_cost = .015;

class WindowedPitch {
    int length_, lower_, upper_;
    std::vector<float> window_, segment_, inverse_window_, scores_;
    RealAutocorrelation ac_;
    int previous_start_ = -1;

    static int fft_size(int span) {
        int size = 2;
        while (size < span) size *= 2;
        return size;
    }

public:
    explicit WindowedPitch(const Reduced& r)
        : length_(std::min(static_cast<int>(r.y.size()), 2*static_cast<int>(1.5*r.rate/fmin)+1)),
          lower_(std::max(2, static_cast<int>(r.rate/fmax)-1)),
          upper_(std::min((length_-1)*2/5, static_cast<int>(std::ceil(r.rate/fmin))+1)),
          window_(length_), segment_(length_), inverse_window_(upper_+1), scores_(upper_+1),
          ac_(fft_size(length_+upper_)) {
        for (int i = 0; i < length_; ++i)
            window_[i] = static_cast<float>(.5-.5*std::cos(2.*pi*i/std::max(1, length_-1)));
        ac_.run(window_);
        for (int lag = 0; lag <= upper_; ++lag) {
            const float energy = ac_.lag(lag);
            if (energy > 0.f) inverse_window_[lag] = ac_.lag(0)/energy;
        }
    }

    void score(LatticeFrame& f, int center, const Sums& sums, double rate, double floor) {
        if (upper_ <= lower_+1) return;
        const int local_half = static_cast<int>(std::ceil(.0046*rate));
        const int from = std::max(0, center-local_half);
        const int until = std::min(static_cast<int>(sums.y.size()), center+local_half+1);
        if (sums.power(from, until-from) <= floor) return;

        // Analyze the closest complete window at a clip boundary. A short
        // recording uses its available samples and only supported lags.
        const int start = std::clamp(center-length_/2, 0, static_cast<int>(sums.y.size())-length_);
        if (start != previous_start_) {
            const auto mean = static_cast<float>((sums.mass[start+length_]-sums.mass[start])/length_);
            for (int i = 0; i < length_; ++i) segment_[i] = window_[i]*(sums.y[start+i]-mean);
            ac_.run(segment_);
            const float energy = ac_.lag(0);
            if (energy <= 0.f) return;
            const float inv_energy = 1.f/energy;
            // Dividing by the window's ACF makes lag scores comparable.
            for (int lag = lower_; lag <= upper_; ++lag)
                scores_[lag] = ac_.lag(lag)*inv_energy*inverse_window_[lag];
            previous_start_ = start;
        }

        for (int lag = lower_+1; lag < upper_; ++lag) {
            const double l = scores_[lag-1], m = scores_[lag], r = scores_[lag+1];
            if (m < .45 || m <= l || m < r) continue;
            const auto [off, peak] = refine(l, m, r);
            const double hz = rate/(lag+off);
            if (hz < fmin*.985 || hz > fmax*1.018) continue;
            const double strength = std::min(peak, 1./std::max(1., m));
            // Check local periodicity without using it to rank lag candidates.
            const int halfwin = std::max(local_half, (lag+1)/2);
            if (center-lag/2-halfwin >= 0 &&
                center+(lag+1)/2+halfwin < static_cast<int>(sums.y.size()) &&
                sums.nccf(center, lag, halfwin, floor) < .44f) continue;
            const double h = std::clamp(hz, fmin, fmax), lh = std::log2(h);
            const double reward = strength-.64-octave_cost*(std::log2(fmax)-lh);
            add(f, {h, strength, reward, lh});
        }
    }
};

double transition(const Hypothesis& from, const Hypothesis& to) {
    if (from.hz == 0.) return to.hz == 0. ? 0. : .24;
    if (to.hz == 0.) return .24;
    const double change = std::abs(from.log_hz - to.log_hz);
    return std::min(.85, .64*change);
}

// Integrate the linearly interpolated F0 within an island, in cycles and
// original sample coordinates. The outer half-frame margins are constant.
struct PhaseCurve {
    std::vector<double> x, rate, phase;

    PhaseCurve(const Track& tr, std::size_t b, std::size_t e, int sr, double begin, double end) {
        x.reserve(e-b+2); rate.reserve(e-b+2);
        // Clip interpolation to this island's centers, even if a rejected
        // out-of-range positive Hz (rather than 0 Hz) precedes or follows it.
        const auto hz = [&](double sample_pos) {
            const double t = std::clamp(sample_pos/tr.step, double(b), double(e-1));
            const auto k = static_cast<std::size_t>(t);
            if (k+1 == e) return tr.frames[e-1].hz;
            return tr.frames[k].hz + (t-k) * (tr.frames[k+1].hz-tr.frames[k].hz);
        };
        x.push_back(begin); rate.push_back(hz(begin) / sr);
        for (auto i = b; i < e; ++i) {
            const double p = double(i) * tr.step;
            if (p <= begin || p >= end) continue;
            x.push_back(p); rate.push_back(tr.frames[i].hz / sr);
        }
        x.push_back(end);
        rate.push_back(hz(end) / sr);
        phase.resize(x.size());
        for (std::size_t i = 1; i < x.size(); ++i)
            phase[i] = phase[i-1] + .5 * (x[i]-x[i-1]) * (rate[i]+rate[i-1]);
    }

    double period_at(double u) const {
        if (u <= x.front()) return 1. / rate.front();
        if (u >= x.back()) return 1. / rate.back();
        const auto k = std::size_t(std::upper_bound(x.begin(), x.end(), u) - x.begin() - 1);
        const double w = (u-x[k]) / (x[k+1]-x[k]);
        return 1. / (rate[k] + (rate[k+1]-rate[k])*w);
    }

    double at(double u) const {
        if (u <= x.front()) return 0.;
        if (u >= x.back()) return phase.back();
        const auto k = std::size_t(std::upper_bound(x.begin(), x.end(), u) - x.begin() - 1);
        const double d = u - x[k];
        const double slope = (rate[k+1]-rate[k]) / (x[k+1]-x[k]);
        return phase[k] + d * (rate[k] + .5*slope*d);
    }

    double time(double ph) const {
        if (ph <= 0.) return x.front();
        if (ph >= phase.back()) return x.back();
        const auto k = std::size_t(std::upper_bound(phase.begin(), phase.end(), ph) - phase.begin() - 1);
        const double delta = ph - phase[k];
        const double slope = (rate[k+1]-rate[k]) / (x[k+1]-x[k]);
        const double denom = rate[k] + std::sqrt(std::max(0., rate[k]*rate[k] + 2.*slope*delta));
        return x[k] + 2.*delta/denom;
    }
};

// Resampling and grain placement are sensitive to sudden phase changes.
// Follow the F0 prediction, correcting only when two neighboring bandlimited
// waveform neighborhoods match; raw high-frequency peaks are never chased.
class MarkLock {
    const Reduced& s_;
    const std::vector<float>& y_;
    std::vector<double> reference_;

    static double cubic(double a, double b, double c, double d, double u) {
        return b + .5*u*(c-a + u*(2.*a-5.*b+4.*c-d + u*(3.*(b-c)+d-a)));
    }

    // The samples in a trial window are unit spaced. Compute the tap
    // indices once and only do per-sample clipping at a recording edge.
    template<class Consume>
    void read_window(double p, int half, Consume consume) const {
        const int k = static_cast<int>(p), first = k-half, last = k+half;
        const int n = static_cast<int>(y_.size());
        const double u = p-k;
        // If the fractional offset is exact at the largest coordinate,
        // its binary precision is sufficient at all smaller coordinates.
        if (first >= 1 && last+2 < n && p+half-last == u) {
            const float* taps = y_.data() + first-1;
            for (int j = 0; j <= 2*half; ++j)
                consume(j, cubic(taps[j], taps[j+1], taps[j+2], taps[j+3], u));
        } else {
            for (int j = 0; j <= 2*half; ++j) {
                const double at = std::clamp(p+(j-half), 0., double(n-1));
                const int i = static_cast<int>(at);
                consume(j, cubic(y_[std::max(0,i-1)], y_[i], y_[std::min(n-1,i+1)], y_[std::min(n-1,i+2)], at-i));
            }
        }
    }

public:
    explicit MarkLock(const Reduced& s) : s_(s), y_(s.y) { reference_.reserve(300); }

    double seed(double center, double period, double begin, double end) const {
        const double p = center/s_.dec, r = .54 * period/s_.dec;
        const int first = std::max(0, static_cast<int>(std::ceil(std::max(begin/s_.dec,p-r))));
        const int last = std::min(static_cast<int>(y_.size())-1,
                                  static_cast<int>(std::floor(std::min(end/s_.dec,p+r))));
        if (last <= first) return center;
        double mean = 0.;
        for (int i = first; i <= last; ++i) mean += y_[i];
        mean /= last-first+1;
        double highest = 0.;
        int best = std::clamp(static_cast<int>(std::round(p)), first, last);
        for (int i = first; i <= last; ++i) {
            const double v = std::abs(y_[i] - mean) * (1. - .06*std::min(1.,std::abs(i-p)/r));
            if (v > highest) { best = i; highest = v; }
        }
        if (highest < 1e-7) return center;
        double at = best;
        if (best > first && best < last) {
            const double a = std::abs(y_[best-1]-mean), b = std::abs(y_[best]-mean), c = std::abs(y_[best+1]-mean);
            if (2.*b-a-c > 1e-12) at += std::clamp(.5*(c-a)/(2.*b-a-c), -.5, .5);
        }
        return std::clamp(at*s_.dec, begin, std::nextafter(end, begin));
    }

    double align(double prev, double predicted, double begin, double end) {
        const double a = prev/s_.dec, b = predicted/s_.dec, period = std::abs(b-a);
        const int half = std::max(2, static_cast<int>(std::round(.41*period)));
        const double radius = .15*period;
        if (std::min(a, b-radius)-half < begin/s_.dec || std::max(a, b+radius)+half >= end/s_.dec)
            return predicted;
        const int n = 2*half+1;
        reference_.resize(n);
        double m0 = 0., e0 = 0.;
        read_window(a, half, [&](int j, double v) { m0 += reference_[j] = v; });
        m0 /= n;
        for (double& v : reference_) { v -= m0; e0 += v*v; }
        if (e0 < n*1e-13) return predicted;
        auto similarity = [&](double f) {
            double dot = 0., e = 0., m = 0.;
            read_window(b+f*radius, half, [&](int j, double v) {
                dot += reference_[j]*v; m += v; e += v*v;
            });
            e -= m*m/n;
            return e > n*1e-13 ? std::clamp(dot / std::sqrt(e*e0), -1., 1.) : 0.;
        };
        constexpr int k = 4;
        std::array<double, 2*k+1> values{};
        int best = k;
        for (int i = 0; i <= 2*k; ++i) {
            const double f = (i-k)/double(k);
            values[i] = similarity(f) - .08*f*f;
        }
        for (int i = 0; i <= 2*k; ++i) if (values[i] > values[best]) best = i;
        double f = (best-k)/double(k);
        if (best > 0 && best < 2*k) {
            const double a0 = values[best-1], a1 = values[best], a2 = values[best+1];
            if (2.*a1-a0-a2 > 1e-10) f += std::clamp(.5*(a2-a0)/(2.*a1-a0-a2), -.5, .5) / k;
        }
        const double quality = values[best] + .08*std::pow((best-k)/double(k), 2);
        if (values[best] < values[k] + .001 || quality < .65) return predicted;
        const double trust = std::clamp((quality - .62) / .32, 0., 1.);
        return predicted + trust*f*radius*s_.dec;
    }
};

}  // namespace

Frame sample(const Track& tr, double source_sample) {
    if (tr.step <= 0 || tr.frames.empty() || std::isnan(source_sample)) return {};
    const double t = source_sample / tr.step;
    if (t <= 0.) return clean(tr.frames.front());
    if (t >= double(tr.frames.size()-1)) return clean(tr.frames.back());
    const auto i = static_cast<std::size_t>(t);
    const double w = t - double(i);
    const Frame a = clean(tr.frames[i]), b = clean(tr.frames[i+1]);
    if ((a.hz == 0.) != (b.hz == 0.)) return w < .5 ? a : b;
    return {a.hz + (b.hz-a.hz)*w, a.confidence + (b.confidence-a.confidence)*w};
}

Track estimate(const std::vector<float>& x, int sr) {
    Track tr;
    if (!valid_rate(sr) || x.empty() || x.size() >= std::size_t(std::numeric_limits<int>::max()-4096)) return tr;
    tr.step = (sr+100)/200;
    tr.frames.resize(1 + (x.size()-1) / tr.step);
    // Only F0 estimation uses this lower rate. Pulse alignment and the
    // harmonic separator still work with their original bandwidths.
    const Reduced r(x, sr, 6300, 2500., true);
    if (r.peak <= 1e-7) return tr;
    const Sums sums(r.y);
    const double level = std::max(1e-7, r.peak*.0018);
    const double floor = level*level;
    WindowedPitch window(r);
    std::vector<LatticeFrame> lattice(tr.frames.size());
    std::array<double, hypotheses> cost{};
    for (std::size_t i = 0; i < lattice.size(); ++i) {
        auto& current = lattice[i];
        const auto center = static_cast<int>((i*tr.step + r.dec/2) / r.dec);
        window.score(current, center, sums, r.rate, floor);
        std::array<double, hypotheses> next{};
        if (i == 0) {
            for (int j = 1; j < current.count; ++j) next[j] = current.h[j].reward - .06;
        } else {
            const auto& previous = lattice[i-1];
            for (int j = 0; j < current.count; ++j) {
                int best = 0;
                double value = cost[0] - transition(previous.h[0], current.h[j]);
                for (int k = 1; k < previous.count; ++k) {
                    const double v = cost[k] - transition(previous.h[k], current.h[j]);
                    if (v > value) { value = v; best = k; }
                }
                current.prev[j] = static_cast<unsigned char>(best);
                next[j] = value + current.h[j].reward;
            }
        }
        cost = next;
    }
    int best = 0;
    const auto& end = lattice.back();
    for (int j = 1; j < end.count; ++j) {
        if (cost[j] - .06 > cost[best] - (best ? .06 : 0.)) best = j;
    }
    for (std::size_t i = lattice.size(); i-- > 0;) {
        const auto& h = lattice[i].h[best];
        tr.frames[i] = {h.hz, h.strength};
        best = lattice[i].prev[best];
    }
    return tr;
}

Track candidate_clock(const std::vector<float>& x, int sr, const Track& anchors, bool supplied) {
    Track tr = anchors;
    if (tr.frames.empty()) return tr;
    const auto use = [&](Frame f) { return f.hz > 0 && (supplied || f.confidence >= .82); };
    if (std::all_of(tr.frames.begin(), tr.frames.end(), use)) return tr;

    // The normal tracker prefers the dominant broadband periodicity. Recheck
    // weak/zero frames on a low-band view so a fricative above a quiet pulse
    // train does not erase the trial period. This is not a second UV gate.
    const Reduced low(x, sr, 4200, 930.);
    const Sums sums(low.y);
    const int lower = std::max(3, static_cast<int>(std::floor(low.rate / 860.)));
    const int upper = static_cast<int>(std::ceil(low.rate / fmin)) + 1;
    const double floor = std::max(1e-14, low.peak * low.peak * 2e-7);
    const auto n = tr.frames.size();
    std::vector<double> upcoming(n);
    double next = 0., prev = 0.;
    for (std::size_t i = n; i-- > 0;) {
        if (use(anchors.frames[i])) next = anchors.frames[i].hz;
        upcoming[i] = next;
    }
    std::vector<float> lags(static_cast<std::size_t>(upper + 1));
    for (std::size_t i = 0; i < n; ++i) {
        auto& f = tr.frames[i];
        if (use(f)) { prev = f.hz; continue; }
        const int center = static_cast<int>((i * tr.step + low.dec / 2) / low.dec);
        const int a = std::clamp(center - 16, 0, static_cast<int>(low.y.size()));
        const int b = std::clamp(center + 17, a, static_cast<int>(low.y.size()));
        if (sums.power(a, b - a) <= floor) { f = {}; continue; }
        for (int lag = lower; lag <= upper; ++lag)
            lags[lag] = sums.nccf(center, lag, std::max(20, lag * 9 / 10), floor);

        // The priors disambiguate equivalent multiples of a period, not its
        // presence. Long gaps can also supply their own low-band F0.
        const double guide = f.hz > 0 ? f.hz : prev > 0 ? prev : upcoming[i];
        double best = -1.;
        Frame trial{};
        for (int lag = lower + 1; lag < upper; ++lag) {
            const double l = lags[lag-1], m = lags[lag], r = lags[lag+1];
            if (m < .60 || m <= l || m < r) continue;
            const auto [off, q] = refine(l, m, r);
            const double hz = low.rate / (lag + off);
            const double delta = guide > 0 ? std::abs(std::log2(hz / guide)) : 0.;
            const double rank = q - std::min(.23, .28 * delta) + .014 * std::log2(hz / fmin);
            const double required = guide > 0 ? .70 + .19 * std::clamp((delta - .18) / .70, 0., 1.) : .73;
            if (rank > best && q >= required && hz >= fmin * .985) {
                trial = {std::max(fmin, hz), q};
                best = rank;
            }
        }
        if (f.hz > 0 && f.confidence >= .65 &&
            (trial.hz == 0. || (f.confidence + .06 > trial.confidence &&
                               std::abs(std::log2(trial.hz / f.hz)) > .14))) {
            prev = f.hz;
        } else {
            f = trial;
            if (f.hz > 0) prev = f.hz;
        }
    }
    // Bridge unsupported spans in log Hz, without claiming they are voiced.
    // A neutral clock gives even a fully unlabelled clip a trial hypothesis.
    for (std::size_t i = 0; i < n;) {
        if (tr.frames[i].hz > 0) { ++i; continue; }
        const std::size_t begin = i;
        while (i < n && tr.frames[i].hz == 0) ++i;
        double left = begin ? tr.frames[begin-1].hz : (i < n ? tr.frames[i].hz : 150.);
        const double right = i < n ? tr.frames[i].hz : left;
        const double span = static_cast<double>(i - begin + 1);
        for (std::size_t j = begin; j < i; ++j)
            tr.frames[j].hz = left * std::pow(right / left, (j - begin + 1) / span);
    }
    return tr;
}

std::vector<Island> locate(const std::vector<float>& x, int sr, const Track& tr) {
    std::vector<Island> result;
    if (!valid_rate(sr) || x.empty() || tr.step <= 0 || tr.frames.empty() ||
        x.size() >= std::size_t(std::numeric_limits<int>::max()-4096)) return result;
    const double max_hz = std::min(2400., sr/3.);
    const auto voiced = [&](std::size_t f) {
        const double hz = tr.frames[f].hz;
        return hz >= 30. && hz <= max_hz;
    };
    if (std::none_of(tr.frames.begin(), tr.frames.end(), [&](Frame f){return f.hz >= 30. && f.hz <= max_hz;})) return result;
    const Reduced r(x, sr);
    MarkLock lock(r);
    const auto n = tr.frames.size();
    const double length = double(x.size());
    std::vector<double> backward;
    for (std::size_t b = 0, e = 0; b < n; b = e) {
        if (!voiced(b)) { e = b+1; continue; }
        e = b+1;
        while (e < n && voiced(e)) ++e;
        const double begin = b == 0 ? 0. : (double(b) - .5)*tr.step;
        if (begin >= length) break;
        const double end = e == n ? length : std::min(length, (double(e)-.5)*tr.step);
        if (end <= begin) continue;
        const PhaseCurve curve(tr, b, e, sr, begin, end);
        if (curve.phase.back() < 1.) continue;
        const double center = curve.time(.5*curve.phase.back());
        const double seed = lock.seed(center, curve.period_at(center), begin, end);
        Island island{begin, end, {}};
        backward.clear();
        const auto follow = [&](int dir, std::vector<double>& out) {
            double prev = seed;
            double ph = curve.at(prev);
            for (;;) {
                const double goal = ph + dir;
                if (goal < -1e-10 || goal >= curve.phase.back()) break;
                const double predict = curve.time(goal);
                double at = lock.align(prev, predict, begin, end);
                if (!(at >= begin && at < end)) at = predict;
                if (!(at >= begin && at < end && (at-prev)*dir > 0.)) break;
                out.push_back(at);
                prev = at;
                ph = curve.at(at);
            }
        };
        follow(-1, backward);
        island.marks.insert(island.marks.end(), backward.rbegin(), backward.rend());
        island.marks.push_back(seed);
        follow(1, island.marks);
        if (island.marks.size() >= 2) result.push_back(std::move(island));
    }
    return result;
}

}  // namespace tn::pitch
