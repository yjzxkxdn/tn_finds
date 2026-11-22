// SPDX-License-Identifier: GPL-3.0-only
#include "harmonic_phase.hpp"
#include "glottal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace tn::detail {
namespace {
constexpr double pi = 3.14159265358979323846;

float ease(double t) {
    t = std::clamp(t, 0., 1.);
    return static_cast<float>(t*t*(3-2*t));
}

float interpolate(const std::vector<float>& a, double pos) {
    pos = std::clamp(pos, 0., static_cast<double>(a.size() - 1));
    const int k = static_cast<int>(pos);
    const float t = static_cast<float>(pos - k);
    const float b = a[k], c = a[std::min(k + 1, static_cast<int>(a.size()) - 1)];
    if (k < 1 || k + 2 >= static_cast<int>(a.size())) return b + (c-b)*t;
    const float l = a[k-1], r = a[k+2];
    return b + .5f*t*(c-l + t*(2*l-5*b+4*c-r + t*(3*(b-c)+r-l)));
}

// A mark-centered cycle avoids choosing a splice at the two adjacent marks.
// Phase 0.5 is this epoch; ±0.5 cycles are from the left/right neighbor.
float around_epoch(const std::vector<float>& x, const std::vector<double>& p,
                   int i, double phase, int origin) {
    const auto mid = p[i];
    const auto side = p[i + (phase >= 0 ? 1 : -1)];
    return interpolate(x, mid + std::abs(phase) * (side - mid) - origin);
}

struct Anchor {
    int mark;
    std::vector<float> correction;
};

// A fixed offset of the source mark within the glottal cycle must not
// influence the comparison: remove a common time-origin slope k * phi1.
// Weight phase transfer only where both source phasor and wet harmonic are
// supported, and preserve the high band where phase estimates are noisy.
Anchor make_anchor(const std::vector<float>& source, const std::vector<double>& src,
                   const std::vector<float>& wet, const std::vector<double>& dst,
                   double p, int i, int origin, int sr, const Fft& fft,
                   std::vector<Complex>& a, std::vector<Complex>& b) {
    const int length = static_cast<int>(a.size());
    const auto it = std::upper_bound(src.begin(), src.end(), p);
    const int right = std::clamp(static_cast<int>(it - src.begin()), 2, static_cast<int>(src.size()) - 2);
    const int left = right - 1;
    const float w = ease((p - src[left]) / (src[right] - src[left]));
    auto source_cycle = [&](int j, double phase) {
        float sum = 0;
        int weight = 0;
        for (int k = -1; k <= 1; ++k) {
            if (j+k < 1 || j+k+1 >= static_cast<int>(src.size())) continue;
            const int m = k == 0 ? 2 : 1;
            sum += m * around_epoch(source, src, j+k, phase, 0);
            weight += m;
        }
        return sum / weight;
    };
    for (int q = 0; q < length; ++q) {
        const double phase = static_cast<double>(q) / length - .5;
        a[q] = source_cycle(left, phase) * (1-w) + source_cycle(right, phase) * w;
        b[q] = around_epoch(wet, dst, i, phase, origin);
    }
    fft.run(a);
    fft.run(b);
    const int half = length / 2;
    const double f0 = 2. * sr / (dst[i+1] - dst[i-1]);
    const int last = std::max(1, std::min(half - 1, static_cast<int>(std::ceil(std::min(4800., .44 * sr) / f0))));
    float amax = 0, bmax = 0;
    for (int k = 1; k <= last; ++k) {
        amax = std::max(amax, std::abs(a[k]));
        bmax = std::max(bmax, std::abs(b[k]));
    }
    // If H1 is genuinely absent its phase is undefined; this anchor merely
    // stays dry. It cannot inject a broadband random linear-phase slope.
    const float reference = ease((std::abs(a[1]) / std::max(1.e-8f, amax) - .006) / .03) *
                            ease((std::abs(b[1]) / std::max(1.e-8f, bmax) - .006) / .03);
    Anchor out{i, std::vector<float>(length)};
    if (reference < 1.e-6f || amax < 1.e-8f || bmax < 1.e-8f) return out;
    const double phi1 = std::arg(a[1] * std::conj(b[1]));
    b[0] = b[1] = b[half] = b[length-1] = Complex();
    for (int k = 2; k < half; ++k) {
        Complex delta{};
        if (k <= last) {
            const float aa = std::abs(a[k]) / amax, bb = std::abs(b[k]) / bmax;
            const float matched = ease((aa - .008) / .060) * ease((bb - .005) / .040);
            const float band = 1.f - ease((k * f0 - 2800.) / (std::min(4800., .44 * sr) - 2800.));
            const double shift = std::remainder(std::arg(a[k] * std::conj(b[k])) - k*phi1, 2*pi);
            // Rotate wet coefficient without changing its magnitude at the
            // anchor. Corrected and unchanged share the same fundamental.
            const double phi = shift * matched * band * reference;
            delta = b[k] * (Complex(static_cast<float>(std::cos(phi)), static_cast<float>(std::sin(phi))) - Complex(1.f, 0.f));
        }
        b[k] = delta;
        b[length-k] = std::conj(delta);
    }
    fft.run(b, true);
    for (int q = 0; q < length; ++q) out.correction[q] = b[q].real() / static_cast<float>(length);
    return out;
}

float cyclic(const std::vector<float>& table, double phase) {
    const int size = static_cast<int>(table.size());
    const double p = (phase + .5) * size;
    const int k = static_cast<int>(p);
    const int mask = size - 1;
    const float b = table[k & mask], c = table[(k+1) & mask];
    return b + static_cast<float>(p-k) * (c-b);
}
double harmonic_value(const std::vector<double>& v, double order) {
    const int last = static_cast<int>(v.size())-1;
    const double at = std::clamp(order, 1., static_cast<double>(last));
    const int i = static_cast<int>(at), j = std::min(i+1, last);
    return v[i]+(at-i)*(v[j]-v[i]);
}

void tract_minimum_phase(const std::vector<double>& magnitude, std::vector<double>& phase) {
    const int orders = static_cast<int>(magnitude.size())-1;
    int n = 64, bits = 6;
    while (n < 4*orders) { n *= 2; ++bits; }
    // Target orders are bounded by 12 kHz / 30 Hz. Reuse the small FFT
    // plans and buffer across frames, independently for each render thread.
    static thread_local std::array<std::optional<Fft>,6> plans;
    static thread_local std::vector<Complex> cepstrum;
    auto& fft = plans[bits-6];
    if (!fft) fft.emplace(n);
    cepstrum.resize(n);
    const double spacing = double(n/2)/(orders+1);
    for (int k = 0; k <= n/2; ++k) {
        const auto level = static_cast<float>(harmonic_value(magnitude, k/spacing));
        cepstrum[k] = level;
        if (k > 0 && k < n/2) cepstrum[n-k] = level;
    }
    fft->run(cepstrum, true);
    for (int i = 1; i < n/2; ++i) cepstrum[i] *= 2.f;
    std::fill(cepstrum.begin()+n/2+1, cepstrum.end(), Complex());
    fft->run(cepstrum);
    phase.resize(orders+1);
    phase[0] = 0.;
    for (int k = 1; k <= orders; ++k) {
        const double at = k*spacing;
        const int i = static_cast<int>(at);
        phase[k] = (cepstrum[i].imag()+(at-i)*(cepstrum[i+1].imag()-cepstrum[i].imag()))/n;
    }
}

// Short-time radiation filter for a 1.5-cm mouth radius. The linear
// frequency term in its magnitude cancels integration of the LF derivative.
double radiation_argument(double hz) {
    constexpr double resistance = 128./(9.*pi*pi);
    constexpr double inductance = 8.*.015/(3.*pi*340.);
    return 2.*pi*hz*inductance/resistance;
}

double radiation_phase(double hz) {
    return pi*.5-std::atan(radiation_argument(hz));
}
} // namespace

void source_filter_spectrum(const std::vector<Complex>& spectrum, int source_orders,
                            double source_hz, double target_hz, double rd, double formant_ratio,
                            std::vector<double>& phase, std::vector<float>* amplitude) {
    phase.clear();
    if (amplitude) amplitude->clear();
    if (rd <= 0. || source_orders < 4) return;
    const int orders = static_cast<int>(std::ceil(source_orders*source_hz/target_hz));
    if (orders < 2) return;
    float peak = 0.f;
    for (int k = 1; k <= source_orders; ++k) peak = std::max(peak, std::abs(spectrum[k]));
    if (peak <= 1.e-10f) return;

    std::vector<std::complex<double>> lf;
    lf_spectrum(rd, std::max(source_orders, orders), lf);
    std::vector<double> tract(source_orders+1), target_tract(orders+1);
    for (int k = 1; k <= source_orders; ++k) {
        const double measured = std::max(double(std::abs(spectrum[k])/peak), 1.e-5);
        const double input_lip = radiation_argument(k*source_hz);
        tract[k] = std::log(measured) - std::log(std::max(1.e-14, std::abs(lf[k]))) +
                   .5*std::log1p(input_lip*input_lip);
    }
    const int amplitude_orders = amplitude ? std::min(orders, static_cast<int>(8000./target_hz)) : 0;
    if (amplitude) amplitude->resize(amplitude_orders+1);
    for (int k = 1; k <= orders; ++k) {
        const double target_frequency = k*target_hz;
        double frequency = target_frequency;
        // On upward shifts Moresampler anchors the bottom of the tract
        // envelope at the new H1, tapering the frequency map to identity.
        if (target_hz > source_hz) {
            const double onset_warp = std::clamp((target_hz+500.-frequency)/500., 0., 1.);
            frequency -= onset_warp*(target_hz-source_hz);
        }
        target_tract[k] = harmonic_value(tract, frequency/(source_hz*formant_ratio));
        if (k <= amplitude_orders) {
            // Flow integration cancels the frequency term in lip radiation.
            // A common gain is handled by the synthesis table's power match.
            const double lip = radiation_argument(target_frequency);
            (*amplitude)[k] = static_cast<float>(peak*std::exp(target_tract[k])*
                                                std::abs(lf[k])/std::sqrt(1.+lip*lip));
        }
    }

    std::vector<double> source_minimum;
    tract_minimum_phase(tract, source_minimum);
    tract_minimum_phase(target_tract, phase);
    const double source_h1 = std::arg(spectrum[1])-source_minimum[1]-radiation_phase(source_hz);
    const double target_filter_h1 = phase[1]+radiation_phase(target_hz);
    for (int k = 2; k <= orders; ++k) {
        double excitation = std::arg(lf[k]);
        if (k <= source_orders) {
            excitation = std::arg(spectrum[k])-source_minimum[k]-radiation_phase(k*source_hz)-k*source_h1;
        }
        phase[k] = std::remainder(excitation+phase[k]+radiation_phase(k*target_hz)-k*target_filter_h1, 2.*pi);
    }
    phase[1] = 0.;
}

void lock_harmonic_phase(const std::vector<float>& source, const std::vector<double>& src,
                         const std::vector<double>& dst, const std::vector<double>& at,
                         int sr, int lo, std::vector<float>& wet) {
    if (src.size() < 4 || dst.size() < 5 || dst.size() != at.size() || wet.empty()) return;
    constexpr int length = 512;
    const Fft fft(length);
    std::vector<Complex> a(length), b(length);
    // Delay-origin reference is H1; interpolate successive corrections in
    // one shared target-cycle coordinate. Retain just two anchor tables. Read
    // from the uncorrected wave so the next FFT does not re-analyze a splice.
    const auto original = wet;
    const int nz = static_cast<int>(dst.size()), last = nz - 2;
    const auto next = [&](int i) {
        int j = i + 1;
        while (j < last && dst[j] - dst[i] < .008 * sr) ++j;
        return std::min(j, last);
    };
    const auto build = [&](int i) { return make_anchor(source, src, original, dst, at[i], i, lo, sr, fft, a, b); };
    Anchor l = build(1), r = build(next(1));
    int period = 0;
    const double fade = std::min(.012 * sr, .5 * (dst.back() - dst.front()));
    for (int j = 0; j < static_cast<int>(wet.size()); ++j) {
        const auto n = static_cast<double>(lo + j);
        while (period + 2 < nz && n >= dst[period+1]) ++period;
        while (r.mark < last && n > dst[r.mark]) {
            l = std::move(r);
            r = build(next(l.mark));
        }
        const double phase = (n - dst[period]) / (dst[period+1] - dst[period]);
        const double phi = std::clamp(phase, 0., 1.) - (phase >= .5 ? 1. : 0.);
        const float a0 = cyclic(l.correction, phi), a1 = cyclic(r.correction, phi);
        const float w = ease((n-dst[l.mark]) / (dst[r.mark] - dst[l.mark]));
        const float change = a0 + w * (a1 - a0);
        wet[j] += ease(std::min(n-dst.front(), dst.back()-n) / fade) * change;
    }
}
} // namespace tn::detail
