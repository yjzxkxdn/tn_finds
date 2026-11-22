// SPDX-License-Identifier: GPL-3.0-only
#include "formant.hpp"
#include "dsp_fft.hpp"
#include "engine.hpp"
#include "glottal.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace tn::detail {
namespace {
constexpr double pi = 3.14159265358979323846;

void smooth_power(std::vector<float>& power, std::vector<float>& tmp, int half) {
    const int last = static_cast<int>(power.size()) - 1;
    const auto reflected = [&](int k) { return power[k < 0 ? -k : k > last ? 2*last-k : k]; };
    const double scale = 1. / (2 * half + 1.);
    double sum = 0.;
    for (int k = -half; k <= half; ++k) sum += reflected(k);
    tmp[0] = static_cast<float>(sum * scale);
    for (int k = 1; k <= last; ++k) {
        sum += static_cast<double>(reflected(k + half)) - reflected(k - half - 1);
        tmp[k] = static_cast<float>(sum * scale);
    }
    power.swap(tmp);
}
} // namespace

void shift_formant_envelope(std::vector<float>& x, int sr, double ratio, const pitch::Track& track,
                            double min_smoothing_hz) {
    if (x.empty() || ratio == 1.) return;
    const int n = std::clamp(1 << static_cast<int>(std::lround(std::log2(sr * .044))), 512, 8192);
    const int last = n / 2, hop = n / 4, len = static_cast<int>(x.size());
    const double df = static_cast<double>(sr) / n;
    const Fft fft(n);
    std::vector<Complex> spec(n);
    std::vector<float> win(n), envelope(last + 1), tmp(last + 1), correction(len);
    for (int i = 0; i < n; ++i) win[i] = static_cast<float>(std::sin(pi * i / n));

    // At 75% overlap the sum of squared sqrt-Hann windows is 2.
    for (int center = -hop; center - last < len; center += hop) {
        const int start = center - last;
        for (int i = 0; i < n; ++i) {
            const int k = start + i;
            spec[i] = k >= 0 && k < len ? Complex(x[k] * win[i], 0.f) : Complex();
        }
        fft.run(spec);
        for (int k = 0; k <= last; ++k) envelope[k] = std::norm(spec[k]);

        // Follow individual formants while averaging out harmonic lines or
        // the fine spectral fluctuations of the aperiodic residual.
        const double half_width = std::max(.5 * pitch::sample(track, center).hz, min_smoothing_hz);
        const int width = std::clamp(static_cast<int>(std::lround(half_width / df)), 1, last - 1);
        smooth_power(envelope, tmp, width);
        smooth_power(envelope, tmp, width);
        const float floor = std::max(1.e-20f, *std::max_element(envelope.begin(), envelope.end()) * 1.e-9f);
        for (auto& power : envelope) power = std::sqrt(std::max(power, floor));

        // Move only the envelope, E_new(f) = E_old(f / ratio). Multiplying
        // by the real ratio of envelopes keeps each FFT bin on its own grid.
        spec[0] = Complex();
        for (int k = 1; k <= last; ++k) {
            const double at = k / ratio;
            const int j = std::min(last - 1, static_cast<int>(at));
            const float frac = static_cast<float>(std::min(1., at - j));
            const float mapped = envelope[j] + frac * (envelope[j+1] - envelope[j]);
            // Fade only the part that maps above the source Nyquist limit.
            const float edge = static_cast<float>(std::clamp(1. + (last - at) / width, 0., 1.));
            spec[k] *= mapped * edge / envelope[k] - 1.f;
            if (k < last) spec[n-k] = std::conj(spec[k]);
        }
        fft.run(spec, true);
        const float norm = .5f / n;
        for (int i = std::max(0, -start); i < std::min(n, len - start); ++i)
            correction[start + i] += norm * win[i] * spec[i].real();
    }
    for (int k = 0; k < len; ++k) x[k] += correction[k];
}
namespace {
int voice_fft_size(int sr) {
    return std::clamp(1 << static_cast<int>(std::lround(std::log2(sr * .044))), 512, 8192);
}

double mouth_log_gain(double hz, double f0) {
    // Log of the squared distance from a repeated pole near target H1.
    // The effective pole width follows Moresampler 0.8.4; its frequency axis
    // is independent of the output WAV rate.
    constexpr double axis_rate = 49226., bandwidth = 153.;
    static const double radius = std::exp(-2.*pi*bandwidth/axis_rate);
    static const double scale = 4.*radius/((1.-radius)*(1.-radius));
    const double angle = pi*std::min(std::abs(hz-f0), axis_rate*.5)/axis_rate;
    const double sine = std::sin(angle);
    return std::log1p(scale*sine*sine);
}

// Sqrt-Hann analysis/synthesis at 75% overlap has squared-window sum 2.
template<class Edit>
std::vector<float> voice_windows(const std::vector<float>& x, int n, Edit edit) {
    const int last = n/2, hop = n/4, len = static_cast<int>(x.size());
    const float norm = .5f/n;
    const Fft fft(n);
    std::vector<Complex> spec(n);
    std::vector<float> win(n), out(len);
    for (int i = 0; i < n; ++i) win[i] = static_cast<float>(std::sin(pi*i/n));
    for (int center = -hop; center-last < len; center += hop) {
        const int start = center-last;
        for (int i = 0; i < n; ++i) {
            const int index = start+i;
            spec[i] = index >= 0 && index < len ? Complex(x[index]*win[i], 0.f) : Complex();
        }
        fft.run(spec);
        edit(spec, center);
        for (int k = 1; k < last; ++k) spec[n-k] = std::conj(spec[k]);
        fft.run(spec, true);
        for (int i = std::max(0, -start); i < std::min(n, len-start); ++i)
            out[start+i] += norm*win[i]*spec[i].real();
    }
    return out;
}

int half_width(double hz, double df, int last) {
    return std::clamp(static_cast<int>(std::lround(hz/df)), 1, last-1);
}

void envelope_power(const std::vector<Complex>& spec, std::vector<float>& power, std::vector<float>& tmp,
                    int width) {
    for (std::size_t k = 0; k < power.size(); ++k) power[k] = std::norm(spec[k]);
    smooth_power(power, tmp, width);
    smooth_power(power, tmp, width);
}

float power_at(const std::vector<float>& power, double at) {
    const int i = std::min(static_cast<int>(power.size())-2, static_cast<int>(at));
    const float u = static_cast<float>(std::min(1., at-i));
    return power[i] + u*(power[i+1]-power[i]);
}
} // namespace

void shape_voice_envelope(std::vector<float>& x, int sr, const pitch::Track& track,
                          const Options& o, bool harmonic, const GlottalTrack* glottal) {
    const double tension = harmonic && glottal && !glottal->rd.empty() ? o.tension : 0.;
    const double mouth = o.mouth/300., ring = o.resonance*.01, contrast = o.formant_emphasis*.01;
    if (x.empty() || (tension == 0. && mouth == 0. && ring == 0. && contrast == 0.)) return;
    const int n = voice_fft_size(sr), last = n/2;
    const double df = double(sr)/n;
    constexpr double db_per_octave = 6.020599913279624;
    constexpr double db_per_neper = 8.685889638065037;
    const bool adaptive = tension != 0. || mouth != 0. || contrast != 0.;
    std::vector<double> db(last+1), lf;
    std::vector<float> gain(last+1), narrow(last+1), broad(last+1), tmp(last+1);
    for (int k = 1; k <= last; ++k) {
        const double r = (k*df-3000.)/405.;
        db[k] = ring*10.*std::exp(-.5*r*r);
        gain[k] = static_cast<float>(std::exp2(db[k]/db_per_octave)-1.);
    }
    const auto correction = voice_windows(x, n, [&](auto& spec, int center) {
        if (adaptive) {
            const double hz = pitch::sample(track, center).hz;
            double floor = 0.;
            if (contrast != 0.) {
                const double width = std::max(harmonic ? 45. : 205., .5*hz);
                envelope_power(spec, narrow, tmp, half_width(width, df, last));
                floor = std::max(1.e-20f, *std::max_element(narrow.begin(), narrow.end())*1.e-8f);
                broad = narrow;
                const int wider = half_width(std::max(480., 1.65*width), df, last);
                smooth_power(broad, tmp, wider);
                smooth_power(broad, tmp, wider);
            }
            if (tension != 0. && hz > 0.) {
                const int orders = static_cast<int>(std::ceil(std::min(sr*.5, 24613.)/hz));
                tension_gain(glottal->sample(center), tension, orders, lf);
            } else {
                lf.clear();
            }
            for (int k = 1; k <= last; ++k) {
                double add = db[k];
                if (contrast != 0.) {
                    const double power = std::max(floor, double(power_at(narrow, k)));
                    const double reference = std::max(floor, double(power_at(broad, k)));
                    add += contrast*.85*std::clamp(10.*std::log10(power/reference), -9., 9.);
                    add = std::clamp(add, -18., 18.);
                }
                double color = 0.;
                if (mouth != 0. && hz > 0.) {
                    const double opening = mouth*mouth_log_gain(k*df, hz);
                    color += harmonic ? opening : std::min(1., opening);
                }
                if (lf.size() > 2) {
                    const double at = std::clamp(k*df/hz, 1., double(lf.size()-1));
                    const int j = std::min(static_cast<int>(at), static_cast<int>(lf.size())-2);
                    color += lf[j]+(at-j)*(lf[j+1]-lf[j]);
                }
                gain[k] = static_cast<float>(std::exp2((add+db_per_neper*color)/db_per_octave)-1.);
            }
        }
        spec[0] = Complex();
        for (int k = 1; k <= last; ++k) spec[k] *= gain[k];
    });
    for (std::size_t i = 0; i < x.size(); ++i) x[i] += correction[i];
}
void color_harmonics(HarmonicFrame& frame, int sr, const Options& o, double rd,
                     double voice_mix, double level_hz) {
    const auto count = static_cast<int>(frame.harmonics.size());
    if (count == 0) return;
    if (o.tension == 0. && o.mouth == 0. && o.resonance == 0. && o.formant_emphasis == 0. &&
        o.opening == 0. && o.dynamics == 0. && o.breathiness <= 0.) return;
    const double hz = frame.hz;
    std::vector<double> lf;
    if (rd > 0. && o.tension != 0.) tension_gain(rd, o.tension, count, lf);
    std::vector<float> narrow, broad;
    if (o.formant_emphasis != 0. && count >= 3) {
        narrow.resize(count+1);
        for (int k = 1; k <= count; ++k) narrow[k] = std::norm(frame.harmonics[k-1]);
        std::vector<float> tmp(narrow.size());
        smooth_power(narrow, tmp, 1);
        smooth_power(narrow, tmp, 1);
        broad = narrow;
        const int wider = half_width(std::max(480., .825*hz), hz, count);
        smooth_power(broad, tmp, wider);
        smooth_power(broad, tmp, wider);
    }
    const double shelf = std::exp2(o.opening/125.)-1.;
    const double pole = std::exp(-2*pi*1700./sr);
    const double dynamics = level_hz > 0. ? 1. + o.dynamics*.0012*
        std::clamp(std::log2(hz/level_hz)*4., -1., 1.) : 1.;
    const double whisper = 1.-std::max(0., o.breathiness*.01)*voice_mix;
    for (int k = 1; k <= count; ++k) {
        const double f = k*hz;
        const double r = (f-3000.)/405.;
        double db = o.resonance*.1*std::exp(-.5*r*r);
        if (!broad.empty()) {
            constexpr double floor = 1.e-20;
            db += o.formant_emphasis*.0085*std::clamp(
                10.*std::log10(std::max(floor, double(narrow[k]))/std::max(floor, double(broad[k]))), -9., 9.);
            db = std::clamp(db, -18., 18.);
        }
        const double color = o.mouth/300.*mouth_log_gain(f, hz) +
                             (k < static_cast<int>(lf.size()) ? lf[k] : 0.);
        double gain = 1. + voice_mix*(std::exp(db*std::log(10.)/20. + color)-1.);
        if (shelf != 0.) {
            const auto lp = (1.-pole)/(1.-pole*std::polar(1., -2*pi*f/sr));
            gain += std::abs(1.+shelf*(1.-lp))-1.;
        }
        frame.harmonics[k-1] *= static_cast<float>(std::max(0., gain)*dynamics*whisper);
    }
}
} // namespace tn::detail
