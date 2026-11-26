// SPDX-License-Identifier: GPL-3.0-only
#include "harmonic_synth.hpp"
#include "dsp_fft.hpp"
#include "glottal.hpp"
#include "harmonic_phase.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <type_traits>
#include <vector>

namespace tn::detail {
namespace {
constexpr double pi = 3.14159265358979323846;

double smooth(double x) {
    x = std::clamp(x, 0., 1.);
    return x*x*(3.-2.*x);
}

float cubic(float a, float b, float c, float d, float t) {
    return b + .5f*t*(c-a + t*(2*a-5*b+4*c-d + t*(3*(b-c)+d-a)));
}

// The cycle fill calls this twice per table sample; keeping the four-tap
// interpolation in line removes the call overhead from that inner loop.
#if defined(_MSC_VER)
#define TN_FORCE_INLINE __forceinline
#else
#define TN_FORCE_INLINE inline __attribute__((always_inline))
#endif

template<bool interior>
TN_FORCE_INLINE float sample(const std::vector<float>& x, double at) {
    if constexpr (interior) {
        const int j = static_cast<int>(at); // All four taps are inside x, at >= 0.
        return cubic(x[j-1], x[j], x[j+1], x[j+2], static_cast<float>(at-j));
    } else {
        const int j = static_cast<int>(std::floor(at));
        const auto read = [&](int i) { return i < 0 || i >= static_cast<int>(x.size()) ? 0.f : x[i]; };
        return cubic(read(j-1), read(j), read(j+1), read(j+2), static_cast<float>(at-j));
    }
}

bool interior_span(const std::vector<float>& x, double first, double last, double span) {
    // Leave an extra sample of margin for rounding the interpolated positions.
    return first-span >= 2. && last+span <= double(x.size())-3.;
}

struct Cycle {
    double center, left, right;
    double period() const { return .5*(left+right); }
};

Cycle cycle_at(const std::vector<double>& marks, int index) {
    const double center = marks[index];
    const double left = index > 0 ? center-marks[index-1] : marks[1]-center;
    const double right = index+1 < static_cast<int>(marks.size()) ? marks[index+1]-center : left;
    return {center, left, right};
}

float periodic(const std::vector<float>& cycle, double phase) {
    const int size = static_cast<int>(cycle.size()), mask = size-1;
    const double at = phase*size;
    const int j = static_cast<int>(at);
    return cubic(cycle[(j-1)&mask], cycle[j&mask], cycle[(j+1)&mask], cycle[(j+2)&mask],
                 static_cast<float>(at-j));
}

class SourceCycles {
    struct Cached {
        int index = -1;
        std::vector<float> wave;
    };
    std::array<Cached, 4> cache_;
    const std::vector<float>& source_;
    const std::vector<double>& marks_;
    const HarmonicSplit& split_;
    int rate_;
    int left_ = 0;
    double mix_ = 0.;

    const std::vector<float>& cycle(int index, int size) {
        auto& entry = cache_[index % cache_.size()];
        if (entry.index == index) return entry.wave;
        entry.index = index;
        entry.wave.resize(size);
        const Cycle c = cycle_at(marks_, index);
        const double period = c.period();
        const auto fill = [&](auto inside) {
            for (int j = 0; j < size; ++j) {
                const double u = double(j)/size, back = 1.-u;
                // Both sides meet at the same sample with the same phase slope.
                const double ahead = u*(c.right+(period-c.right)*back*back*back*(1.+3.*u));
                const double behind = back*(c.left+(period-c.left)*u*u*u*(1.+3.*back));
                const float a = sample<inside.value>(source_, c.center+ahead);
                const float b = sample<inside.value>(source_, c.center-behind);
                const float blend = static_cast<float>(u*u*u*(10.+u*(-15.+6.*u)));
                entry.wave[j] = a + blend*(b-a);
            }
        };
        if (interior_span(source_, c.center, c.center, std::max(c.left, c.right))) fill(std::true_type{});
        else fill(std::false_type{});
        return entry.wave;
    }

    double neighbor_weight(int index) const {
        if (index < 1 || index+1 >= static_cast<int>(marks_.size())) return 0.;
        const double a = marks_[index-1], b = marks_[index+1];
        const double voiced = smooth(std::min(split_.activity(a), split_.activity(b))/.25);
        // The symmetric 1:2:1 kernel rejects alternating cycles. Fade its
        // support near voice boundaries and avoid smearing very low notes.
        return .25*voiced*(1.-smooth(((b-a)/rate_-.012)/.012));
    }

public:
    SourceCycles(const std::vector<float>& source, const std::vector<double>& marks,
                 const HarmonicSplit& split, int rate) : source_(source), marks_(marks), split_(split), rate_(rate) {}

    double read(double at, std::vector<float>& out) {
        const int size = static_cast<int>(out.size());
        auto next = std::upper_bound(marks_.begin(), marks_.end(), at);
        const int left = std::clamp(static_cast<int>(next-marks_.begin())-1, 0, static_cast<int>(marks_.size())-2);
        const double mix = std::clamp((at-marks_[left])/(marks_[left+1]-marks_[left]), 0., 1.);
        left_ = left;
        mix_ = mix;
        const double a = neighbor_weight(left), b = neighbor_weight(left+1);
        const std::array<double,4> gains = {
            (1.-mix)*a, (1.-mix)*(1.-2.*a)+mix*b,
            mix*(1.-2.*b)+(1.-mix)*a, mix*b
        };
        double period = 0.;
        std::fill(out.begin(), out.end(), 0.f);
        for (int j = 0; j < 4; ++j) {
            if (gains[j] <= 0.) continue;
            const int index = left+j-1;
            period += gains[j]*cycle_at(marks_, index).period();
            const auto& wave = cycle(index, size);
            const auto gain = static_cast<float>(gains[j]);
            for (int k = 0; k < size; ++k) out[k] += gain*wave[k];
        }
        return rate_/period;
    }

    // Shift two successive, phase-aligned source pulses by one target period
    // and overlap their half cosine bells. Use only the harmonic recording.
    bool pulse(double period, const std::vector<float>& bell, std::vector<Complex>& out) const {
        const auto a = cycle_at(marks_, left_), b = cycle_at(marks_, left_+1);
        const double width = std::min({a.left, a.right, b.left, b.right});
        if (period > width || std::min(split_.activity(a.center), split_.activity(b.center)) < .35)
            return false;
        // With both source spans long enough the complementary windows sum
        // to one. Reject shorter spans instead of filling a gap in the probe.
        const double step = period/bell.size();
        const auto u = static_cast<float>(mix_);
        const auto fill = [&](auto inside) {
            for (std::size_t j = 0; j < bell.size(); ++j) {
                const double dt = j*step, previous = dt-period;
                const float a0 = sample<inside.value>(source_, a.center+dt), b0 = sample<inside.value>(source_, b.center+dt);
                const float a1 = sample<inside.value>(source_, a.center+previous), b1 = sample<inside.value>(source_, b.center+previous);
                const float forward = a0+u*(b0-a0), backward = a1+u*(b1-a1);
                out[j] = backward+bell[j]*(forward-backward);
            }
        };
        if (interior_span(source_, a.center, b.center, period)) fill(std::true_type{});
        else fill(std::false_type{});
        return true;
    }
};

// Two-times oversampled tables leave room for the interpolation transition
// band. Interpolated Lagrange weights suppress replicas that would otherwise
// fold into the spaces between target harmonics at fractional periods.
class TableReader {
    static constexpr int phases = 256;
    static constexpr int taps = 8;
    std::array<std::array<float, taps>, phases+1> weights_{};
public:
    TableReader() {
        for (int i = 0; i <= phases; ++i) {
            const double x = double(i)/phases + 3.;
            for (int j = 0; j < taps; ++j) {
                double w = 1.;
                for (int k = 0; k < taps; ++k) if (k != j) w *= (x-k)/(j-k);
                weights_[i][j] = static_cast<float>(w);
            }
        }
    }

    float read(const std::vector<float>& left, const std::vector<float>& right,
               double phase, float left_gain, float right_gain) const {
        const int size = static_cast<int>((left_gain > 0.f ? left : right).size()), mask = size-1;
        const double at = phase*size;
        const int i = static_cast<int>(at);
        const double f = (at-i)*phases;
        const int row = static_cast<int>(f);
        const auto& a = weights_[row];
        const auto& b = weights_[row+1];
        const float mix = static_cast<float>(f-row);
        std::array<float, taps> coeff;
        for (int j = 0; j < taps; ++j) coeff[j] = a[j]+mix*(b[j]-a[j]);
        const auto interpolate = [&](const std::vector<float>& wave) {
            float value = 0.f;
            for (int j = 0; j < taps; ++j) value += wave[(i+j-3)&mask]*coeff[j];
            return value;
        };
        return (left_gain > 0.f ? left_gain*interpolate(left) : 0.f) +
               (right_gain > 0.f ? right_gain*interpolate(right) : 0.f);
    }
};

struct Table {
    std::vector<float> wave;
    float weight = 0.f;
};

constexpr double no_phase = std::numeric_limits<double>::quiet_NaN();

class Tables {
    Fft fft_;
    std::optional<RealInverseFft> output_fft_;
    std::vector<Complex> spec_, folded_, output_spec_, unshaped_spec_;
    std::vector<float> cycle_, original_cycle_, bell_, source_power_, source_angle_, log_amp_, linear_amp_,
                       phase_, model_amplitude_;
    std::vector<double> arcs_, output_phase_, fold_difference_, model_phase_;
    SourceFilterWork filter_work_;
    double spectrum_at_ = no_phase, spectrum_hz_ = 0., folded_hz_ = no_phase, old_power_ = 0.;
    float peak_ = 0.f;
    double source_origin_ = 0., output_origin_ = 0.;
    bool voiced_ = false;

    void reset_phase() {
        if (!voiced_) return;
        std::fill(arcs_.begin(), arcs_.end(), no_phase);
        std::fill(output_phase_.begin(), output_phase_.end(), no_phase);
        std::fill(fold_difference_.begin(), fold_difference_.end(), no_phase);
        voiced_ = false;
    }

    void fold_pulse(double width) {
        const int size = static_cast<int>(folded_.size());
        for (int j = 0; j < size; ++j) {
            const double u = double(j)/size, back = 1.-u, blend = smooth(u);
            const double v = (1.-blend)*smooth(back/.16), w = blend*smooth(u/.16);
            const double a = periodic(cycle_, u*width), b = periodic(cycle_, 1.-back*width);
            folded_[j] = Complex(static_cast<float>((v*a+w*b)/std::max(.32, v+w)), 0.f);
        }
        fft_.run(folded_);
    }

    double continuous_phase(int k, double measured, double phase_step) {
        auto& previous = output_phase_[k];
        if (std::isfinite(previous)) {
            // Keep adjacent table coefficients far from antiphase. The
            // frequency-relative limit also bounds extra pitch modulation.
            const double limit = std::min(pi/9., k*phase_step);
            measured = previous + std::clamp(std::remainder(measured-previous, 2*pi), -limit, limit);
        }
        previous = std::remainder(measured, 2*pi);
        return previous;
    }

public:
    explicit Tables(int size) : fft_(size), spec_(size), folded_(size), output_spec_(size+1),
        cycle_(size), bell_(size), source_power_(size/2), source_angle_(size/2), log_amp_(size/2),
        linear_amp_(size/2), phase_(size/2), arcs_(size/2, no_phase),
        output_phase_(size/2, no_phase), fold_difference_(size/2, no_phase) {
        for (int j = 0; j < size; ++j)
            bell_[j] = static_cast<float>(.5+.5*std::cos(pi*j/size));
    }

    void build(Table& out, SourceCycles& source, double at, int sr, double hz, float weight, int elapsed,
               const PhaseReference* phase_reference, SourceCycles* original_cycles, HarmonicFrame* frame = nullptr) {
        out.weight = weight;
        if (frame) frame->hz = static_cast<float>(hz);
        if (weight <= 0.f) { reset_phase(); return; }
        const int size = static_cast<int>(spec_.size()), half = size/2;
        // A held sustain can read exactly the same source point for many
        // output frames. Reuse stationary measurements, not target phase.
        const bool new_source = at != spectrum_at_;
        if (new_source) {
            spectrum_at_ = at;
            spectrum_hz_ = source.read(at, cycle_);
            for (int j = 0; j < size; ++j) spec_[j] = Complex(cycle_[j], 0.f);
            fft_.run(spec_);
            unshaped_spec_.clear();
            folded_hz_ = no_phase;
        }
        const double source_hz = spectrum_hz_;
        const int old_last = std::clamp(static_cast<int>(.485*sr/source_hz), 1, half-1);
        const int new_last = std::clamp(static_cast<int>(.485*sr/hz), 1, half-1);
        // These measurements depend only on the source FFT. Keep the raw
        // angles here; the moving origin and temporal arcs are updated below.
        if (new_source) {
            old_power_ = 0.;
            peak_ = 0.f;
            for (int k = 1; k <= old_last; ++k) {
                const float power = std::norm(spec_[k]);
                source_power_[k] = power;
                old_power_ += power;
                peak_ = std::max(peak_, power);
            }
            if (old_power_ >= 1.e-16) {
                const float floor = std::max(1.e-20f, peak_*1.e-12f);
                for (int k = 1; k <= old_last; ++k) {
                    const float bounded = std::max(floor, source_power_[k]);
                    log_amp_[k] = .5f*std::log(bounded);
                    linear_amp_[k] = std::sqrt(bounded);
                    source_angle_[k] = std::arg(spec_[k]);
                }
            }
        }
        const double old_power = old_power_;
        const float peak = peak_;
        if (old_power < 1.e-16) { out.weight = 0.f; reset_phase(); return; }

        const double fundamental_power = source_power_[1];
        const float fundamental_amp = std::abs(spec_[1]);
        // A weak H1 is an unreliable source phase reference. Hold the last
        // reference through a weak patch and keep the output origin fixed
        // to the target clock throughout the voiced region.
        const double reference = smooth((std::sqrt(fundamental_power/peak)-.02)/.08);
        const double origin = source_angle_[1];
        if (!voiced_) {
            source_origin_ = output_origin_ = reference > 0. ? origin : 0.;
            voiced_ = true;
        } else {
            source_origin_ += reference*std::remainder(origin-source_origin_, 2*pi);
        }
        // Ease out of the sine case. In a deep downshift this also keeps
        // the low partials sampled below source H1 from appearing all at once.
        const float overtone_mix = static_cast<float>(smooth(((old_power-fundamental_power)/old_power-2.e-4)/.006));
        const float floor = std::max(1.e-20f, peak*1.e-12f);
        const float amplitude_floor = std::sqrt(floor);
        const float phase_floor = peak*1.e-6f;
        const double source_origin = std::remainder(source_origin_, 2*pi);
        for (int k = 1; k <= old_last; ++k) {
            phase_[k] = static_cast<float>(source_angle_[k]-k*source_origin);
            if (k == 1) continue;
            auto& arc = arcs_[k-1];
            const double difference = double(phase_[k])-phase_[k-1];
            // The shortest arc is ambiguous around +/-pi. Unwrap the arc
            // over time so a branch crossing does not flip the interpolant.
            if (source_power_[k] > phase_floor && source_power_[k-1] > phase_floor) {
                arc = std::isfinite(arc) ? arc + std::remainder(difference-arc, 2*pi)
                                         : std::remainder(difference, 2*pi);
            } else {
                arc = no_phase;
            }
        }
        std::fill(arcs_.begin()+old_last, arcs_.end(), no_phase);

        // Use overlap-added pulses as a phase reference while retaining the
        // harmonic magnitudes and the target pitch clock.
        const double pitch_ratio = hz/source_hz;
        double pulse_weight = overtone_mix*smooth((pitch_ratio-1.025)/.23);
        double pulse_origin = 0., pulse_peak = 0.;
        if (pulse_weight > 0.) {
            // Only the unrotated measurement is reusable. The phase
            // difference and bandwidth limits below still advance each frame.
            if (hz != folded_hz_) {
                if (source.pulse(sr/hz, bell_, folded_))
                    fft_.run(folded_);
                else
                    fold_pulse(1./pitch_ratio);
                folded_hz_ = hz;
            }
            float peak_power = 0.f;
            for (int k = 1; k <= new_last; ++k) peak_power = std::max(peak_power, std::norm(folded_[k]));
            pulse_peak = std::sqrt(peak_power);
            const double reference = std::abs(folded_[1])/std::max(1.e-10, pulse_peak);
            pulse_weight *= smooth((reference-.025)/.10);
            pulse_origin = std::arg(folded_[1]);
        }
        if (pulse_weight <= 0.) std::fill(fold_difference_.begin(), fold_difference_.end(), no_phase);

        // Reuse the phase model's tract envelope and glottal spectrum.
        // Keep the baseline H1 phase and the continuous target cycle clock.
        model_phase_.clear();
        model_amplitude_.clear();
        const auto* phase_spectrum = &spec_;
        double phase_blend = 0., amplitude_blend = 0., model_source_peak = 0.;
        int model_source_orders = 0;
        if (phase_reference) {
            const double rd = phase_reference->glottal.sample(at);
            const double pitch_change = std::abs(std::log2(pitch_ratio));
            const double support = reference*overtone_mix*smooth((weight-.35)/.55);
            phase_blend = support*smooth((pitch_change+
                                         std::abs(std::log2(phase_reference->formant_ratio)))/.065);
            // Flat low-order spectra give little evidence for replacing the
            // measured envelope with a speech source-filter model.
            const int examined = std::min(old_last, static_cast<int>(2400./source_hz));
            if (examined >= 4) {
                const auto range = std::minmax_element(log_amp_.begin()+1, log_amp_.begin()+examined+1);
                const double contrast = *range.second-*range.first;
                amplitude_blend = support*smooth(pitch_change/.065)*smooth((contrast-.34)/.90);
            }
            if (phase_blend > 0. && rd > 0.) {
                if (original_cycles) {
                    if (unshaped_spec_.empty()) {
                        original_cycle_.resize(size);
                        original_cycles->read(at, original_cycle_);
                        unshaped_spec_.resize(size);
                        for (int j = 0; j < size; ++j) unshaped_spec_[j] = Complex(original_cycle_[j], 0.f);
                        fft_.run(unshaped_spec_);
                    }
                    phase_spectrum = &unshaped_spec_;
                }
                const int orders = std::min({80, old_last, static_cast<int>(12000./source_hz)});
                model_source_peak = source_filter_spectrum(*phase_spectrum, orders, source_hz, hz, rd,
                    phase_reference->formant_ratio, model_phase_,
                    amplitude_blend > 0. ? &model_amplitude_ : nullptr, filter_work_);
                model_source_orders = orders;
            }
        }

        std::fill(output_spec_.begin(), output_spec_.end(), Complex());
        double new_power = 0.;
        const double phase_step = 2*pi*.004*hz*elapsed/sr;
        const int last = overtone_mix > 0.f ? new_last : 1;
        for (int k = 1; k <= last; ++k) {
            float amp = fundamental_amp;
            double phase = output_origin_;
            if (overtone_mix > 0.f) {
                const double source_order = k*hz/source_hz;
                const double at = std::clamp(source_order, 1., double(old_last));
                const int a = static_cast<int>(at), b = std::min(a+1, old_last);
                const float u = static_cast<float>(at-a);
                const float log_amp = log_amp_[a] + u*(log_amp_[b]-log_amp_[a]);
                const float wa = (1.f-u)*linear_amp_[a];
                const float wb = u*linear_amp_[b];
                const double arc = std::isfinite(arcs_[a]) ? arcs_[a] :
                                   std::remainder(double(phase_[b])-phase_[a], 2*pi);
                phase = phase_[a] + wb/(wa+wb)*arc + k*output_origin_;
                const float edge = static_cast<float>(1.-smooth(source_order-old_last));
                amp = edge*std::max(0.f, std::exp(log_amp)-amplitude_floor)*overtone_mix;
                if (k == 1) {
                    amp += fundamental_amp*(1.f-overtone_mix);
                    phase = output_origin_ + overtone_mix*std::remainder(phase-output_origin_, 2*pi);
                }
            }
            if (k < static_cast<int>(model_amplitude_.size()) && amp > amplitude_floor) {
                const double band = 1.-smooth((k*hz-6000.)/2000.);
                const float model = std::clamp(model_amplitude_[k], .25f*amp, 4.f*amp);
                amp += static_cast<float>(amplitude_blend*band)*(model-amp);
            }
            if (amp <= amplitude_floor) continue;
            if (pulse_weight > 0.) {
                const double strength = std::abs(folded_[k])/std::max(1.e-10, pulse_peak);
                const double gain = pulse_weight*smooth((strength-.015)/.07);
                auto& difference = fold_difference_[k];
                if (gain > 0.) {
                    const double measured = std::arg(folded_[k]) + k*(output_origin_-pulse_origin);
                    const double change = measured-phase;
                    difference = std::isfinite(difference) ? difference + std::remainder(change-difference, 2*pi)
                                                           : std::remainder(change, 2*pi);
                    phase += gain*difference;
                } else {
                    difference = no_phase;
                }
            }
            if (k > 1 && k < static_cast<int>(model_phase_.size()) && std::isfinite(output_phase_[1])) {
                const double source_strength = k <= model_source_orders ?
                    filter_work_.source_magnitude[k]/model_source_peak : 1.;
                const double band = 1.-smooth((k*hz-6000.)/2000.);
                const double confidence = phase_blend*band*smooth((source_strength-.0015)/.012);
                const double desired = k*output_phase_[1]+model_phase_[k];
                phase += confidence*std::remainder(desired-phase, 2.*pi);
            }
            phase = continuous_phase(k, phase, phase_step);
            output_spec_[k] = Complex(amp*static_cast<float>(std::cos(phase)), amp*static_cast<float>(std::sin(phase)));
            new_power += double(amp)*amp;
        }
        const double match = std::sqrt(old_power/std::max(1.e-16, new_power));
        // A sparse target grid can miss a narrow resonance. Limit how far
        // level matching may amplify the remaining weak partials.
        const double limit = 2.*std::sqrt(std::max(1., hz/source_hz));
        const float gain = static_cast<float>(std::min(match, limit)/size);
        if (frame) {
            if (new_power <= 0.) return;
            frame->support = weight;
            frame->harmonics.resize(last);
            const double origin = std::arg(output_spec_[1]);
            for (int k = 1; k <= last; ++k) {
                const auto rotate = std::polar(1.f, static_cast<float>(-k*origin));
                frame->harmonics[k-1] = output_spec_[k] * rotate * (weight*gain);
            }
            frame->harmonics[0] = {std::abs(output_spec_[1])*weight*gain, 0.f};
            return;
        }
        // Zero padding changes the sample grid, not the input coefficient scale.
        if (!output_fft_) output_fft_.emplace(2*size);
        output_fft_->run(output_spec_, out.wave, gain);
    }
};
int table_size(const std::vector<double>& marks, const std::vector<PitchEpoch>& epochs, int sr) {
    double longest = 1.;
    for (const auto& p : epochs) longest = std::max(longest, sr/p.hz);
    for (std::size_t j = 1; j < marks.size(); ++j)
        longest = std::max(longest, marks[j]-marks[j-1]);
    int size = 128;
    while (size < std::min(8192., 1.20*longest)) size *= 2;
    return size;
}

class TableTrack {
    const std::vector<PitchEpoch>& epochs_;
    const HarmonicSplit& split_;
    int rate_;
    Tables tables_;
    SourceCycles cycles_;
    std::optional<SourceCycles> original_;
    const PhaseReference* reference_;

public:
    std::size_t index = 0;

    TableTrack(const std::vector<float>& source, const std::vector<double>& marks,
               const std::vector<PitchEpoch>& epochs, const HarmonicSplit& split, int sr,
               const PhaseReference* reference)
        : epochs_(epochs), split_(split), rate_(sr), tables_(table_size(marks, epochs, sr)),
          cycles_(source, marks, split, sr), reference_(reference) {
        if (reference && &reference->harmonic != &source)
            original_.emplace(reference->harmonic, marks, split, sr);
    }

    double phase_at(int n) {
        while (index+2 < epochs_.size() && epochs_[index+1].output <= n) ++index;
        return (n-epochs_[index].output)/(epochs_[index+1].output-epochs_[index].output);
    }

    void make(Table& table, int n, int elapsed, HarmonicFrame* frame = nullptr) {
        const double u = std::clamp(phase_at(n), 0., 1.);
        const auto& a = epochs_[index];
        const auto& b = epochs_[index+1];
        const double at = a.source + u*(b.source-a.source), hz = a.hz + u*(b.hz-a.hz);
        const auto weight = static_cast<float>(smooth(split_.activity(at)/.25));
        tables_.build(table, cycles_, at, rate_, hz, weight, elapsed, reference_,
                      original_ ? &*original_ : nullptr, frame);
    }
};
} // namespace

GlottalTrack estimate_glottal_track(const std::vector<float>& source, const std::vector<double>& marks,
                                   const HarmonicSplit& split, int sr) {
    GlottalTrack result;
    if (marks.size() < 2 || source.empty()) return result;
    result.hop = std::max(1, static_cast<int>(std::lround(.010*sr)));
    const int length = static_cast<int>(source.size());
    const int count = 1+(length+result.hop-1)/result.hop;
    std::vector<double> rd(count);
    double longest = 1.;
    for (std::size_t j = 1; j < marks.size(); ++j)
        longest = std::max(longest, marks[j]-marks[j-1]);
    int size = 128;
    while (size < std::min(2048., 1.20*longest)) size *= 2;
    SourceCycles cycles(source, marks, split, sr);
    const Fft fft(size);
    std::vector<float> cycle(size);
    std::vector<Complex> spec(size);
    for (int i = 0; i < count; ++i) {
        const double at = std::min(i*result.hop, length-1);
        if (split.activity(at) < .08) continue;
        const double f0 = cycles.read(at, cycle);
        for (int j = 0; j < size; ++j) spec[j] = Complex(cycle[j], 0.f);
        fft.run(spec);
        rd[i] = estimate_rd(spec, std::min(size/2-1, static_cast<int>(.485*sr/f0)), f0);
    }
    result.rd = rd;
    for (int i = 0; i < count; ++i) {
        if (rd[i] <= 0.) continue;
        // Symmetric triangular smoothing, only within the same voiced run.
        double sum = 3.*rd[i], weight = 3.;
        for (int sign : {-1, 1}) {
            for (int k = 1; k <= 2; ++k) {
                const int j = i+sign*k;
                if (j < 0 || j >= count || rd[j] <= 0.) break;
                sum += (3.-k)*rd[j];
                weight += 3.-k;
            }
        }
        result.rd[i] = sum/weight;
    }
    return result;
}

void synthesize_harmonics(const std::vector<float>& source, const std::vector<double>& marks,
                          const std::vector<PitchEpoch>& epochs, const HarmonicSplit& split,
                           int sr, std::vector<float>& wave, const TargetClock& clock,
                           const PhaseReference* phase_reference) {
    if (wave.empty() || marks.size() < 2 || epochs.size() < 2) return;
    TableTrack track(source, marks, epochs, split, sr, phase_reference);
    static const TableReader reader;
    const int len = static_cast<int>(wave.size());
    const int hop = std::max(1, static_cast<int>(std::lround(.006*sr)));
    Table a, b;
    track.make(a, 0, 0);
    track.index = 0;
    // Integrate the prescribed pitch, not the average between two pulses.
    // The midpoint step gives the exact integral of each linear interval.
    int knot = 0, next_knot = 0;
    double phase = clock.phase, advance = 0., slope = 0.;
    for (int start = 0; start < len; start += hop) {
        const int end = std::min(len, start+hop);
        const auto saved = track.index;
        track.make(b, end, end-start);
        track.index = saved;
        const bool audible = a.weight > 0.f || b.weight > 0.f;
        if (audible || clock.hop > 0) {
            for (int n = start; n < end; ++n) {
                if (audible) {
                    const double position = clock.hop > 0 ? phase : std::clamp(track.phase_at(n), 0., 1.);
                    const float w = static_cast<float>(smooth(double(n-start)/(end-start)));
                    const float l = (1.f-w)*a.weight, r = w*b.weight;
                    wave[n] = reader.read(a.wave, b.wave, position, l, r);
                }
                if (clock.hop > 0) {
                    if (n == next_knot) {
                        next_knot = std::min(n+clock.hop, len);
                        const double f = clock.f0[knot++];
                        slope = (clock.f0[knot]-f)/(sr*double(next_knot-n));
                        advance = f/sr+.5*slope;
                    }
                    phase += advance;
                    if (phase >= 1.) phase -= 1.;
                    advance += slope;
                }
            }
        }
        std::swap(a,b);
    }
}
void prepare_harmonics(const std::vector<float>& source, const std::vector<double>& marks,
                       const std::vector<PitchEpoch>& epochs, const HarmonicSplit& split,
                       int sr, int active, const PhaseReference* reference, NoteModel& note) {
    TableTrack track(source, marks, epochs, split, sr, reference);
    note.hop = std::max(1, static_cast<int>(std::lround(.006*sr)));
    const int count = 1 + (active+note.hop-1)/note.hop;
    note.frames.resize(count);
    Table unused;
    for (int j = 0, previous = 0; j < count; ++j) {
        const int at = std::min(active, j*note.hop);
        track.make(unused, at, at-previous, &note.frames[j]);
        previous = at;
    }
}

} // namespace tn::detail
