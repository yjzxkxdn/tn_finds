// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "pitch.hpp"

namespace tn { struct Options; struct HarmonicFrame; }

namespace tn::detail {
// Shift a smooth spectral envelope with positive real STFT gains.
// The minimum half-width allows more averaging for aperiodic audio.
void shift_formant_envelope(std::vector<float>& x, int sr, double ratio, const pitch::Track& track,
                            double min_smoothing_hz = 0.);

// Change only spectral magnitude. Mt is harmonic-only; the caller limits
// both stems' corrections to voiced support.
struct GlottalTrack;
void shape_voice_envelope(std::vector<float>& x, int sr, const pitch::Track& track,
                          const Options& options, bool harmonic, const GlottalTrack* glottal = nullptr);

void color_harmonics(HarmonicFrame& frame, int sample_rate, const Options& options,
                     double rd, double voice_mix, double level_hz);
} // namespace tn::detail
