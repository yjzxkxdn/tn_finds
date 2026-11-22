// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "dsp_fft.hpp"

namespace tn::detail {
// Relative phase and, when requested, amplitudes on the target harmonic grid.
// Analyze the coherent cycle before formant EQ so its source spectrum is unedited.
void source_filter_spectrum(const std::vector<Complex>& spectrum, int source_orders,
                            double source_hz, double target_hz, double rd, double formant_ratio,
                            std::vector<double>& phase, std::vector<float>* amplitude);
// Both mark trains have the same phase origin in source time. Target-center
// positions can be stretched, with one source coordinate per target mark.
// Transform the coherent waveform only; the caller mixes the same aperiodic
// path and wet/DC gate as the reference mode.
void lock_harmonic_phase(const std::vector<float>& source, const std::vector<double>& source_marks,
                         const std::vector<double>& target_marks, const std::vector<double>& source_at_target,
                         int sr, int lo, std::vector<float>& wet);
} // namespace tn::detail
