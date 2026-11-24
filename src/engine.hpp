// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "note_model.hpp"

#include <cstddef>
#include <vector>

namespace tn {
struct F0Point { double time_s, hz; };

inline constexpr double default_mask_return = .75;

struct Options {
    double offset_ms = 0, cutoff_ms = 0, length_ms = -1, consonant_ms = 0, velocity = 100;
    double target_hz = 0, pitch_st = 0, modulation = 0, tempo = 120;
    std::vector<int> bends;  // cents at 1/96 beat, starting at output sample 0
    double formant_st = 0, forced_hz = 0, consonant_boost = 0;
    double opening = 0, dynamics = 0, volume = 100;
    bool loop = false, peak_normalize = false;
    bool grain_formant = false; // opt out of phase-preserving formant EQ
    bool phase_lock = false;    // optional source-relative harmonic phase transfer
    std::vector<F0Point> f0; // time of original (not trimmed) input; 0 Hz supplies no period hint
    double mask_return = default_mask_return; // 0..1, classic Tnr divided by 100
    double residual_mix = 0.; // Tns: post-fixed spectral residual mix, percent
    double fixed_residual_mix = 0.; // Tnc: fixed-region spectral residual mix, percent
    // Moresampler-style controls; all neutral at zero.
    double tension = 0, breathiness = 0, mouth = 0, resonance = 0, formant_emphasis = 0;
    bool force_loop = false; // prefer cycling even for short extensions
};

struct Report {
    double analysis_ms = 0, synthesis_ms = 0, source_pitch_hz = 0;
    std::size_t input_samples = 0, output_samples = 0, islands = 0;
    bool bypass = false, loop_used = false;
    const char* f0_backend = "standard";
};

// Mono finite samples. Output uses input sample rate; no required external libraries.
// Throws std::invalid_argument for unsupported settings or an empty region.
std::vector<float> render(const std::vector<float>& samples, int sample_rate,
                          const Options& options, Report* report = nullptr);

// Prepare a classic note without evaluating its periodic waveform.
NoteModel prepare_note(const std::vector<float>& samples, int sample_rate, const Options& options);
}  // namespace tn
