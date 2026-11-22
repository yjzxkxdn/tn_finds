// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "harmonic_mask.hpp"

#include <vector>

namespace tn::detail {
struct GlottalTrack;
GlottalTrack estimate_glottal_track(const std::vector<float>& harmonic, const std::vector<double>& marks,
                                   const HarmonicSplit& split, int sample_rate);

struct PitchEpoch {
    double output, source, hz;
};

struct TargetClock {
    int hop = 0;
    double phase = 0.;
    std::vector<double> f0;
};

struct PhaseReference {
    const std::vector<float>& harmonic;
    const GlottalTrack& glottal;
    double formant_ratio;
};

// Resample the measured harmonic envelope onto the target F0 grid, using
// the same continuous target phase for both upward and downward edits.
void synthesize_harmonics(const std::vector<float>& source, const std::vector<double>& marks,
                          const std::vector<PitchEpoch>& epochs, const HarmonicSplit& split,
                          int sample_rate, std::vector<float>& wave, const TargetClock& clock,
                          const PhaseReference* phase_reference);

// Export the same modeled harmonics before they are evaluated on a sample clock.
void prepare_harmonics(const std::vector<float>& source, const std::vector<double>& marks,
                       const std::vector<PitchEpoch>& epochs, const HarmonicSplit& split,
                       int sample_rate, int active, const PhaseReference* phase_reference, NoteModel& note);
} // namespace tn::detail
