// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "engine.hpp"
#include "pitch.hpp"

#include <vector>

namespace tn::detail {

struct HarmonicSplit {
    std::vector<float> harmonic;
    // A little of the residual may be attenuated in the *strongly periodic*
    // lower bands only. It is warped with the exact same sample weights as H.
    std::vector<float> cleanup;
    std::vector<float> strength, periodic_share;
    int hop = 1;
    double activity(double sample) const;
    double share(double sample) const;
    std::size_t regions() const;
};

// Supply the original automatic F0 track when available. The continuous
// track and marks are used for low-frequency protection and ridge evidence.
HarmonicSplit split_harmonic(const std::vector<float>& x, int sr,
                             const std::vector<double>& marks, const pitch::Track& track,
                             bool suppress_low_residue, double amount = default_mask_return,
                             const pitch::Track* automatic_track = nullptr);

} // namespace tn::detail
