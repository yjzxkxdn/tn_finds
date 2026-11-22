// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "dsp_fft.hpp"

namespace tn::detail {
struct GlottalTrack {
    int hop = 1;
    std::vector<double> rd;
    double sample(double at) const;
};

// Fit LF shape to the first <= 8 kHz of a phase-aligned source cycle.
double estimate_rd(const std::vector<Complex>& spectrum, int harmonics, double f0);

// Positive-frequency coefficients of the LF flow derivative, indexed by order.
void lf_spectrum(double rd, int harmonics, std::vector<std::complex<double>>& out);

// H1-referenced log-amplitude ratio. Orders zero and one are neutral.
void tension_gain(double rd, double mt, int harmonics, std::vector<double>& log_gain);
} // namespace tn::detail
