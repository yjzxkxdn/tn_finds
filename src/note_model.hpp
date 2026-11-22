// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <complex>
#include <vector>

namespace tn {
// A positive-frequency coefficient contributes 2*Re(H[k]*exp(2*pi*i*k*phase)).
// H1 is real and nonnegative so every note uses the same phase reference.
struct HarmonicFrame {
    float hz = 0;
    float support = 0;
    std::vector<std::complex<float>> harmonics; // H1, H2, ...
};

struct NoteModel {
    int sample_rate = 0;
    int hop = 0;
    int pitch_hop = 0;
    int active = 0;
    std::vector<float> residual;
    std::vector<float> f0;
    std::vector<HarmonicFrame> frames;
};
} // namespace tn
