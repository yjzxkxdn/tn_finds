// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace tn::detail {
// Fractional read with a windowed low-pass sinc. A ratio above one
// decimates. The taps outside the chosen source interval read zero.
class Sinc {
    static constexpr double pi = 3.14159265358979323846;
    static constexpr int taps = 32, half = taps / 2, phases = 512;
    std::vector<std::array<float, taps>> table_;
    std::vector<float> padded_;
public:
    Sinc(const std::vector<float>& x, double r)
        : table_(phases), padded_(x.size() + 2 * half + 2) {
        std::copy(x.begin(), x.end(), padded_.begin() + half);
        const double cutoff = std::min(1., 1. / r);
        for (int f = 0; f < phases; ++f) {
            double sum = 0;
            for (int j = 0; j < taps; ++j) {
                const double d = f / static_cast<double>(phases) - (j - half + 1);
                const double sinc = d == 0 ? cutoff : std::sin(pi * cutoff * d) / (pi * d);
                const double w = .42 + .5 * std::cos(pi * d / half) + .08 * std::cos(2 * pi * d / half);
                table_[f][j] = static_cast<float>(sinc * w);
                sum += table_[f][j];
            }
            for (float& coeff : table_[f]) coeff = static_cast<float>(coeff / sum);
        }
        if (cutoff == 1) {
            table_[0].fill(0);
            table_[0][half - 1] = 1;
        }
    }

    float read(double p, int begin, int end) const {
        if (p <= begin - half || p >= end + half) return 0;
        const auto z = static_cast<std::int64_t>(std::floor(p * phases + .50000001));
        const auto phase = static_cast<int>((z % phases + phases) % phases);
        const int first = static_cast<int>((z - phase) / phases) - half + 1;
        const float* row = table_[phase].data();
        if (first < begin || first + taps > end) {
            float y = 0;
            for (int j = std::max(0, begin - first); j < std::min(taps, end - first); ++j)
                y += row[j] * padded_[half + first + j];
            return y;
        }
        const float* s = padded_.data() + half + first;
        float a = 0, b = 0, c = 0, d = 0;
        for (int j = 0; j < taps; j += 4) {
            a += row[j] * s[j]; b += row[j+1] * s[j+1];
            c += row[j+2] * s[j+2]; d += row[j+3] * s[j+3];
        }
        return (a + c) + (b + d);
    }
};
} // namespace tn::detail
