// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "pitch.hpp"

#include <vector>

namespace tn::detail {

// Attenuate breath-like parts of the five-cycle estimate. The input supplies
// spectral evidence, and only the estimate supplies audio to remove.
void return_breath(const std::vector<float>& input, int sr, const std::vector<double>& marks,
                   const std::vector<pitch::Island>& voiced, std::vector<float>& average, double amount);
} // namespace tn::detail
