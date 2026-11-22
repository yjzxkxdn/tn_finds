// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

namespace tn::wavtool {
// The fourth positional argument is a duration. The preceding argument is
// the numeric input offset, where the resampler expects a note name.
bool is_classic_call(const std::vector<std::filesystem::path>& args);
enum class JoinMode { classic_phrase, explicit_phrase, pcm };
void append(const std::vector<std::filesystem::path>& args, std::size_t first, JoinMode mode);
} // namespace tn::wavtool
