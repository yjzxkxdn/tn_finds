// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <filesystem>

namespace tn::wavtool {
// Recognize the stock Windows temp.bat service contract. Return zero for an
// unknown batch shape. The script is read as data and is never invoked here.
int classic_calls(const std::filesystem::path& output);
// Identify a PCM cache hit for a note written by this dual-role build.
bool is_dual_tool_cache(const std::filesystem::path& note);
} // namespace tn::wavtool
