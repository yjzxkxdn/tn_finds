// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "pitch.hpp"

#include <cstddef>
#include <filesystem>

namespace tn::pitch {

// Resolve a stock OpenUtau source alias only if its exact path fingerprint and
// decoded audio match an explicitly registered voicebank recording.
Track lookup_bound_mrq(const std::filesystem::path& input, const std::vector<float>& samples, int rate);

// Register a WAV or WAVs in a chosen voicebank. The index is tool-owned;
// originals and OpenUtau's cache are never modified.
std::size_t bind_mrq_sources(const std::filesystem::path& original);

} // namespace tn::pitch
