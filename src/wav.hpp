// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

namespace tn::wav {

inline constexpr std::size_t max_samples = 16'000'000;

struct Audio {
    std::vector<float> samples;
    int sample_rate;
};

// Read mono or stereo RIFF/WAVE and return finite mono samples.
Audio read(const std::filesystem::path& path);

// Write mono PCM16 or IEEE float32 at the given sample rate.
void write(const std::filesystem::path& path, const std::vector<float>& samples,
           int sample_rate, bool pcm16);

// Follow links and compare file identity where possible. Throw on ambiguity.
void require_distinct(const std::filesystem::path& input,
                      const std::filesystem::path& output);

} // namespace tn::wav
