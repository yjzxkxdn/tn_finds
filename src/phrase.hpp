// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "note_model.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace tn::wavtool {
inline constexpr int rate = 44100;
inline constexpr double frames_per_ms = rate / 1000.;

struct Clip {
    int left = 0, end = 0;
    double offset = 0;
    std::array<double, 7> positions{}, gains{};
    NoteModel note;

    double gain(int at) const;
};

// Mix the recording and synthesize the combined harmonic model on one clock.
std::vector<float> synthesize_phrase(const std::vector<Clip>& clips, int size, int sample_rate = rate);

// The classic cache path holds a versioned model in deferred mode. PCM or a
// missing path is reported as no model. A damaged native model is an error.
bool read_model(const std::filesystem::path& path, NoteModel& model);
void write_model(const std::filesystem::path& path, const NoteModel& model);

struct Journal {
    int count = 0;
    int expected = 0; // zero for explicit LAST_NOTE protocol
    double end = 0;
    std::uint64_t bytes = 0;
};
Journal read_journal(const std::filesystem::path& path);
std::vector<Clip> read_clips(const std::filesystem::path& path, const Journal& journal);
void save_clip(const std::filesystem::path& path, const Journal& journal,
               int expected, double end, const Clip& clip);
} // namespace tn::wavtool
