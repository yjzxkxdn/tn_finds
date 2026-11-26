// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "pitch.hpp"

#include <cstddef>
#include <filesystem>
#include <string>

namespace tn::pitch {
// Read the exact UTF-16 basename in an MRQ v1/v2 descriptor, on its native
// zero-origin source grid. A missing, busy or invalid entry returns {};
// an all-unvoiced entry is a successful nonempty track. Never writes the MRQ.
Track read_mrq(const std::filesystem::path& mrq, const std::u16string& basename,
               int rate, std::size_t samples);

// Producer side of the same format, byte-compatible with Moresampler's
// mrq_write_entry: a matching entry is rewritten in place when its old payload
// has room for the new trajectory, otherwise the old record is tombstoned and a
// fresh one is appended. The descriptor is locked on its version field for the
// duration so concurrent readers never see a half-written record. Returns false
// when the file cannot be created, locked or written; the existing content is
// left untouched in that case as far as the OS allows.
bool write_mrq(const std::filesystem::path& mrq, const std::u16string& basename,
               const Track& track, int rate);
} // namespace tn::pitch
