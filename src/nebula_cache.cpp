// SPDX-License-Identifier: GPL-3.0-only
#include "pitch.hpp"
#include "cache_win.hpp"

#include <cmath>
#include <cstdint>

namespace tn::pitch {
namespace {
using namespace tn::detail::cache;
using Analyze = Track (*)(const std::vector<float>&, int);
constexpr std::array<char, 8> magic{'T','N','N','E','B','F','0','1'};
constexpr char estimator_id[] = TN_FINDS_NEBULA_CACHE_ID;
static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559);

struct Header {
    std::array<char, 8> tag = magic;
    Digest audio{};
    Digest payload{};
    std::uint32_t step = 0;
    std::uint32_t frames = 0;
};
static_assert(sizeof(Header) == 80);

Track read_entry(const fs::path& path, const Digest& key, const Sha256& sha, int step, std::size_t samples) {
    std::ifstream in(path, std::ios::binary);
    Header header;
    if (!in.read(reinterpret_cast<char*>(&header), sizeof(header)) || header.tag != magic || header.audio != key ||
        header.step != static_cast<std::uint32_t>(step) || header.frames == 0 ||
        header.frames > 1 + samples/step) return {};
    std::vector<float> values(header.frames);
    const auto size = values.size()*sizeof(float);
    Digest digest;
    if (!in.read(reinterpret_cast<char*>(values.data()), size) || in.peek() != std::ifstream::traits_type::eof() ||
        !sha.hash(digest, {bytes(values)}) || header.payload != digest) return {};
    Track track;
    track.step = step;
    track.frames.reserve(values.size());
    for (const float hz : values) {
        if (hz == 0.f) track.frames.push_back({});
        else if (std::isfinite(hz) && hz >= 30.f && hz <= 2400.f) track.frames.push_back({hz, 1.});
        else return {};
    }
    return track;
}

void write_entry(const fs::path& path, const Digest& key, const Sha256& sha, const Track& track) {
    if (track.frames.empty()) return;
    std::vector<float> values;
    values.reserve(track.frames.size());
    for (const auto& frame : track.frames) values.push_back(static_cast<float>(frame.hz));
    Header header;
    header.audio = key;
    header.step = track.step;
    header.frames = static_cast<std::uint32_t>(values.size());
    if (!sha.hash(header.payload, {bytes(values)})) return;
    publish(path, {bytes(header), bytes(values)});
}
} // namespace

Track cached_nebula(const std::vector<float>& input, int rate, int step, Analyze analyze) {
    using namespace tn::detail::cache;
    const auto dir = directory(L"TN_FINDS_NEBULA_CACHE", L"nebula-f0-v1");
    if (dir.empty()) return analyze(input, rate);
    const Sha256 sha;
    Digest key;
    // Use the decoded audio (not filename, size or mtime) and all estimator
    // inputs. A changed source cannot hit even with restored timestamps.
    if (!sha.hash(key, {bytes(estimator_id), bytes(rate), bytes(input)})) return analyze(input, rate);
    const auto name = hex(key);
    const auto path = dir / (name + L".f0");
    auto hit = read_entry(path, key, sha, step, input.size());
    if (!hit.frames.empty()) return hit;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) return analyze(input, rate);
    const EntryLock lock(L"Local\\tn_finds_nebula_" + name);
    if (!lock.owned) return analyze(input, rate);
    hit = read_entry(path, key, sha, step, input.size());
    if (!hit.frames.empty()) return hit;
    auto track = analyze(input, rate);
    write_entry(path, key, sha, track);
    return track;
}
} // namespace tn::pitch
