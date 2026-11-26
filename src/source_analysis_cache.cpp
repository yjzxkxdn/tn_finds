// SPDX-License-Identifier: GPL-3.0-only
#include "source_analysis.hpp"
#include "cache_win.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace tn::detail {
namespace {
using namespace cache;
using Analyze = SourceAnalysis (*)(const SourceAnalysisInput&);
constexpr char analysis_id[] = TN_FINDS_ANALYSIS_CACHE_ID;
constexpr std::array<char, 8> magic{'T','N','O','T','O','A','N','1'};
// Very long OTO regions still render, without adding huge cache entries.
constexpr std::uint64_t max_payload = 16u*1024*1024;
static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559);
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
static_assert(sizeof(pitch::Frame) == 16 && std::is_trivially_copyable_v<pitch::Frame>);
static_assert(sizeof(F0Point) == 16 && std::is_trivially_copyable_v<F0Point>);

// Windows on supported architectures is little-endian. Use fixed-width
// fields and bulk arrays, without serializing any compiler-owned containers.
struct Layout {
    std::uint32_t clock_hop, frames, marks, split_hop, harmonic, cleanup, curve, rd_hop, rd;
    std::uint32_t reserved = 0;
    double first_pulse_region, end_pulse_region;

    std::uint64_t size() const {
        return 16ull*frames + 8ull*(std::uint64_t(marks)+rd) +
               4ull*(std::uint64_t(harmonic)+cleanup+2ull*curve);
    }
};
struct Header {
    std::array<char, 8> tag = magic;
    Digest key{};
    Digest payload{};
    Layout layout{};
};
static_assert(sizeof(Layout) == 56 && sizeof(Header) == 128);

bool analysis_key(const Sha256& sha, Digest& key, const SourceAnalysisInput& in) {
    const auto* pcm = &in.segment;
#ifdef TN_FINDS_NEBULA_F0
    // For automatic Nebula, even samples outside the OTO segment can affect
    // the native pitch lattice. Supplied F0 and the standard backend are local.
    if (in.f0.empty() && in.mrq_f0.frames.empty() && in.forced_hz == 0.) pcm = &in.recording;
#endif
    const auto first = static_cast<std::uint32_t>(in.first);
    const auto& mrq = in.mrq_f0;
    const std::array<std::uint32_t, 9> dims{
        static_cast<std::uint32_t>(in.rate), first, first+static_cast<std::uint32_t>(in.segment.size()),
        static_cast<std::uint32_t>(in.cleanup) | (static_cast<std::uint32_t>(in.need_rd) << 1),
        static_cast<std::uint32_t>(in.f0.size()), static_cast<std::uint32_t>(pcm->size()),
        static_cast<std::uint32_t>(mrq.frames.empty() ? 0 : in.recording.size()),
        static_cast<std::uint32_t>(mrq.frames.empty() ? 0 : mrq.step), static_cast<std::uint32_t>(mrq.frames.size())};
    // The chosen MRQ track (including UV and hop), not its timestamp, governs reuse.
    return sha.hash(key, {bytes(analysis_id), bytes(dims), bytes(in.mask_return), bytes(in.forced_hz),
                          bytes(in.f0), bytes(mrq.frames), bytes(*pcm)});
}

template<class T> bool read(std::ifstream& stream, std::vector<T>& part, std::uint32_t count) {
    part.resize(count);
    return !count || bool(stream.read(reinterpret_cast<char*>(part.data()), count*sizeof(T)));
}

bool inside(double value, double low, double high) { return value >= low && value <= high; }

bool read_entry(const fs::path& path, const Digest& key, const Sha256& sha,
                const SourceAnalysisInput& in, SourceAnalysis& out) {
    std::ifstream stream(path, std::ios::binary | std::ios::ate);
    const auto size = static_cast<std::streamoff>(stream.tellg());
    if (size < static_cast<std::streamoff>(sizeof(Header)) ||
        size > static_cast<std::streamoff>(sizeof(Header)+max_payload)) return false;
    stream.seekg(0);
    Header head;
    const auto& l = head.layout;
    const auto n = in.segment.size();
    if (!stream.read(reinterpret_cast<char*>(&head), sizeof(head)) || head.tag != magic || head.key != key ||
        l.reserved || l.size()+sizeof(head) != static_cast<std::uint64_t>(size) || !l.clock_hop || !l.frames) return false;
    if (in.mrq_f0.frames.empty()) {
        if (l.clock_hop > 1024 || l.frames > 1+(n+l.clock_hop-1)/l.clock_hop) return false;
    } else {
        if (l.clock_hop != static_cast<unsigned>(in.mrq_f0.step)) return false;
        const auto frames = in.first == 0 && n == in.recording.size()
            ? in.mrq_f0.frames.size() : 1+(n-1+l.clock_hop-1)/l.clock_hop;
        if (l.frames != frames) return false;
    }
    if (l.split_hop == 0 || l.split_hop > 2048 || l.curve > (n-1)/l.split_hop+2 ||
        (l.harmonic && (l.harmonic != n || l.marks < 2 || l.curve != (n-1)/l.split_hop+2)) ||
        l.cleanup != (in.cleanup ? l.harmonic : 0) || l.marks == 1 || l.marks > n ||
        (l.marks && (!inside(l.first_pulse_region, 0., double(n)) ||
                     !inside(l.end_pulse_region, l.first_pulse_region, double(n)))) ||
        (l.rd && (!in.need_rd || !l.harmonic || l.rd_hop != static_cast<unsigned>(std::lround(.010*in.rate)) ||
                  l.rd != 1+(n+l.rd_hop-1)/l.rd_hop)) ||
        (!l.rd && (l.rd_hop != 1 || (in.need_rd && l.harmonic)))) return false;

    SourceAnalysis a;
    auto& s = a.split;
    if (!read(stream, a.clock.frames, l.frames) || !read(stream, a.pulses.marks, l.marks) ||
        !read(stream, s.harmonic, l.harmonic) || !read(stream, s.cleanup, l.cleanup) ||
        !read(stream, s.strength, l.curve) || !read(stream, s.periodic_share, l.curve) ||
        !read(stream, a.rd.rd, l.rd)) return false;
    Digest actual;
    if (!sha.hash(actual, {bytes(l), bytes(a.clock.frames), bytes(a.pulses.marks), bytes(s.harmonic),
                           bytes(s.cleanup), bytes(s.strength), bytes(s.periodic_share), bytes(a.rd.rd)}) ||
        actual != head.payload) return false;
    for (const auto& f : a.clock.frames)
        if (!inside(f.hz, 0., 2401.) || !std::isfinite(f.confidence)) return false;
    double previous = -1.;
    for (const auto mark : a.pulses.marks) {
        if (!inside(mark, l.first_pulse_region, l.end_pulse_region) || !(mark > previous)) return false;
        previous = mark;
    }
    const auto finite = [](const auto& v) {
        return std::all_of(v.begin(), v.end(), [](auto x) { return std::isfinite(x); });
    };
    if (!finite(s.harmonic) || !finite(s.cleanup) || !finite(s.strength) || !finite(s.periodic_share) || !finite(a.rd.rd))
        return false;
    a.clock.step = static_cast<int>(l.clock_hop);
    a.pulses.begin = l.first_pulse_region;
    a.pulses.end = l.end_pulse_region;
    a.split.hop = static_cast<int>(l.split_hop);
    a.rd.hop = static_cast<int>(l.rd_hop);
    out = std::move(a);
    return true;
}

void write_entry(const fs::path& path, const Digest& key, const Sha256& sha, const SourceAnalysis& a) {
    const auto u32 = [](auto value) { return static_cast<std::uint32_t>(value); };
    const auto& s = a.split;
    Header head;
    head.key = key;
    head.layout = {u32(a.clock.step), u32(a.clock.frames.size()), u32(a.pulses.marks.size()), u32(s.hop),
                   u32(s.harmonic.size()), u32(s.cleanup.size()), u32(s.strength.size()),
                   u32(a.rd.hop), u32(a.rd.rd.size()), 0, a.pulses.begin, a.pulses.end};
    if (!head.layout.frames || head.layout.size() > max_payload) return;
    if (!sha.hash(head.payload, {bytes(head.layout), bytes(a.clock.frames), bytes(a.pulses.marks), bytes(s.harmonic),
                                 bytes(s.cleanup), bytes(s.strength), bytes(s.periodic_share), bytes(a.rd.rd)})) return;
    publish(path, {bytes(head), bytes(a.clock.frames), bytes(a.pulses.marks), bytes(s.harmonic), bytes(s.cleanup),
                   bytes(s.strength), bytes(s.periodic_share), bytes(a.rd.rd)});
}
} // namespace

SourceAnalysis cached_source_analysis(const SourceAnalysisInput& in, Analyze analyze, bool& hit) {
    if (in.forced_hz < 0.) return {};
    const auto dir = directory(L"TN_FINDS_ANALYSIS_CACHE", L"oto-analysis-v1");
    if (dir.empty()) return analyze(in);
    const Sha256 sha;
    Digest key;
    if (!analysis_key(sha, key, in)) return analyze(in);
    const auto name = hex(key);
    const auto path = dir / (name + L".ana");
    SourceAnalysis a;
    if ((hit = read_entry(path, key, sha, in, a))) return a;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) return analyze(in);
    const EntryLock lock(L"Local\\tn_finds_oto_" + name);
    if (!lock.owned) return analyze(in);
    if ((hit = read_entry(path, key, sha, in, a))) return a;
    a = analyze(in);
    write_entry(path, key, sha, a);
    return a;
}
} // namespace tn::detail
