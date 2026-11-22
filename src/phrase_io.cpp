// SPDX-License-Identifier: GPL-3.0-only
#include "phrase.hpp"
#include "wav.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace tn::wavtool {
namespace {
namespace fs = std::filesystem;
constexpr char note_tag[8] = {'T','N','F','M','O','D','1',0};
constexpr char journal_tag[8] = {'T','N','F','P','H','R','1',0};
constexpr std::uint64_t max_bytes = 384ull*1024*1024;
constexpr std::uint64_t journal_header = 32;
constexpr int max_calls = 20000;
static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559);
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);

[[noreturn]] void bad(const fs::path& path, const char* message) {
    throw std::runtime_error("phrase '"+path.u8string()+"': "+message);
}

bool little_endian() {
    const std::uint16_t one = 1;
    return *reinterpret_cast<const unsigned char*>(&one) == 1;
}

class Reader {
    std::istream& stream_;
    const fs::path& path_;
    std::uint64_t remaining_;
public:
    Reader(std::istream& stream, const fs::path& path, std::uint64_t size)
        : stream_(stream), path_(path), remaining_(size) {}

    std::uint64_t remaining() const { return remaining_; }
    void check(bool condition, const char* message) const { if (!condition) bad(path_,message); }
    void bytes(void* data, std::size_t size) {
        check(size <= remaining_, "truncated parameter file");
        if (size && !stream_.read(static_cast<char*>(data), static_cast<std::streamsize>(size)))
            bad(path_, "parameter read failed");
        remaining_ -= size;
    }
    template<class T> T value() {
        T x;
        unsigned char data[sizeof(T)];
        bytes(data, sizeof(data));
        if (!little_endian()) std::reverse(std::begin(data), std::end(data));
        std::memcpy(&x, data, sizeof(x));
        return x;
    }
    void floats(std::vector<float>& v, std::size_t n) {
        check(n <= remaining_/sizeof(float), "truncated sample data");
        v.resize(n);
        if (little_endian()) bytes(v.data(), n*sizeof(float));
        else for (auto& x : v) x = value<float>();
    }
};

class Writer {
    std::ostream& stream_;
    const fs::path& path_;
public:
    Writer(std::ostream& stream, const fs::path& path) : stream_(stream), path_(path) {}
    void bytes(const void* data, std::size_t size) {
        if (size && !stream_.write(static_cast<const char*>(data), static_cast<std::streamsize>(size)))
            bad(path_, "parameter write failed");
    }
    template<class T> void value(T x) {
        unsigned char data[sizeof(T)];
        std::memcpy(data, &x, sizeof(x));
        if (!little_endian()) std::reverse(std::begin(data), std::end(data));
        bytes(data, sizeof(data));
    }
    void floats(const std::vector<float>& v) {
        if (little_endian()) bytes(v.data(), v.size()*sizeof(float));
        else for (const auto x : v) value(x);
    }
};

bool in_range(double x, double low, double high) { return x >= low && x <= high; }
int frame_count(int n, int hop) { return 1+(n+hop-1)/hop; }

NoteModel model_body(Reader& in) {
    NoteModel m;
    m.sample_rate = in.value<std::int32_t>();
    m.active = in.value<std::int32_t>();
    m.hop = in.value<std::int32_t>();
    m.pitch_hop = in.value<std::int32_t>();
    const auto samples = in.value<std::uint32_t>();
    const auto frames = in.value<std::uint32_t>();
    in.check(in_range(m.sample_rate, 8000, 192000) && samples <= wav::max_samples &&
             in_range(m.active, 0, samples), "invalid parameter dimensions");
    if (frames) {
        in.check(m.active > 0 && m.hop == std::lround(.006*m.sample_rate) &&
                 m.pitch_hop == std::lround(.001*m.sample_rate) &&
                 frames == static_cast<unsigned>(frame_count(m.active, m.hop)), "invalid frame grid");
        in.floats(m.f0, frame_count(m.active, m.pitch_hop));
        for (const auto hz : m.f0) in.check(in_range(hz, 30., 2400.), "invalid target F0");
    }
    in.floats(m.residual, samples);
    for (const auto x : m.residual) in.check(std::isfinite(x), "invalid residual sample");
    in.check(frames <= in.remaining()/12, "truncated harmonic frames");
    m.frames.resize(frames);
    for (auto& frame : m.frames) {
        frame.hz = in.value<float>();
        frame.support = in.value<float>();
        const auto orders = in.value<std::uint32_t>();
        in.check(in_range(frame.hz, 29., 2401.) && in_range(frame.support, 0., 1.) &&
                 orders <= 4095 && orders <= in.remaining()/8, "invalid harmonic frame");
        frame.harmonics.resize(orders);
        for (auto& coefficient : frame.harmonics) {
            const float real = in.value<float>(), imag = in.value<float>();
            in.check(in_range(real, -1.e6, 1.e6) && in_range(imag, -1.e6, 1.e6), "invalid harmonic coefficient");
            coefficient = {real, imag};
        }
        in.check(!orders || (frame.harmonics[0].real() >= 0.f && frame.harmonics[0].imag() == 0.f),
                 "invalid H1 phase reference");
    }
    return m;
}

void model_body(Writer& out, const NoteModel& m) {
    out.value<std::int32_t>(m.sample_rate);
    out.value<std::int32_t>(m.active);
    out.value<std::int32_t>(m.hop);
    out.value<std::int32_t>(m.pitch_hop);
    out.value<std::uint32_t>(static_cast<std::uint32_t>(m.residual.size()));
    out.value<std::uint32_t>(static_cast<std::uint32_t>(m.frames.size()));
    out.floats(m.f0);
    out.floats(m.residual);
    for (const auto& frame : m.frames) {
        out.value(frame.hz);
        out.value(frame.support);
        out.value<std::uint32_t>(static_cast<std::uint32_t>(frame.harmonics.size()));
        for (const auto coefficient : frame.harmonics) {
            out.value(coefficient.real());
            out.value(coefficient.imag());
        }
    }
}

std::uint64_t model_bytes(const NoteModel& m) {
    std::uint64_t bytes = 24u + 4ull*(m.residual.size()+m.f0.size());
    for (const auto& frame : m.frames) bytes += 12u + 8ull*frame.harmonics.size();
    return bytes;
}

void close(std::ofstream& stream, const fs::path& path) {
    stream.close();
    if (!stream) bad(path, "could not finish writing parameters");
}
} // namespace

bool read_model(const std::filesystem::path& path, NoteModel& model) {
    std::ifstream stream(path, std::ios::binary | std::ios::ate);
    const auto size = static_cast<std::streamoff>(stream.tellg());
    if (!stream || size < 8) return false;
    stream.seekg(0);
    char tag[sizeof(note_tag)];
    if (!stream.read(tag, sizeof(tag)) || std::memcmp(tag, note_tag, 6) != 0) return false;
    if (std::memcmp(tag, note_tag, sizeof(tag)) != 0 || size > static_cast<std::streamoff>(max_bytes))
        bad(path, "unsupported note model version or size");
    Reader in(stream, path, size-sizeof(tag));
    model = model_body(in);
    in.check(in.remaining() == 0, "extra data after note model");
    return true;
}

void write_model(const std::filesystem::path& path, const NoteModel& model) {
    if (model_bytes(model)+sizeof(note_tag) > max_bytes) bad(path, "note model exceeds size limit");
    // Write a complete model beside the cache before making it visible to the host.
    auto directory = path;
    directory += ".tnf-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    if (!fs::create_directory(directory)) bad(path, "could not reserve model staging directory");
    const auto temp = directory/"note";
    try {
        std::ofstream stream(temp, std::ios::binary);
        Writer out(stream, temp);
        out.bytes(note_tag, sizeof(note_tag));
        model_body(out, model);
        close(stream, temp);
        fs::rename(temp, path);
    } catch (...) {
        std::error_code ignored;
        fs::remove(temp, ignored);
        fs::remove(directory, ignored);
        throw;
    }
    std::error_code ignored;
    fs::remove(directory, ignored);
}

Journal read_journal(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary | std::ios::ate);
    const auto size = static_cast<std::streamoff>(stream.tellg());
    if (!stream || size < static_cast<std::streamoff>(journal_header)) bad(path, "missing phrase journal");
    stream.seekg(0);
    Reader in(stream, path, size);
    char tag[sizeof(journal_tag)];
    in.bytes(tag, sizeof(tag));
    in.check(std::memcmp(tag, journal_tag, sizeof(tag)) == 0, "unsupported phrase journal version");
    Journal j;
    j.expected = in.value<std::int32_t>();
    j.count = in.value<std::int32_t>();
    j.end = in.value<double>();
    j.bytes = in.value<std::uint64_t>();
    in.check(in_range(j.expected, 0, max_calls) && in_range(j.count, 1, max_calls) &&
             (!j.expected || j.count < j.expected) && in_range(j.end, 0, wav::max_samples) &&
             j.bytes > journal_header && j.bytes <= max_bytes && j.bytes <= static_cast<std::uint64_t>(size),
             "invalid phrase journal header");
    return j;
}

std::vector<Clip> read_clips(const std::filesystem::path& path, const Journal& journal) {
    std::vector<Clip> clips;
    if (journal.count == 0) return clips;
    std::ifstream stream(path, std::ios::binary);
    stream.seekg(journal_header);
    Reader in(stream, path, journal.bytes-journal_header);
    clips.reserve(journal.count);
    for (int i = 0; i < journal.count; ++i) {
        Clip clip;
        clip.left = in.value<std::int32_t>();
        clip.end = in.value<std::int32_t>();
        clip.offset = in.value<double>();
        for (auto& t : clip.positions) t = in.value<double>();
        for (auto& v : clip.gains) v = in.value<double>();
        in.check(in_range(clip.left, -441.e6, wav::max_samples) && in_range(clip.end, 0, wav::max_samples) &&
                 clip.left <= clip.end && in_range(clip.offset, -1.e7, 1.e7) &&
                 std::is_sorted(clip.positions.begin(), clip.positions.end()), "invalid clip timing");
        for (auto t : clip.positions) in.check(in_range(t, 0, wav::max_samples/frames_per_ms), "invalid envelope time");
        for (auto v : clip.gains) in.check(in_range(v, 0, 10.), "invalid envelope gain");
        clip.note = model_body(in);
        clips.push_back(std::move(clip));
    }
    in.check(in.remaining() == 0, "extra data after phrase clips");
    return clips;
}

void save_clip(const std::filesystem::path& path, const Journal& journal,
               int expected, double end, const Clip& clip) {
    const auto pos = journal.count ? journal.bytes : journal_header;
    const auto bytes = pos+128u+model_bytes(clip.note);
    if (bytes > max_bytes || journal.count+1 >= max_calls) bad(path, "phrase journal exceeds size limit");
    auto mode = std::ios::binary | std::ios::out | (journal.count ? std::ios::in : std::ios::trunc);
    std::ofstream stream(path, mode);
    stream.seekp(static_cast<std::streamoff>(pos));
    Writer out(stream, path);
    out.value<std::int32_t>(clip.left);
    out.value<std::int32_t>(clip.end);
    out.value(clip.offset);
    for (auto t : clip.positions) out.value(t);
    for (auto v : clip.gains) out.value(v);
    model_body(out, clip.note);
    // Advance the committed cursor after writing the complete entry.
    stream.seekp(0);
    out.bytes(journal_tag, sizeof(journal_tag));
    out.value<std::int32_t>(expected);
    out.value<std::int32_t>(journal.count+1);
    out.value(end);
    out.value(bytes);
    close(stream, path);
}
} // namespace tn::wavtool
