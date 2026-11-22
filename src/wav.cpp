// SPDX-License-Identifier: GPL-3.0-only
#include "wav.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>

namespace tn::wav {
namespace {
namespace fs = std::filesystem;
using Byte = unsigned char;

static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559,
              "WAV I/O requires IEEE binary32");
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559,
              "WAV I/O requires IEEE binary64");

constexpr int min_rate = 8000;
constexpr int max_rate = 192000;
constexpr std::size_t block_frames = 4096;

[[noreturn]] void bad(const fs::path& path, const std::string& reason) {
    throw std::runtime_error("WAV '" + path.u8string() + "': " + reason);
}

std::uint32_t u16(const Byte* p) {
    return std::uint32_t(p[0]) | (std::uint32_t(p[1]) << 8);
}

std::uint32_t u32(const Byte* p) {
    return u16(p) | (u16(p + 2) << 16);
}

std::uint64_t u64(const Byte* p) {
    return u32(p) | (std::uint64_t(u32(p + 4)) << 32);
}

bool id(const Byte* p, const char (&tag)[5]) {
    return std::memcmp(p, tag, 4) == 0;
}

void read_bytes(std::ifstream& stream, Byte* p, std::size_t count, const fs::path& path) {
    if (!stream.read(reinterpret_cast<char*>(p), static_cast<std::streamsize>(count))) {
        bad(path, "truncated or unreadable file");
    }
}

void seek(std::ifstream& stream, std::uint64_t position, const fs::path& path) {
    if (!stream.seekg(static_cast<std::streamoff>(position), std::ios::beg)) {
        bad(path, "could not seek");
    }
}

struct Format {
    unsigned channels = 0;
    unsigned rate = 0;
    unsigned width = 0;
    unsigned align = 0;
    unsigned valid_bits = 0;
    bool fp = false;
};

Format parse_format(std::ifstream& stream, std::uint32_t size, const fs::path& path) {
    if (size < 16 || size == 17) bad(path, "invalid fmt chunk size");
    std::array<Byte, 40> bytes{};
    read_bytes(stream, bytes.data(), std::min<std::size_t>(size, bytes.size()), path);
    const Byte* b = bytes.data();
    auto code = u16(b);
    const unsigned bits = u16(b + 14);
    const unsigned extra = size >= 18 ? u16(b + 16) : 0;
    if (size >= 18 && extra > size - 18) bad(path, "fmt extension exceeds chunk bounds");
    Format fmt{u16(b + 2), u32(b + 4), bits / 8, u16(b + 12), bits, false};

    if (code == 0xfffe) {
        if (extra < 22 || size < 40) bad(path, "truncated extensible fmt chunk");
        // Canonical suffix of KSDATAFORMAT_SUBTYPE_PCM/IEEE_FLOAT, serialized LE.
        constexpr std::array<Byte, 12> guid_tail = {
            0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xaa, 0x00, 0x38, 0x9b, 0x71
        };
        if (!std::equal(guid_tail.begin(), guid_tail.end(), b + 28)) {
            bad(path, "unsupported extensible WAV subformat");
        }
        code = u32(b + 24);
        fmt.valid_bits = u16(b + 18);
    }

    if (fmt.channels != 1 && fmt.channels != 2) bad(path, "only mono and stereo are supported");
    if (fmt.rate < min_rate || fmt.rate > max_rate) bad(path, "sample rate must be 8000..192000 Hz");
    switch (code) {
    case 1:
        if (bits != 8 && bits != 16 && bits != 24 && bits != 32) {
            bad(path, "PCM must have 8, 16, 24, or 32 bits per sample");
        }
        if (fmt.valid_bits == 0 || fmt.valid_bits > bits) bad(path, "invalid PCM valid bits");
        break;
    case 3:
        if ((bits != 32 && bits != 64) || fmt.valid_bits != bits) {
            bad(path, "IEEE float WAV must have 32 or 64 bits per sample");
        }
        fmt.fp = true;
        break;
    default:
        bad(path, "unsupported WAVE encoding (expected PCM or IEEE float)");
    }
    if (fmt.align != fmt.channels * fmt.width || u32(b + 8) != fmt.rate * fmt.align) {
        bad(path, "invalid block alignment or byte rate");
    }
    return fmt;
}

double decode(const Byte* p, const Format& fmt) {
    if (fmt.fp) {
        if (fmt.width == 4) {
            const std::uint32_t bits = u32(p);
            float value;
            std::memcpy(&value, &bits, sizeof(value));
            return value;
        }
        const std::uint64_t bits = u64(p);
        double value;
        std::memcpy(&value, &bits, sizeof(value));
        return value;
    }
    std::uint32_t bits = p[0];
    if (fmt.width >= 2) bits = u16(p);
    if (fmt.width >= 3) bits |= std::uint32_t(p[2]) << 16;
    if (fmt.width == 4) bits |= std::uint32_t(p[3]) << 24;
    // Extensible PCM with fewer valid bits uses left-aligned sample values.
    bits >>= 8 * fmt.width - fmt.valid_bits;
    const std::int64_t half = std::int64_t{1} << (fmt.valid_bits - 1);
    std::int64_t s = bits;
    if (fmt.width == 1) s -= half; // 8-bit PCM is unsigned
    else if (s >= half) s -= 2 * half;
    return static_cast<double>(s) / static_cast<double>(half);
}

void put16(Byte* p, std::uint32_t value) {
    p[0] = static_cast<Byte>(value);
    p[1] = static_cast<Byte>(value >> 8);
}

void put32(Byte* p, std::uint32_t value) {
    put16(p, value);
    put16(p + 2, value >> 16);
}

void write_bytes(std::ofstream& stream, const Byte* bytes, std::size_t count, const fs::path& path) {
    if (!stream.write(reinterpret_cast<const char*>(bytes), static_cast<std::streamsize>(count))) {
        bad(path, "write failed");
    }
}
} // namespace

Audio read(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary | std::ios::ate);
    if (!stream) bad(path, "could not open for reading");
    const auto physical_end = static_cast<std::streamoff>(stream.tellg());
    if (physical_end < 12) bad(path, "missing RIFF/WAVE header");
    seek(stream, 0, path);
    std::array<Byte, 12> h{};
    read_bytes(stream, h.data(), h.size(), path);
    if (!id(h.data(), "RIFF") || !id(h.data() + 8, "WAVE")) bad(path, "not a little-endian RIFF/WAVE file");
    const std::uint64_t riff_end = std::uint64_t{8} + u32(h.data() + 4);
    if (riff_end < 12 || riff_end > static_cast<std::uint64_t>(physical_end)) {
        bad(path, "invalid RIFF length");
    }

    Format fmt;
    bool have_fmt = false;
    bool have_data = false;
    std::uint64_t data_pos = 0;
    std::uint32_t data_size = 0;
    for (std::uint64_t pos = 12; pos < riff_end;) {
        if (riff_end - pos < 8) bad(path, "incomplete chunk header");
        seek(stream, pos, path);
        read_bytes(stream, h.data(), 8, path);
        const std::uint32_t size = u32(h.data() + 4);
        const std::uint64_t payload = pos + 8;
        const std::uint64_t end = payload + size + (size & 1u);
        if (end > riff_end) bad(path, "chunk or its padding exceeds RIFF bounds");
        if (id(h.data(), "fmt ")) {
            if (have_fmt) bad(path, "multiple fmt chunks");
            fmt = parse_format(stream, size, path);
            have_fmt = true;
        } else if (id(h.data(), "data")) {
            if (have_data) bad(path, "multiple data chunks");
            data_pos = payload;
            data_size = size;
            have_data = true;
        }
        pos = end;
    }
    if (!have_fmt || !have_data) bad(path, "missing fmt or data chunk");
    if (data_size % fmt.align != 0) bad(path, "data does not contain whole audio frames");
    const std::size_t frames = data_size / fmt.align;
    if (frames > max_samples) bad(path, "more than 16,000,000 mono samples");
    Audio audio{std::vector<float>(frames), static_cast<int>(fmt.rate)};
    seek(stream, data_pos, path);
    std::array<Byte, block_frames * 16> bytes;
    for (std::size_t begin = 0; begin < frames;) {
        const std::size_t count = std::min(block_frames, frames - begin);
        read_bytes(stream, bytes.data(), count * fmt.align, path);
        const Byte* p = bytes.data();
        for (std::size_t i = 0; i < count; ++i, p += fmt.align) {
            double value = decode(p, fmt);
            if (fmt.channels == 2) value = 0.5 * value + 0.5 * decode(p + fmt.width, fmt);
            if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<float>::max()) {
                bad(path, "nonfinite sample or sample outside the finite float32 range");
            }
            audio.samples[begin + i] = static_cast<float>(value);
        }
        begin += count;
    }
    return audio;
}

void write(const std::filesystem::path& path, const std::vector<float>& samples,
           int sample_rate, bool pcm16) {
    if (sample_rate < min_rate || sample_rate > max_rate) bad(path, "invalid output sample rate");
    if (samples.size() > max_samples) bad(path, "render exceeds the 16,000,000 mono sample limit");
    if (!std::all_of(samples.begin(), samples.end(), [](float v) { return std::isfinite(v); })) {
        bad(path, "render contains nonfinite samples");
    }
    const unsigned width = pcm16 ? 2u : 4u;
    const auto data_size = static_cast<std::uint32_t>(samples.size() * width);
    std::array<Byte, 58> header{};
    std::memcpy(header.data(), "RIFF", 4);
    put32(header.data() + 4, (pcm16 ? 36u : 50u) + data_size);
    std::memcpy(header.data() + 8, "WAVEfmt ", 8);
    put32(header.data() + 16, pcm16 ? 16u : 18u);
    put16(header.data() + 20, pcm16 ? 1u : 3u);
    put16(header.data() + 22, 1);
    put32(header.data() + 24, static_cast<unsigned>(sample_rate));
    put32(header.data() + 28, static_cast<unsigned>(sample_rate) * width);
    put16(header.data() + 32, width);
    put16(header.data() + 34, 8u * width);
    std::size_t data_header = 36;
    if (!pcm16) {
        // Non-PCM WAVEFORMATEX with cbSize=0 and a fact chunk.
        std::memcpy(header.data() + 38, "fact", 4);
        put32(header.data() + 42, 4);
        put32(header.data() + 46, static_cast<std::uint32_t>(samples.size()));
        data_header = 50;
    }
    std::memcpy(header.data() + data_header, "data", 4);
    put32(header.data() + data_header + 4, data_size);

    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream) bad(path, "could not open for writing");
    write_bytes(stream, header.data(), data_header + 8, path);
    std::array<Byte, block_frames * 4> bytes;
    for (std::size_t begin = 0; begin < samples.size();) {
        const std::size_t count = std::min(block_frames, samples.size() - begin);
        Byte* p = bytes.data();
        for (std::size_t i = 0; i < count; ++i, p += width) {
            const float value = samples[begin + i];
            if (pcm16) {
                const double v = std::clamp(double(value) * 32768.0, -32768.0, 32767.0);
                put16(p, static_cast<std::uint32_t>(std::lround(v)));
            } else {
                std::uint32_t bits;
                std::memcpy(&bits, &value, sizeof(bits));
                put32(p, bits);
            }
        }
        write_bytes(stream, bytes.data(), count * width, path);
        begin += count;
    }
    stream.close();
    if (!stream) bad(path, "could not finish writing output");
}

void require_distinct(const std::filesystem::path& input, const std::filesystem::path& output) {
    namespace fs = std::filesystem;
    if (input.empty() || output.empty()) throw std::runtime_error("empty input or output path");
    if (!fs::is_regular_file(fs::status(input))) {
        throw std::runtime_error("input is not a regular file: " + input.u8string());
    }
    std::error_code ec;
    const fs::file_status output_link = fs::symlink_status(output, ec);
    if (ec == std::errc::no_such_file_or_directory || output_link.type() == fs::file_type::not_found) return;
    if (ec) throw fs::filesystem_error("cannot inspect output", output, ec);
    if (!fs::is_regular_file(fs::status(output))) {
        throw std::runtime_error("output exists but is not a regular file: " + output.u8string());
    }
    if (fs::equivalent(input, output)) throw std::runtime_error("input and output refer to the same file");
}

} // namespace tn::wav
