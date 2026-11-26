// SPDX-License-Identifier: GPL-3.0-only
#include "mrq.hpp"
#include "wav.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <limits>
#include <memory>
#include <vector>

#include <fcntl.h>
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <io.h>
#include <windows.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace tn::pitch {
namespace {
using Byte = unsigned char;
constexpr std::uint64_t max_file_bytes = 384ull*1024*1024;
constexpr std::uint32_t max_entries = 1'000'000, max_frames = 1'000'000;
constexpr std::uint32_t max_name_units = 32'767;
constexpr std::uint32_t track_header = 12; // nf0, fs, nhop
constexpr std::uint64_t min_record_bytes = 4+4+track_header; // empty name, length, track
static_assert(max_file_bytes <= std::numeric_limits<long>::max());
static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559);

struct FileCloser {
    void operator()(std::FILE* file) const { std::fclose(file); }
};
using File = std::unique_ptr<std::FILE, FileCloser>;

File open_mrq(const std::filesystem::path& path, std::uint64_t& size) {
#ifdef _WIN32
    const auto handle = CreateFileW(path.c_str(), GENERIC_READ,
        FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, nullptr, OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL, nullptr);
    if (handle == INVALID_HANDLE_VALUE) return {};
    OVERLAPPED lock{};
    lock.Offset = 4;
    LARGE_INTEGER length{};
    // OpenMoresampler uses an exclusive lock on the version at bytes 4..7.
    // Readers can share it; a concurrent writer must not block rendering.
    if (GetFileType(handle) != FILE_TYPE_DISK ||
        !LockFileEx(handle, LOCKFILE_FAIL_IMMEDIATELY, 0, 4, 0, &lock) ||
        !GetFileSizeEx(handle, &length) || length.QuadPart < 0 ||
        static_cast<std::uint64_t>(length.QuadPart) > max_file_bytes) {
        CloseHandle(handle);
        return {};
    }
    size = static_cast<std::uint64_t>(length.QuadPart);
    const int fd = _open_osfhandle(reinterpret_cast<std::intptr_t>(handle), _O_RDONLY | _O_BINARY);
    if (fd == -1) { CloseHandle(handle); return {}; }
    // The stream owns the descriptor (and hence the lock) on success.
    File stream(_fdopen(fd, "rb"));
    if (!stream) _close(fd);
#else
    // Don't block on a FIFO at the expected descriptor path; use the actual
    // opened descriptor to check file kind and size.
    const int fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC | O_NONBLOCK);
    if (fd == -1) return {};
    struct stat info{};
    if (fstat(fd, &info) != 0 || !S_ISREG(info.st_mode) || info.st_size < 0 ||
        static_cast<std::uint64_t>(info.st_size) > max_file_bytes) {
        ::close(fd);
        return {};
    }
    size = static_cast<std::uint64_t>(info.st_size);
    File stream(fdopen(fd, "rb"));
    if (!stream) ::close(fd);
#endif
    return stream;
}

std::uint32_t le32(const Byte* b) {
    return std::uint32_t(b[0]) | (std::uint32_t(b[1]) << 8) |
           (std::uint32_t(b[2]) << 16) | (std::uint32_t(b[3]) << 24);
}

struct Reader {
    std::FILE* file;
    std::uint64_t remaining;

    bool read(Byte* destination, std::size_t size) {
        if (size > remaining || std::fread(destination, 1, size, file) != size) return false;
        remaining -= size;
        return true;
    }
    bool read(std::uint32_t& value) {
        Byte bytes[4];
        if (!read(bytes, sizeof(bytes))) return false;
        value = le32(bytes);
        return true;
    }
    bool skip(std::uint64_t size) {
        if (size > remaining || std::fseek(file, static_cast<long>(size), SEEK_CUR) != 0) return false;
        remaining -= size;
        return true;
    }
};

Track read_track(Reader& in, std::uint32_t size, int rate, std::uint64_t samples) {
    Byte header[track_header];
    if (!in.read(header, sizeof(header))) return {};
    const auto count = le32(header), fs = le32(header+4), hop = le32(header+8);
    if (count == 0 || count > max_frames || count > (size-track_header)/4 ||
        fs != static_cast<unsigned>(rate) || hop == 0 || hop > fs) return {};
    // The last center may be the recording endpoint. Some producers count
    // floor(N/hop) frames, while Nebula rounds the resampled 8 kHz length
    // before counting frames. Allow that boundary rounding and float32 math,
    // but don't accept a trajectory whose duration belongs to another clip.
    const std::uint64_t tolerance = 1 + (fs+15'999)/16'000 + samples/2'000'000;
    // A rounded one-frame file can be somewhat shorter than its hop, but
    // should still fill at least half the analysis interval.
    if (hop > samples+tolerance && (count != 1 || hop > 2*samples+tolerance)) return {};
    if ((std::uint64_t(count)-1)*hop > samples+tolerance ||
        (std::uint64_t(count)+1)*hop+tolerance < samples) return {};

    Track result;
    result.step = static_cast<int>(hop);
    result.frames.reserve(count);
    std::array<Byte, 4096> buffer;
    for (std::uint32_t decoded = 0; decoded < count;) {
        const auto block = std::min<std::size_t>(buffer.size()/4, count-decoded);
        if (!in.read(buffer.data(), block*4)) return {};
        for (std::size_t i = 0; i < block; ++i) {
            const auto bits = le32(buffer.data()+4*i);
            float hz;
            std::memcpy(&hz, &bits, sizeof(hz));
            if (hz == 0.f) result.frames.push_back({});
            else if (hz >= 30.f && hz <= 2400.f) result.frames.push_back({hz, 1.});
            else return {}; // including NaN and infinities
        }
        decoded += static_cast<std::uint32_t>(block);
    }
    // Some descriptors contain both layouts, and an in-place rewrite can
    // retain padding. Timestamp and producer/manual flags are not identities.
    return result;
}
} // namespace

Track read_mrq(const std::filesystem::path& mrq, const std::u16string& basename,
               int rate, std::size_t samples) {
    if (basename.empty() || basename.size() > max_name_units ||
        basename.find(u'\0') != std::u16string::npos ||
        rate < 8000 || rate > 192000 || samples == 0 || samples > wav::max_samples) return {};

    std::uint64_t size = 0;
    const auto file = open_mrq(mrq, size);
    if (!file) return {};
    Reader in{file.get(), size};
    Byte header[12];
    if (!in.read(header, sizeof(header)) || std::memcmp(header, "mrq ", 4) != 0) return {};
    const auto version = le32(header+4), entries = le32(header+8);
    if ((version != 1 && version != 2) || entries > max_entries ||
        entries > in.remaining/min_record_bytes) return {};

    std::vector<Byte> candidate_name(basename.size()*2);
    for (std::uint32_t i = 0; i < entries; ++i) {
        std::uint32_t units, payload;
        if (!in.read(units) || units > max_name_units) return {};
        bool match = units == basename.size();
        if (match) {
            if (!in.read(candidate_name.data(), candidate_name.size())) return {};
            for (std::size_t j = 0; match && j < units; ++j) {
                const auto stored = unsigned(candidate_name[2*j]) | (unsigned(candidate_name[2*j+1]) << 8);
                match = stored == basename[j]; // exact UTF-16LE; also rejects zeroed/tombstoned names
            }
        } else if (!in.skip(std::uint64_t(units)*2)) return {};
        if (!in.read(payload) || payload < track_header || payload%4 != 0 || payload > in.remaining) return {};
        if (match) return read_track(in, payload, rate, samples);
        if (!in.skip(payload)) return {};
    }
    return {};
}

namespace {
// Moresampler's own producer stamps every generated entry with this marker, and
// its readers treat 0 as "edited by hand". Matching the value matters: a file
// written here must look like one of its own to the other tool.
constexpr std::int32_t mrq_modified_by_tool = 524;
// v2 payload prefix: nf0, fs, nhop, timestamp, modified. v1 omits the last two,
// which is why an entry can only be rewritten in place when it has room for the
// v2 framing as well as the trajectory itself.
constexpr std::uint32_t mrq_v2_fixed = 20;

void put32(Byte* out, std::uint32_t value) {
    out[0] = Byte(value & 0xFF);
    out[1] = Byte((value >> 8) & 0xFF);
    out[2] = Byte((value >> 16) & 0xFF);
    out[3] = Byte((value >> 24) & 0xFF);
}

struct WriteCloser {
    void operator()(std::FILE* file) const { std::fclose(file); }
};
using WriteFile = std::unique_ptr<std::FILE, WriteCloser>;

bool write_bytes(std::FILE* file, const void* data, std::size_t size) {
    return size == 0 || std::fwrite(data, 1, size, file) == size;
}

bool write32(std::FILE* file, std::uint32_t value) {
    Byte bytes[4];
    put32(bytes, value);
    return write_bytes(file, bytes, sizeof(bytes));
}

bool write_name(std::FILE* file, const std::u16string& name) {
    std::vector<Byte> bytes(name.size()*2);
    for (std::size_t i = 0; i < name.size(); ++i) {
        bytes[2*i] = Byte(std::uint32_t(name[i]) & 0xFF);
        bytes[2*i+1] = Byte((std::uint32_t(name[i]) >> 8) & 0xFF);
    }
    return write_bytes(file, bytes.data(), bytes.size());
}

bool read32(std::FILE* file, std::uint32_t& value) {
    Byte bytes[4];
    if (std::fread(bytes, 1, sizeof(bytes), file) != sizeof(bytes)) return false;
    value = le32(bytes);
    return true;
}

// Moresampler locks the version quadword for the duration of a write and its
// readers probe the same range, so both directions stay consistent.
#ifdef _WIN32
bool lock_for_write(std::FILE* file) {
    const auto handle = reinterpret_cast<HANDLE>(_get_osfhandle(_fileno(file)));
    if (handle == INVALID_HANDLE_VALUE) return false;
    for (int attempt = 0; attempt < 300; ++attempt) {
        OVERLAPPED lock{};
        lock.Offset = 4;
        if (LockFileEx(handle, LOCKFILE_EXCLUSIVE_LOCK | LOCKFILE_FAIL_IMMEDIATELY, 0, 4, 0, &lock))
            return true;
        Sleep(50);
    }
    return false;
}
void unlock_for_write(std::FILE* file) {
    const auto handle = reinterpret_cast<HANDLE>(_get_osfhandle(_fileno(file)));
    if (handle == INVALID_HANDLE_VALUE) return;
    OVERLAPPED lock{};
    lock.Offset = 4;
    UnlockFileEx(handle, 0, 4, 0, &lock);
}
#else
bool lock_for_write(std::FILE*) { return true; }
void unlock_for_write(std::FILE*) {}
#endif

WriteFile open_for_write(const std::filesystem::path& path, bool& created) {
    created = false;
#ifdef _WIN32
    std::FILE* file = _wfopen(path.c_str(), L"r+b");
    if (file) return WriteFile(file);
    file = _wfopen(path.c_str(), L"w+b");
#else
    std::FILE* file = std::fopen(path.c_str(), "r+b");
    if (file) return WriteFile(file);
    file = std::fopen(path.c_str(), "w+b");
#endif
    created = file != nullptr;
    return WriteFile(file);
}
} // namespace

bool write_mrq(const std::filesystem::path& mrq, const std::u16string& basename,
               const Track& track, int rate) {
    if (basename.empty() || basename.size() > max_name_units ||
        basename.find(u'\0') != std::u16string::npos) return false;
    if (rate < 8000 || rate > 192000) return false;
    if (track.step <= 0 || track.step > rate || track.frames.empty() ||
        track.frames.size() > max_frames) return false;
    for (const auto& frame : track.frames) {
        // Zero is unvoiced; anything else has to be a plausible periodic pitch.
        if (!std::isfinite(frame.hz) || frame.hz < 0. || frame.hz > 2400. ||
            (frame.hz > 0. && frame.hz < 30.)) return false;
    }
    const auto count = static_cast<std::uint32_t>(track.frames.size());
    const std::uint64_t payload = std::uint64_t(mrq_v2_fixed) + std::uint64_t(count)*4;
    if (payload + 12 > max_file_bytes) return false;

    bool created = false;
    auto file = open_for_write(mrq, created);
    if (!file) return false;
    if (!lock_for_write(file.get())) return false;

    bool ok = false;
    do {
        if (created) {
            if (!write_bytes(file.get(), "mrq ", 4) || !write32(file.get(), 2) ||
                !write32(file.get(), 0)) break;
        }
        std::uint64_t total = 0;
        if (std::fseek(file.get(), 0, SEEK_END) != 0) break;
        const long end = std::ftell(file.get());
        if (end < 0) break;
        total = static_cast<std::uint64_t>(end);
        if (total < 12) break;
        if (std::fseek(file.get(), 0, SEEK_SET) != 0) break;
        Byte header[12];
        if (std::fread(header, 1, sizeof(header), file.get()) != sizeof(header)) break;
        if (std::memcmp(header, "mrq ", 4) != 0) break;
        const auto version = le32(header+4);
        if (version != 1 && version != 2) break;
        std::uint32_t entries = le32(header+8);
        if (entries > max_entries) break;

        bool matched = false;
        std::uint32_t old_payload = 0;
        std::uint64_t data_pos = 0, name_pos = 0, old_units = 0;
        bool walked = true;
        for (std::uint32_t i = 0; i < entries; ++i) {
            std::uint32_t units = 0, size = 0;
            const long name_at = std::ftell(file.get());
            if (name_at < 0 || !read32(file.get(), units) || units > max_name_units) { walked = false; break; }
            std::vector<Byte> raw(units*2);
            if (!raw.empty() && std::fread(raw.data(), 1, raw.size(), file.get()) != raw.size()) { walked = false; break; }
            const long payload_at = std::ftell(file.get());
            if (payload_at < 0 || !read32(file.get(), size) || size < track_header ||
                size % 4 != 0) { walked = false; break; }
            const long trunk_at = std::ftell(file.get());
            if (trunk_at < 0 || std::uint64_t(size) > total - std::uint64_t(trunk_at)) { walked = false; break; }
            bool same = units == basename.size();
            for (std::size_t j = 0; same && j < units; ++j) {
                const auto stored = std::uint32_t(raw[2*j]) | (std::uint32_t(raw[2*j+1]) << 8);
                same = stored == basename[j];
            }
            if (same) {
                matched = true;
                old_payload = size;
                data_pos = std::uint64_t(trunk_at);
                name_pos = std::uint64_t(name_at) + 4;
                old_units = units;
                break;
            }
            if (std::fseek(file.get(), long(size), SEEK_CUR) != 0) { walked = false; break; }
        }
        if (!walked) break;

        const bool in_place = matched && old_payload >= mrq_v2_fixed + count*4;
        if (matched && !in_place) {
            // Tombstone the obsolete record exactly as Moresampler does, then
            // append the new one; readers skip zeroed names by design.
            std::vector<Byte> zeros(old_units*2, 0);
            if (std::fseek(file.get(), long(name_pos), SEEK_SET) != 0) break;
            if (!write_bytes(file.get(), zeros.data(), zeros.size())) break;
        }
        if (!matched || !in_place) {
            // The entry count only grows when a record is appended.
            if (std::fseek(file.get(), 8, SEEK_SET) != 0) break;
            const std::uint32_t next = entries + 1;
            if (!write32(file.get(), next)) break;
            if (std::fseek(file.get(), 0, SEEK_END) != 0) break;
            if (!write32(file.get(), static_cast<std::uint32_t>(basename.size())) ||
                !write_name(file.get(), basename) ||
                !write32(file.get(), static_cast<std::uint32_t>(payload))) break;
        } else if (std::fseek(file.get(), long(data_pos), SEEK_SET) != 0) {
            break;
        }

        const auto stamp = static_cast<std::int32_t>(std::time(nullptr));
        if (!write32(file.get(), count) || !write32(file.get(), static_cast<std::uint32_t>(rate)) ||
            !write32(file.get(), static_cast<std::uint32_t>(track.step))) break;
        bool trajectory = true;
        std::array<Byte, 4096> block{};
        std::size_t filled = 0;
        for (std::uint32_t i = 0; i < count && trajectory; ++i) {
            const float hz = static_cast<float>(track.frames[i].hz > 0. ? track.frames[i].hz : 0.);
            std::uint32_t bits = 0;
            std::memcpy(&bits, &hz, sizeof(bits));
            put32(block.data() + 4*filled, bits);
            if (++filled == block.size()/4) {
                trajectory = write_bytes(file.get(), block.data(), block.size());
                filled = 0;
            }
        }
        if (trajectory && filled) trajectory = write_bytes(file.get(), block.data(), filled*4);
        if (!trajectory) break;
        if (!write32(file.get(), static_cast<std::uint32_t>(stamp)) ||
            !write32(file.get(), static_cast<std::uint32_t>(mrq_modified_by_tool))) break;
        // A v1 descriptor becomes v2 once it carries the new framing.
        if (std::fseek(file.get(), 4, SEEK_SET) != 0 || !write32(file.get(), 2)) break;
        if (std::fflush(file.get()) != 0) break;
        ok = true;
    } while (false);

    unlock_for_write(file.get());
    return ok;
}
} // namespace tn::pitch
