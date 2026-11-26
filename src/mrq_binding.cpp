// SPDX-License-Identifier: GPL-3.0-only
#include "mrq_binding.hpp"

#ifdef _WIN32
#include "cache_win.hpp"
#include "mrq.hpp"
#include "wav.hpp"

#include <array>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>

namespace tn::pitch {
namespace {
using namespace tn::detail::cache;
constexpr auto kind = L"mrq-bindings-v1";
constexpr auto setting = L"TN_FINDS_MRQ_BINDINGS";
constexpr std::array<char, 8> magic{'T','N','M','R','Q','B','1',0};
constexpr std::size_t max_path = 32767;
constexpr std::size_t max_scan = 250'000; // explicit command only, never on render
constexpr std::size_t max_depth = 24;
static_assert(sizeof(wchar_t) == 2 && sizeof(float) == 4);

std::uint32_t rot(std::uint32_t v, int n) { return (v << n) | (v >> (32-n)); }

// OpenUtau uses K4os XXH32 with seed zero on UTF-8 bytes of UOto.File.
// Hash the spelling, not a case-folded path or a resolved junction target.
std::uint32_t path_fingerprint(const fs::path& path) {
    constexpr std::uint32_t p1=2654435761u, p2=2246822519u, p3=3266489917u,
                            p4=668265263u, p5=374761393u;
    const auto text = path.u8string();
    const auto* p = reinterpret_cast<const unsigned char*>(text.data());
    const auto* end = p+text.size();
    const auto word = [](const unsigned char* p) {
        return std::uint32_t(p[0]) | (std::uint32_t(p[1])<<8) | (std::uint32_t(p[2])<<16) | (std::uint32_t(p[3])<<24);
    };
    std::uint32_t h = p5;
    if (text.size() >= 16) {
        std::uint32_t lanes[]{p1+p2,p2,0,0u-p1};
        do {
            for (auto& lane : lanes) {
                lane = rot(lane+word(p)*p2,13)*p1;
                p += 4;
            }
        } while (end-p >= 16);
        h = rot(lanes[0],1)+rot(lanes[1],7)+rot(lanes[2],12)+rot(lanes[3],18);
    }
    h += static_cast<std::uint32_t>(text.size());
    while (end-p >= 4) { h = rot(h+word(p)*p3,17)*p4; p += 4; }
    while (p != end) h = rot(h+*p++*p5,11)*p1;
    h = (h^(h>>15))*p2;
    h = (h^(h>>13))*p3;
    return h^(h>>16);
}

std::optional<std::uint32_t> alias_fingerprint(std::wstring_view name) {
    if (name.size() != 34 || name.substr(0,4) != L"src-" || name.substr(30) != L".wav" ||
        name[12] != L'-' || name[21] != L'-') return {};
    std::uint32_t value = 0;
    for (std::size_t i = 4; i < 30; ++i) {
        if (i == 12 || i == 21) continue;
        const auto c = name[i];
        unsigned digit;
        if (c >= L'0' && c <= L'9') digit = c-L'0';
        else if (c >= L'a' && c <= L'f') digit = c-L'a'+10;
        else return {};
        if (i >= 22) value = (value<<4)|digit;
    }
    return value;
}

struct Stamp {
    std::uint64_t size = 0, modified = 0;
    bool operator==(const Stamp& o) const { return size == o.size && modified == o.modified; }
};

std::optional<Stamp> stamp(const fs::path& path) {
    WIN32_FILE_ATTRIBUTE_DATA data{};
    if (!GetFileAttributesExW(path.c_str(), GetFileExInfoStandard, &data) ||
        (data.dwFileAttributes & (FILE_ATTRIBUTE_DIRECTORY|FILE_ATTRIBUTE_REPARSE_POINT|FILE_ATTRIBUTE_DEVICE))) return {};
    return Stamp{(std::uint64_t(data.nFileSizeHigh)<<32)|data.nFileSizeLow,
                 (std::uint64_t(data.ftLastWriteTime.dwHighDateTime)<<32)|data.ftLastWriteTime.dwLowDateTime};
}

bool within(const fs::path& child, const fs::path& parent) {
    auto a = child.begin();
    for (auto b = parent.begin(); b != parent.end(); ++b, ++a)
        if (a == child.end() || CompareStringOrdinal(a->c_str(),-1,b->c_str(),-1,TRUE) != CSTR_EQUAL) return false;
    return true;
}

struct Header {
    std::array<char,8> signature = magic;
    Digest key{};
    Stamp source{};
    std::uint32_t length = 0, reserved = 0;
    Digest payload{};
};
static_assert(sizeof(Header) == 96);

// A present record with an empty pathname marks a compound-key collision.
// Keep it ambiguous rather than silently letting a later registration win.
struct Record { fs::path file; Stamp source; };

Bytes utf16(const fs::path& path) { return {path.native().data(), path.native().size()*sizeof(wchar_t)}; }

std::optional<Record> read_record(const Sha256& sha, const fs::path& path, const Digest& key) {
    std::ifstream read(path, std::ios::binary|std::ios::ate);
    Header head;
    const auto size = static_cast<std::streamoff>(read.tellg());
    if (size < static_cast<std::streamoff>(sizeof(head)) ||
        size > static_cast<std::streamoff>(sizeof(head)+max_path*sizeof(wchar_t))) return {};
    read.seekg(0);
    if (!read.read(reinterpret_cast<char*>(&head),sizeof(head)) || head.signature != magic || head.key != key ||
        head.reserved || head.length > max_path ||
        sizeof(head)+head.length*sizeof(wchar_t) != static_cast<std::size_t>(size)) return {};
    std::wstring source(head.length,L'\0');
    const auto len = source.size()*sizeof(wchar_t);
    if (len && !read.read(reinterpret_cast<char*>(source.data()),len)) return {};
    Digest payload;
    if (source.find(L'\0') != source.npos || !sha.hash(payload,{bytes(head.source),{source.data(),len}}) ||
        payload != head.payload) return {};
    fs::path binding(source);
    if (!binding.empty() && !binding.is_absolute()) return {};
    return Record{std::move(binding),head.source};
}

bool save_record(const Sha256& sha, const fs::path& path, const Digest& key, Record reg) {
    Header head;
    head.key = key;
    head.source = reg.source;
    head.length = static_cast<std::uint32_t>(reg.file.native().size());
    if (!sha.hash(head.payload,{bytes(head.source),utf16(reg.file)})) return false;
    publish(path,{bytes(head),utf16(reg.file)});
    const auto stored = read_record(sha,path,key);
    return stored && stored->file.native() == reg.file.native() && stored->source == reg.source;
}

bool make_key(const Sha256& sha, Digest& digest, std::uint32_t fingerprint, const std::vector<float>& x, int rate) {
    return sha.hash(digest,{bytes(fingerprint),bytes(rate),bytes(x)});
}

fs::path bucket(const fs::path& dir, std::uint32_t fingerprint) {
    constexpr wchar_t digits[] = L"0123456789abcdef";
    std::wstring name(8,L'0');
    for (int i = 7; i >= 0; --i, fingerprint >>= 4) name[i] = digits[fingerprint&15];
    return dir/name;
}

bool is_wav(const fs::path& p) {
    return CompareStringOrdinal(p.extension().c_str(),-1,L".wav",-1,TRUE) == CSTR_EQUAL;
}

bool bind_one(const Sha256& sha, const fs::path& dir, const fs::path& src) {
    if (!is_wav(src) || src.native().size() > max_path) return false;
    const auto mrq = src.parent_path()/L"desc.mrq";
    if (!stamp(mrq)) return false;
    const auto before = stamp(src);
    if (!before) return false;
    wav::Audio audio;
    try { audio = wav::read(src); }
    catch (const std::runtime_error&) { return false; }
    if (read_mrq(mrq,src.filename().u16string(),audio.sample_rate,audio.samples.size()).frames.empty()) return false;
    const auto after = stamp(src);
    if (!after || !(*before == *after)) return false;
    const auto fingerprint = path_fingerprint(src);
    Digest key;
    if (!make_key(sha,key,fingerprint,audio.samples,audio.sample_rate)) return false;
    const auto name = hex(key);
    const auto shard = bucket(dir,fingerprint);
    const auto path = shard/(name+L".mrqb");
    const EntryLock lock(L"Local\\tn_finds_mrq_bind_"+name);
    if (!lock.owned) throw std::runtime_error("could not lock MRQ binding");
    fs::create_directories(shard);
    const auto attrs = GetFileAttributesW(shard.c_str());
    if (attrs == INVALID_FILE_ATTRIBUTES || (attrs & FILE_ATTRIBUTE_REPARSE_POINT))
        throw std::runtime_error("MRQ binding bucket is unavailable or redirected");
    const auto prior = read_record(sha,path,key);
    if (prior && prior->file.empty()) return false;
    const auto matches = !prior || prior->file.native() == src.native();
    if (!save_record(sha,path,key,matches ? Record{src,*after} : Record{}))
        throw std::runtime_error("could not publish MRQ binding");
    return matches;
}
} // namespace

Track lookup_bound_mrq(const fs::path& input, const std::vector<float>& samples, int rate) {
    const auto fingerprint = alias_fingerprint(input.filename().native());
    if (!fingerprint || samples.empty()) return {};
    try {
        const auto dir = directory(setting,kind);
        if (dir.empty()) return {};
        const auto shard = bucket(dir,*fingerprint);
        const auto attrs = GetFileAttributesW(shard.c_str());
        if (attrs == INVALID_FILE_ATTRIBUTES || !(attrs & FILE_ATTRIBUTE_DIRECTORY) ||
            (attrs & FILE_ATTRIBUTE_REPARSE_POINT)) return {};
        const Sha256 sha;
        Digest key;
        if (!make_key(sha,key,*fingerprint,samples,rate)) return {};
        const auto found = read_record(sha,shard/(hex(key)+L".mrqb"),key);
        if (!found || found->file.empty() || path_fingerprint(found->file) != *fingerprint) return {};
        const auto source = stamp(found->file);
        if (!source || !(*source == found->source)) return {};
        const auto mrq = found->file.parent_path()/L"desc.mrq";
        if (!stamp(mrq)) return {};
        auto track = read_mrq(mrq,found->file.filename().u16string(),rate,samples.size());
        const auto after = stamp(found->file);
        return after && *after == *source ? std::move(track) : Track{};
    } catch (const std::runtime_error&) {
        return {}; // An unavailable binding is not a render failure.
    }
}

std::size_t bind_mrq_sources(const fs::path& original) {
    const auto dir = directory(setting,kind);
    if (dir.empty()) throw std::runtime_error("MRQ binding requires a writable local cache root");
    auto src = fs::absolute(original).make_preferred();
    const auto attrs = GetFileAttributesW(src.c_str());
    if (attrs == INVALID_FILE_ATTRIBUTES || (attrs & FILE_ATTRIBUTE_REPARSE_POINT))
        throw std::invalid_argument("select an existing original WAV or voicebank directory (not a link)");
    const bool bank = (attrs & FILE_ATTRIBUTE_DIRECTORY) != 0;
    const auto bank_dir = fs::canonical(bank ? src : src.parent_path());
    if (bank_dir == bank_dir.root_path() || within(fs::weakly_canonical(dir),bank_dir))
        throw std::invalid_argument("select a voicebank below its drive root and keep MRQ bindings outside it");
    fs::create_directories(dir);
    const Sha256 sha;
    if (!bank) return bind_one(sha,dir,src) ? 1 : 0;
    std::size_t count = 0, checked = 0;
    std::error_code ec;
    fs::recursive_directory_iterator it(src,fs::directory_options::skip_permission_denied,ec), end;
    for (; it != end; it.increment(ec)) {
        if (++checked > max_scan) throw std::runtime_error("MRQ registration exceeded directory scan limit; select a smaller bank");
        const auto& candidate = it->path();
        const auto a = GetFileAttributesW(candidate.c_str());
        if (a == INVALID_FILE_ATTRIBUTES || (a & FILE_ATTRIBUTE_REPARSE_POINT) ||
            static_cast<std::size_t>(it.depth()) >= max_depth) it.disable_recursion_pending();
        if (a != INVALID_FILE_ATTRIBUTES && !(a & (FILE_ATTRIBUTE_DIRECTORY|FILE_ATTRIBUTE_REPARSE_POINT)) &&
            bind_one(sha,dir,candidate)) ++count;
    }
    if (ec) throw fs::filesystem_error("cannot scan selected MRQ voicebank",src,ec);
    return count;
}
} // namespace tn::pitch

#else
namespace tn::pitch {
Track lookup_bound_mrq(const std::filesystem::path&, const std::vector<float>&, int) { return {}; }
std::size_t bind_mrq_sources(const std::filesystem::path&) { return 0; }
} // namespace tn::pitch
#endif
