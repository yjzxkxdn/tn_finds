// SPDX-License-Identifier: GPL-3.0-only
#include "wavtool.hpp"
#include "fractional_sinc.hpp"
#include "wav.hpp"
#include "bat_host.hpp"
#include "phrase.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace tn::wavtool {
namespace {
namespace fs = std::filesystem;
using Byte = unsigned char;
constexpr int sample_limit = static_cast<int>(wav::max_samples);
constexpr std::size_t io_frames = 4096;
constexpr std::size_t header_size = 60;
using Header = std::array<Byte, header_size>;

[[noreturn]] void io_error(const fs::path& path, const char* reason) {
    throw std::runtime_error("wavtool '" + path.u8string() + "': " + reason);
}

bool number(std::string_view text, double& value) {
    if (text.empty() || text.size() > 128) return false;
    std::istringstream input{std::string(text)};
    input.imbue(std::locale::classic());
    input >> std::noskipws >> value;
    return input && input.eof() && std::isfinite(value);
}

double numeric(std::string_view text, const char* label) {
    double value = 0;
    if (!number(text, value) || std::abs(value) > 1.e7)
        throw std::invalid_argument(std::string("invalid wavtool ") + label);
    return value;
}

double duration(std::string_view text) {
    const auto sep = text.find('@');
    if (sep == text.npos) return numeric(text, "duration");
    const double ticks = numeric(text.substr(0, sep), "tick length");
    const auto suffix = text.substr(sep+1);
    if (suffix.empty() || suffix.size() > 128)
        throw std::invalid_argument("wavtool duration expects ticks@BPM[+/-milliseconds]");
    std::istringstream in{std::string(suffix)};
    in.imbue(std::locale::classic());
    double bpm = 0;
    in >> std::noskipws >> bpm;
    if (!in || bpm < 1. || bpm > 1000. || !std::isfinite(bpm) || ticks < 0.)
        throw std::invalid_argument("invalid wavtool tempo or tick length");
    double correction = 0;
    if (!in.eof()) {
        if ((in.peek() != '+' && in.peek() != '-') || !(in >> correction) || !in.eof() ||
            !std::isfinite(correction) || std::abs(correction) > 1.e7)
            throw std::invalid_argument("invalid wavtool duration correction");
    }
    return ticks * 125. / bpm + correction;
}

struct Note {
    fs::path output, input;
    double offset = 0, length = 0, overlap = 0;
    bool last = false, rest = false;
    std::array<double, 7> positions{}, gains{};
};

Note parse(const std::vector<fs::path>& args, std::size_t first) {
    std::size_t end = args.size();
    const bool last = end > first && args.back() == "LAST_NOTE";
    if (last) --end;
    if (end < first+4 || end > first+15)
        throw std::invalid_argument("wavtool requires output, input, offset, duration and up to 11 envelope values");
    Note note;
    note.output = args[first];
    note.input = args[first+1];
    note.last = last;
    if (note.input.empty() || note.output.empty()) throw std::invalid_argument("empty wavtool path");
    note.offset = numeric(args[first+2].u8string(), "offset");
    note.length = duration(args[first+3].u8string());
    if (note.length < 0. || note.length * frames_per_ms > sample_limit)
        throw std::invalid_argument("wavtool duration is outside the supported range");

    std::array<double, 11> env{0., 0., 0., 100., 100., 100., 100., 0., 0., 0., 100.};
    const auto fields = end-first-4;
    for (std::size_t i = 0; i < fields; ++i) env[i] = numeric(args[first+4+i].u8string(), "envelope value");
    if (fields < 11) env[10] = env[4];
    note.rest = fields == 2 && env[0] == 0. && env[1] == 0.;
    if (!note.rest && args[first+1].filename().u8string() == "R.wav") {
        std::error_code ec;
        note.rest = !fs::exists(note.input, ec) && (!ec || ec == std::errc::no_such_file_or_directory);
    }
    note.overlap = env[7];
    for (int i : {0, 1, 2, 8, 9}) if (env[i] < 0.) throw std::invalid_argument("negative wavtool envelope time");
    for (int i : {3, 4, 5, 6, 10})
        if (env[i] < 0. || env[i] > 1000.) throw std::invalid_argument("wavtool envelope gain must be 0..1000");

    // Keep the five vertices ordered even if a short note cannot fit both ramps.
    const double span = env[0]+env[1]+env[2]+env[8]+env[9];
    const double scale = span > note.length ? note.length/span : 1.;
    const double p1 = scale*env[0], p2 = p1+scale*env[1];
    const double p4 = note.length-scale*env[8];
    const double p3 = std::max(p2, p4-scale*env[2]);
    const double p5 = std::min(p3, p2+scale*env[9]);
    note.positions = {0., p1, p2, p5, p3, p4, note.length};
    note.gains = {0., env[3]*.01, env[4]*.01, env[10]*.01, env[5]*.01, env[6]*.01, 0.};
    return note;
}

void put(Byte* at, std::uint64_t value, int size) {
    for (int i = 0; i < size; ++i) at[i] = static_cast<Byte>(value >> (8*i));
}

std::uint64_t get(const Byte* at, int size) {
    std::uint64_t out = 0;
    for (int i = 0; i < size; ++i) out |= std::uint64_t(at[i]) << (8*i);
    return out;
}

Header wav_header(double cursor, int frames) {
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559,
                  "wavtool cursors use IEEE binary64");
    Header h{};
    std::memcpy(h.data(), "RIFF", 4);
    put(h.data()+4, header_size-8+2u*frames, 4);
    std::memcpy(h.data()+8, "WAVEfmt ", 8);
    put(h.data()+16, 16, 4);
    put(h.data()+20, 1, 2);
    put(h.data()+22, 1, 2);
    put(h.data()+24, rate, 4);
    put(h.data()+28, 2*rate, 4);
    put(h.data()+32, 2, 2);
    put(h.data()+34, 16, 2);
    // Carry the fractional sample clock between processes in a WAV chunk
    // whose payload is ignored on playback. This bounds rounding drift.
    std::memcpy(h.data()+36, "JUNK", 4);
    put(h.data()+40, 8, 4);
    std::uint64_t bits;
    std::memcpy(&bits, &cursor, 8);
    put(h.data()+44, bits, 8);
    std::memcpy(h.data()+52, "data", 4);
    put(h.data()+56, 2u*frames, 4);
    return h;
}

bool exists_regular(const fs::path& path) {
    std::error_code ec;
    const auto status = fs::symlink_status(path, ec);
    if (status.type() == fs::file_type::not_found || ec == std::errc::no_such_file_or_directory) return false;
    if (ec) throw fs::filesystem_error("could not inspect wavtool output", path, ec);
    if (!fs::is_regular_file(status)) io_error(path, "output must be a regular file");
    return true;
}

void ensure_no_alias(const fs::path& a, const fs::path& b, bool a_exists, bool b_exists) {
    if (fs::weakly_canonical(a) == fs::weakly_canonical(b) || (a_exists && b_exists && fs::equivalent(a, b)))
        throw std::invalid_argument("wavtool paths refer to the same file");
}

bool has_marker(const fs::path& path) {
    if (fs::file_size(path) != 44) return false;
    std::array<Byte, 44> marker{};
    std::ifstream in(path, std::ios::binary);
    if (!in.read(reinterpret_cast<char*>(marker.data()), marker.size())) io_error(path, "could not read phrase marker");
    return std::memcmp(marker.data(), "RIFF", 4) == 0 && get(marker.data()+4, 4) == 36 &&
           std::memcmp(marker.data()+8, "WAVEfmt ", 8) == 0 &&
           std::memcmp(marker.data()+36, "data", 4) == 0 && get(marker.data()+40, 4) == 0;
}

struct Cursor {
    double exact = 0;
    int frames = 0;
};

Cursor previous_end(const fs::path& header, const fs::path& data) {
    Header h{};
    if (fs::file_size(header) != h.size()) io_error(header, "invalid phrase header");
    std::ifstream in(header, std::ios::binary);
    if (!in.read(reinterpret_cast<char*>(h.data()), h.size())) io_error(header, "could not read phrase header");
    auto bits = get(h.data()+44, 8);
    Cursor last;
    std::memcpy(&last.exact, &bits, 8);
    if (!std::isfinite(last.exact) || last.exact < 0. || last.exact > sample_limit)
        io_error(header, "invalid phrase length");
    last.frames = static_cast<int>(std::lround(last.exact));
    if (h != wav_header(last.exact, last.frames) || fs::file_size(data) != 2u*last.frames)
        io_error(header, "phrase header and data disagree");
    return last;
}

std::vector<float> envelope(const Note& note, const wav::Audio& audio, int skipped, int size) {
    std::vector<float> output(size);
    if (audio.samples.empty() || size == 0) return output;
    const double ratio = double(audio.sample_rate)/rate;
    const double offset = std::round(note.offset*audio.sample_rate/1000.);
    const int available = static_cast<int>(audio.samples.size());
    std::optional<detail::Sinc> sinc;
    if (audio.sample_rate != rate) sinc.emplace(audio.samples, ratio);
    std::size_t k = 0;
    for (int j = 0; j < size; ++j) {
        const double position = offset+(skipped+j)*ratio;
        if (position < 0. || position >= available) continue;
        const double t = (skipped+j)/frames_per_ms;
        while (k+2 < note.positions.size() && t >= note.positions[k+1]) ++k;
        const double width = note.positions[k+1]-note.positions[k];
        const double mix = width <= 0. ? 0. : std::clamp((t-note.positions[k])/width, 0., 1.);
        const double gain = note.gains[k]+mix*(note.gains[k+1]-note.gains[k]);
        const float sample = sinc ? sinc->read(position, 0, available) : audio.samples[static_cast<int>(position)];
        output[j] = static_cast<float>(gain*sample);
    }
    return output;
}

void write_bytes(std::ostream& stream, const Byte* bytes, std::size_t size, const fs::path& path) {
    if (!stream.write(reinterpret_cast<const char*>(bytes), static_cast<std::streamsize>(size)))
        io_error(path, "write failed");
}

void finish(std::ofstream& stream, const fs::path& path) {
    stream.close();
    if (!stream) io_error(path, "could not finish writing");
}

void update_tail(const fs::path& data, bool active, int old_count, int start,
                 const std::vector<float>& incoming, int end) {
    if (!active) {
        std::ofstream fresh(data, std::ios::binary | std::ios::trunc);
        if (!fresh) io_error(data, "could not start phrase data");
        finish(fresh, data);
    }
    std::fstream stream(data, std::ios::in | std::ios::out | std::ios::binary);
    if (!stream) io_error(data, "could not open phrase data");
    const int overlap = std::max(0, std::min(old_count, end)-start);
    std::vector<float> mixed;
    std::array<Byte, io_frames*2> bytes{};
    if (overlap > 0) {
        mixed.resize(overlap);
        stream.seekg(2LL*start);
        for (int pos = 0; pos < overlap;) {
            const int count = std::min(static_cast<int>(io_frames), overlap-pos);
            if (!stream.read(reinterpret_cast<char*>(bytes.data()), 2LL*count))
                io_error(data, "could not read the previous overlap");
            for (int i = 0; i < count; ++i) {
                int value = static_cast<int>(get(bytes.data()+2*i, 2));
                if (value >= 32768) value -= 65536;
                mixed[pos+i] = static_cast<float>(value/32768.)+incoming[pos+i];
            }
            pos += count;
        }
    }

    // Rewrite only the overlap, then stream the new material or silent gap.
    stream.seekp(2LL*std::min(old_count, start));
    for (int pos = std::min(old_count, start); pos < end;) {
        const int count = std::min(static_cast<int>(io_frames), end-pos);
        for (int i = 0; i < count; ++i) {
            const int j = pos+i-start;
            const double sample = j < 0 ? 0. : j < overlap ? mixed[j] : incoming[j];
            if (!std::isfinite(sample)) io_error(data, "mix overflowed");
            const double scaled = std::clamp(sample*32768., -32768., 32767.);
            put(bytes.data()+2*i, static_cast<std::uint16_t>(std::lround(scaled)), 2);
        }
        write_bytes(stream, bytes.data(), 2u*count, data);
        pos += count;
    }
    stream.close();
    if (!stream) io_error(data, "could not finish writing phrase data");
    if (end < old_count) fs::resize_file(data, 2u*end);
}

void assemble(const fs::path& output, const fs::path& data, const Header& header, int frames) {
    std::ifstream source(data, std::ios::binary);
    if (!source) io_error(data, "could not open phrase data to finalize");
    std::ofstream dest(output, std::ios::binary | std::ios::trunc);
    if (!dest) io_error(output, "could not create final WAV");
    write_bytes(dest, header.data(), header.size(), output);
    std::array<Byte, io_frames*2> bytes;
    for (int left = frames; left > 0;) {
        const int n = std::min(left, static_cast<int>(io_frames));
        if (!source.read(reinterpret_cast<char*>(bytes.data()), 2LL*n)) io_error(data, "truncated phrase data");
        write_bytes(dest, bytes.data(), 2u*n, output);
        left -= n;
    }
    finish(dest, output);
}
bool phrase_marker(const fs::path& path) {
    if (fs::file_size(path) != header_size) return false;
    Header bytes{};
    std::ifstream stream(path, std::ios::binary);
    return stream.read(reinterpret_cast<char*>(bytes.data()), bytes.size()) && bytes == wav_header(-1., 0);
}

Clip load_clip(const Note& note, int left, int end, bool host_batch) {
    Clip clip{left, end, note.offset, note.positions, note.gains, {}};
    if (note.rest) clip.note.sample_rate = rate;
    else if (!read_model(note.input, clip.note)) {
        if (host_batch && is_dual_tool_cache(note.input))
            io_error(note.input, "missing native model for this dual-role note; clear the host note cache and render again");
        auto audio = wav::read(note.input);
        clip.note.sample_rate = audio.sample_rate;
        clip.note.active = static_cast<int>(audio.samples.size());
        clip.note.residual = std::move(audio.samples);
    }
    return clip;
}

void finish_phrase(const Note& note, const fs::path& journal_path, const fs::path& header_path,
                   const fs::path& data_path, bool active, int expected) {
    const Journal journal = active ? read_journal(journal_path) : Journal{};
    if (active && expected != journal.expected && !note.last)
        io_error(journal_path, "classic phrase call count changed");
    const double start = active ? journal.end-note.overlap*frames_per_ms : 0.;
    const double end = std::max(0., start+note.length*frames_per_ms);
    if (end > sample_limit) throw std::invalid_argument("assembled wavtool output is too long");
    const int frames = static_cast<int>(std::lround(end));
    auto clip = load_clip(note, static_cast<int>(std::lround(start)), frames, journal.expected || expected);
    const bool final = note.last || (expected && journal.count+1 == expected);
    if (!final) {
        save_clip(journal_path, journal, expected, end, clip);
        if (!active) {
            // Standard batches delete the output marker at startup. Discard
            // split fragments from an earlier aborted render before publishing
            // the new marker, so its footer cannot pick up an old result.
            fs::remove(header_path);
            fs::remove(data_path);
            const auto h = wav_header(-1., 0);
            std::ofstream marker(note.output, std::ios::binary | std::ios::trunc);
            write_bytes(marker, h.data(), h.size(), note.output);
            finish(marker, note.output);
        }
        return;
    }
    auto clips = read_clips(journal_path, journal);
    clips.push_back(std::move(clip));
    const auto result = synthesize_phrase(clips, frames);
    if (!std::all_of(result.begin(), result.end(), [](float v) { return std::isfinite(v); }))
        io_error(note.output, "mix overflowed");
    const auto head = wav_header(end, frames);
    std::ofstream stream(header_path, std::ios::binary | std::ios::trunc);
    if (!stream) io_error(header_path, "could not open phrase header for writing");
    update_tail(data_path, false, 0, 0, result, frames);
    write_bytes(stream, head.data(), head.size(), header_path);
    finish(stream, header_path);
    if (note.last) assemble(note.output, data_path, head, frames);
    if (active) {
        std::error_code ec;
        fs::remove(journal_path, ec);
    }
}
} // namespace

bool is_classic_call(const std::vector<fs::path>& args) {
    double offset = 0;
    return args.size() >= 4 && number(args[2].u8string(), offset);
}

void append(const std::vector<fs::path>& args, std::size_t first, JoinMode mode) {
    const Note note = parse(args, first);
    auto header_path = note.output, data_path = note.output;
    header_path += ".whd";
    data_path += ".dat";
    const bool out = exists_regular(note.output), whd = exists_regular(header_path), dat = exists_regular(data_path);
    ensure_no_alias(note.output, header_path, out, whd);
    ensure_no_alias(note.output, data_path, out, dat);
    ensure_no_alias(header_path, data_path, whd, dat);
    for (const auto& target : {note.output, header_path, data_path}) {
        if (fs::is_regular_file(note.input)) wav::require_distinct(note.input, target);
        else if (fs::weakly_canonical(note.input) == fs::weakly_canonical(target))
            throw std::invalid_argument("wavtool input and output paths coincide");
    }

    const bool pending_phrase = out && phrase_marker(note.output);
    if (mode != JoinMode::pcm) {
        const int expected = mode == JoinMode::classic_phrase ? classic_calls(note.output) : 0;
        if (mode == JoinMode::classic_phrase && !expected && !note.last)
            io_error(note.output, "could not find a supported classic batch footer; use --phrase with LAST_NOTE");
        if (out && has_marker(note.output)) io_error(note.output, "cannot switch a live PCM join to phrase mode");
        auto journal_path = note.output;
        journal_path += ".tnphrase";
        const bool journal = exists_regular(journal_path);
        for (const auto& target : {note.output, header_path, data_path})
            ensure_no_alias(journal_path, target, journal, fs::exists(target));
        if (fs::is_regular_file(note.input)) wav::require_distinct(note.input, journal_path);
        else if (fs::weakly_canonical(note.input) == fs::weakly_canonical(journal_path))
            throw std::invalid_argument("wavtool input and journal paths coincide");
        finish_phrase(note, journal_path, header_path, data_path, pending_phrase, expected);
        return;
    }
    if (pending_phrase) io_error(note.output, "cannot switch a live phrase to PCM wavtool mode");
    const bool active = out && has_marker(note.output);
    if (active && (!whd || !dat)) io_error(note.output, "phrase data is missing");
    const Cursor old = active ? previous_end(header_path, data_path) : Cursor{};
    const double start_exact = active ? old.exact-note.overlap*frames_per_ms : 0.;
    const double end_exact = std::max(0., start_exact+note.length*frames_per_ms);
    if (end_exact > sample_limit) throw std::invalid_argument("assembled wavtool output is too long");
    const int left = static_cast<int>(std::lround(start_exact));
    const int start = std::max(0, left);
    const int end = static_cast<int>(std::lround(end_exact));
    // This command only accepts audio. Native models are consumed by the
    // phrase renderer, so selecting the wrong tool cannot render them again.
    const wav::Audio audio = note.rest ? wav::Audio{{}, rate} : wav::read(note.input);
    const auto input = envelope(note, audio, start-left, end-start);
    if (!std::all_of(input.begin(), input.end(), [](float v) { return std::isfinite(v); }))
        io_error(note.input, "enveloped audio exceeds the finite float32 range");

    const Header header = wav_header(end_exact, end);
    // Open an active header for update before touching its payload. A
    // failed open leaves the previous pair intact, as does invalid input.
    const auto header_mode = std::ios::binary | std::ios::out |
                             (active ? std::ios::in : std::ios::trunc);
    std::ofstream output_header(header_path, header_mode);
    if (!output_header) io_error(header_path, "could not open phrase header for writing");
    update_tail(data_path, active, old.frames, start, input, end);
    write_bytes(output_header, header.data(), header.size(), header_path);
    finish(output_header, header_path);
    if (note.last) assemble(note.output, data_path, header, end);
    else if (!active) wav::write(note.output, {}, rate, true);
}
} // namespace tn::wavtool
