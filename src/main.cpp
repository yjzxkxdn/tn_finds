// SPDX-License-Identifier: GPL-3.0-only
#include "engine.hpp"
#include "wav.hpp"

#ifndef TN_FINDS_DUAL_ROLE
#define TN_FINDS_DUAL_ROLE 0
#endif

#if TN_FINDS_DUAL_ROLE
#include "phrase.hpp"
#include "wavtool.hpp"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifndef TN_FINDS_VERSION
#define TN_FINDS_VERSION "0.4.7"
#endif

namespace {
namespace fs = std::filesystem;
using Args = std::vector<fs::path>;
constexpr std::size_t max_bends = 1'000'000;
constexpr std::size_t max_f0_bytes = 64u * 1024u * 1024u;
constexpr std::size_t max_csv_line = 1024;

struct Job {
    fs::path input, output, f0_csv;
    tn::Options options;
    bool pcm16 = false;
    bool report = false;
};

[[noreturn]] void bad(const std::string& text) { throw std::invalid_argument(text); }

bool space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }
bool digit(char c) { return c >= '0' && c <= '9'; }
bool letter(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }

std::string_view trim(std::string_view s) {
    while (!s.empty() && space(s.front())) s.remove_prefix(1);
    while (!s.empty() && space(s.back())) s.remove_suffix(1);
    return s;
}

bool parse_number(std::string_view text, double& value) {
    if (text.empty() || text.size() > 128) return false;
    std::istringstream in{std::string(text)};
    in.imbue(std::locale::classic());
    in >> std::noskipws >> value;
    return in && in.eof() && std::isfinite(value);
}

double number(std::string_view s, std::string_view label) {
    double value = 0;
    if (!parse_number(s, value)) bad("invalid numeric value for " + std::string(label));
    return value;
}

double note_hz(std::string_view note) {
    if (note.size() < 2) bad("invalid note; expected e.g. C#4 or Bb3");
    const char name = note[0] >= 'a' && note[0] <= 'g' ? static_cast<char>(note[0] - 'a' + 'A') : note[0];
    // Pitch class for A..G in a MIDI C-based octave.
    constexpr std::array<int, 7> pitch_class{9, 11, 0, 2, 4, 5, 7};
    if (name < 'A' || name > 'G') bad("invalid note; expected e.g. C#4 or Bb3");
    int key = pitch_class[static_cast<std::size_t>(name - 'A')];
    note.remove_prefix(1);
    while (!note.empty() && (note.front() == '#' || note.front() == 'b')) {
        key += note.front() == '#' ? 1 : -1;
        if (key < -12 || key > 23) bad("invalid note accidental");
        note.remove_prefix(1);
    }
    int sign = 1;
    if (!note.empty() && (note.front() == '-' || note.front() == '+')) {
        if (note.front() == '-') sign = -1;
        note.remove_prefix(1);
    }
    if (note.empty() || note.size() > 2) bad("invalid note octave");
    int octave = 0;
    for (char c : note) {
        if (!digit(c)) bad("invalid note octave");
        octave = 10 * octave + c - '0';
    }
    const int midi = (octave * sign + 1) * 12 + key;
    if (midi < 0 || midi > 127) bad("note is outside MIDI range 0..127");
    return 440. * std::exp2((midi - 69) / 12.);
}

int b64(char c) {
    if (c >= 'A' && c <= 'Z') return c - 'A';
    if (c >= 'a' && c <= 'z') return c - 'a' + 26;
    if (digit(c)) return c - '0' + 52;
    if (c == '+') return 62;
    if (c == '/') return 63;
    bad("invalid Base64 pitch bend character");
}

std::vector<int> pitch_bends(std::string_view text) {
    std::vector<int> result;
    for (std::size_t i = 0; i < text.size();) {
        if (text[i] == '#') {
            if (result.empty()) bad("pitch repeat requires a preceding value");
            ++i;
            if (i == text.size() || !digit(text[i])) bad("invalid #n# pitch repeat");
            std::size_t count = 0;
            const std::size_t left = max_bends - result.size();
            while (i < text.size() && digit(text[i])) {
                const auto d = static_cast<std::size_t>(text[i++] - '0');
                if (count > left / 10 || (count == left / 10 && d > left % 10)) {
                    bad("pitch bend exceeds 1,000,000 values");
                }
                count = 10 * count + d;
            }
            if (i == text.size() || text[i] != '#') bad("unterminated #n# pitch repeat");
            ++i;
            const int last = result.back();
            result.insert(result.end(), count, last);
        } else {
            if (text.size() - i < 2) bad("incomplete Base64 pitch bend pair");
            if (result.size() == max_bends) bad("pitch bend exceeds 1,000,000 values");
            const int code = 64 * b64(text[i]) + b64(text[i+1]);
            result.push_back(code >= 2048 ? code - 4096 : code);
            i += 2;
        }
    }
    return result;
}

double forced_hz(double value) {
    if (value < 0.) return -1.;
    if (value < 50.) return 0.;
    return std::min(value, 1000.);
}

struct VoiceFlag {
    std::string_view classic, native;
    double tn::Options::*member;
};

constexpr VoiceFlag voice_flags[] = {
    {"Mt", "--tension", &tn::Options::tension},
    {"Mb", "--breathiness", &tn::Options::breathiness},
    {"Mo", "--mouth", &tn::Options::mouth},
    {"Mr", "--resonance", &tn::Options::resonance},
    {"ME", "--formant-emphasis", &tn::Options::formant_emphasis},
};

bool bare_flag(char c) { return c == 'e' || c == 'N' || c == 'u'; }

bool known_after_me(std::string_view suffix) {
    // Me has no value and can be adjacent to another flag. Do not split an
    // unrelated name like Metal80, or expose its trailing single-letter flags.
    while (!suffix.empty()) {
        if (suffix.substr(0, 2) == "Me") suffix.remove_prefix(2);
        else if (bare_flag(suffix.front())) suffix.remove_prefix(1);
        else break;
    }
    if (suffix.empty() || suffix == "Tnr") return true;
    // Reserve the other Moresampler names without implementing their effects.
    if (suffix == "Md" || suffix == "MC" || suffix == "MG" || suffix == "MD" || suffix == "Ms" || suffix == "Mm")
        return true;
    // Keep MeB50 valid in older projects while ignoring B50 itself.
    if (suffix == "B" || (suffix.size() == 1 && std::string_view("gtbAOW").find(suffix.front()) != suffix.npos)) return true;
    return std::any_of(std::begin(voice_flags), std::end(voice_flags),
                       [&](const VoiceFlag& f) { return suffix == f.classic; });
}

void flags(std::string_view text, tn::Options& out) {
    // Values are concatenated in classic commands. Consume whole named flags
    // so that an unsupported multi-letter token isn't mistaken for a control.
    for (std::size_t i = 0; i < text.size();) {
        if (!letter(text[i]) && text[i] != '/') { ++i; continue; }
        const auto start = i++;
        const char flag = text[start];
        if (!bare_flag(flag)) {
            while (i < text.size() && (letter(text[i]) || text[i] == '/')) ++i;
            const auto word = text.substr(start, i - start);
            if (word.size() > 2 && word.substr(0, 2) == "Me" && known_after_me(word.substr(2))) i = start + 2;
        }
        const auto key = text.substr(start, i - start);
        const auto begin = i;
        if (i < text.size() && (text[i] == '+' || text[i] == '-')) ++i;
        const auto digits_begin = i;
        while (i < text.size() && digit(text[i])) ++i;
        bool any_digit = i != digits_begin;
        if (i < text.size() && text[i] == '.') {
            const auto dot_end = ++i;
            while (i < text.size() && digit(text[i])) ++i;
            any_digit = any_digit || i != dot_end;
        }
        if (key == "e" && !any_digit) { out.loop = out.force_loop = false; continue; }
        if (key == "Me" && !any_digit) { out.loop = out.force_loop = true; continue; }
        // OpenUtau's Options expressions append their option text verbatim.
        // e0/e1 and Me0/Me1 keep that choice separate from a numeric flag.
        if (key == "e" || key == "Me") {
            const double option = number(text.substr(begin, i - begin), "UTAU flag " + std::string(key));
            if (option == 1.) out.loop = out.force_loop = key == "Me";
            continue;
        }
        if (!any_digit) continue;
        const auto voice = std::find_if(std::begin(voice_flags), std::end(voice_flags),
                                        [&](const VoiceFlag& f) { return f.classic == key; });
        if (voice != std::end(voice_flags)) {
            out.*(voice->member) = std::clamp(number(text.substr(begin, i - begin),
                                                    "UTAU flag " + std::string(key)), -100., 100.);
            continue;
        }
        if (key == "Tnr") {
            out.mask_return = std::clamp(number(text.substr(begin, i - begin), "UTAU flag Tnr"), 0., 100.) / 100.;
            continue;
        }
        if (key.size() != 1) continue;
        switch (flag) {
        case 'g': case 't': case 'b': case 'A': case 'O': case 'W': break;
        default: continue;
        }
        const double v = number(text.substr(begin, i - begin), "UTAU flag " + std::string(key));
        switch (flag) {
        case 'g': out.formant_st = -std::clamp(v, -100., 100.) * (12. * std::log2(10.) / 200.); break;
        case 't': out.pitch_st += v / 10.; break;
        case 'b': out.consonant_boost = std::clamp(v, 0., 100.); break;
        case 'A': out.dynamics = std::clamp(v, 0., 100.); break;
        case 'O': out.opening = std::clamp(v, -100., 100.); break;
        case 'W': out.forced_hz = forced_hz(v); break;
        }
    }
}

std::vector<tn::F0Point> read_f0(const fs::path& path) {
    const auto fail = [&](std::size_t line, std::string_view why) {
        bad("F0 CSV '" + path.u8string() + "', line " + std::to_string(line) + ": " + std::string(why));
    };
    if (!fs::is_regular_file(fs::status(path))) fail(1, "not a regular file");
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    const auto size = static_cast<std::streamoff>(in.tellg());
    if (!in || size < 0 || size > static_cast<std::streamoff>(max_f0_bytes)) fail(1, "cannot read or exceeds 64 MiB");
    in.seekg(0);
    std::vector<tn::F0Point> points;
    std::array<char, max_csv_line + 2> buffer{};
    std::size_t line_no = 0;
    bool seen_header = false;
    double previous = -1.;
    for (;;) {
        ++line_no;
        if (!in.getline(buffer.data(), static_cast<std::streamsize>(buffer.size()))) {
            if (in.eof() && !in.bad() && in.gcount() == 0) break;
            fail(line_no, "line too long or unreadable");
        }
        const auto bytes = static_cast<std::size_t>(in.gcount()) - (in.eof() ? 0u : 1u);
        if (bytes > max_csv_line) fail(line_no, "line too long");
        std::string_view text(buffer.data(), bytes);
        if (line_no == 1 && text.substr(0, 3) == "\xef\xbb\xbf") text.remove_prefix(3);
        text = trim(text);
        if (text.empty()) continue;
        const auto comma = text.find(',');
        if (comma == text.npos || text.find(',', comma + 1) != text.npos) fail(line_no, "expected two comma-separated fields");
        const auto time = trim(text.substr(0, comma));
        const auto hz = trim(text.substr(comma + 1));
        if (!seen_header) {
            if (time != "time_s" || hz != "f0_hz") fail(line_no, "expected header time_s,f0_hz");
            seen_header = true;
            continue;
        }
        tn::F0Point p{0., 0.};
        if (!parse_number(time, p.time_s) || !parse_number(hz, p.hz)) fail(line_no, "invalid or nonfinite number");
        if (p.time_s < 0. || p.time_s <= previous) fail(line_no, "times must be nonnegative and strictly increasing");
        if (p.hz != 0. && (p.hz < 30. || p.hz > 2400.)) fail(line_no, "frequencies must be 0 or 30..2400 Hz");
        if (points.size() == max_bends) fail(line_no, "more than 1,000,000 F0 points");
        points.push_back(p);
        previous = p.time_s;
    }
    if (points.empty()) fail(line_no, "no F0 points");
    return points;
}

Job classic(const Args& args) {
    if (args.size() < 11 || args.size() > 13) bad("expected 11 classic arguments, plus optional !tempo and pitchbend");
    Job job;
    job.input = args[0];
    job.output = args[1];
    job.pcm16 = true;
    auto& o = job.options;
    o.loop = true;
    o.peak_normalize = true;
    o.target_hz = note_hz(args[2].u8string());
    o.velocity = number(args[3].u8string(), "velocity");
    flags(args[4].u8string(), o);
    o.offset_ms = number(args[5].u8string(), "offset");
    o.length_ms = number(args[6].u8string(), "length");
    if (!(o.length_ms > 0.)) bad("classic length must be positive");
    o.consonant_ms = number(args[7].u8string(), "consonant");
    o.cutoff_ms = number(args[8].u8string(), "cutoff");
    o.volume = number(args[9].u8string(), "volume");
    o.modulation = number(args[10].u8string(), "modulation");
    if (args.size() >= 12) {
        const auto tempo = args[11].u8string();
        if (tempo.empty() || tempo[0] != '!') bad("classic tempo must be written as !BPM");
        o.tempo = number(std::string_view(tempo).substr(1), "tempo");
        if (!(o.tempo >= 1. && o.tempo <= 1000.)) bad("tempo must be 1..1000 BPM");
    }
    if (args.size() == 13) o.bends = pitch_bends(args[12].u8string());
    return job;
}

Job native(const Args& args) {
    if (args.size() < 3) bad("--render requires input and output paths");
    Job job;
    job.input = args[1];
    job.output = args[2];
    auto& o = job.options;
    struct Numeric { std::string_view flag; double tn::Options::*member; };
    static constexpr Numeric numerics[] = {
        {"--pitch-st", &tn::Options::pitch_st}, {"--target-hz", &tn::Options::target_hz},
        {"--formant-st", &tn::Options::formant_st}, {"--length-ms", &tn::Options::length_ms},
        {"--mask-return", &tn::Options::mask_return}, {"--harmonic-amount", &tn::Options::mask_return},
        {"--offset-ms", &tn::Options::offset_ms}, {"--consonant-ms", &tn::Options::consonant_ms},
        {"--cutoff-ms", &tn::Options::cutoff_ms}, {"--velocity", &tn::Options::velocity},
        {"--volume", &tn::Options::volume},
        {"--consonant-boost", &tn::Options::consonant_boost}, {"--opening", &tn::Options::opening},
        {"--dynamics", &tn::Options::dynamics}, {"--modulation", &tn::Options::modulation},
        {"--tempo", &tn::Options::tempo},
    };
    for (std::size_t i = 3; i < args.size(); ++i) {
        const auto option = args[i].u8string();
        const auto value = [&]() -> const fs::path& {
            if (++i == args.size()) bad("missing value for " + option);
            return args[i];
        };
        if (option == "--f0") {
            if (!job.f0_csv.empty()) bad("--f0 may only be specified once");
            job.f0_csv = value();
            if (job.f0_csv.empty()) bad("--f0 requires a nonempty path");
        } else if (option == "--phase-lock") {
            o.phase_lock = true;
        } else if (option == "--grain-formant") {
            o.grain_formant = true;
        } else if (option == "--loop") {
            o.loop = true;
        } else if (option == "--force-loop") {
            o.loop = o.force_loop = true;
        } else if (option == "--pcm16") {
            job.pcm16 = true;
        } else if (option == "--report") {
            job.report = true;
        } else if (option == "--peak-normalize") {
            o.peak_normalize = true;
        } else if (option == "--forced-hz") {
            o.forced_hz = forced_hz(number(value().u8string(), option));
        } else if (option == "--pitchbend") {
            o.bends = pitch_bends(value().u8string());
        } else {
            const auto found = std::find_if(std::begin(numerics), std::end(numerics),
                                            [&](const Numeric& n) { return n.flag == option; });
            if (found != std::end(numerics)) {
                o.*(found->member) = number(value().u8string(), option);
            } else {
                const auto voice = std::find_if(std::begin(voice_flags), std::end(voice_flags),
                                                [&](const VoiceFlag& f) { return f.native == option; });
                if (voice == std::end(voice_flags)) bad("unknown option " + option);
                o.*(voice->member) = number(value().u8string(), option);
            }
        }
    }
    return job;
}

void print_help(std::ostream& stream) {
    stream << "tn_finds UTAU resampler";
#if TN_FINDS_DUAL_ROLE
    stream << " and phrase wavtool";
#endif
    stream << R"(.

Usage:
  tn_finds input output note velocity flags offset length consonant cutoff volume modulation [!tempo [pitchbend]]
  tn_finds --render input.wav output.wav [options]
  tn_finds --help | --version

)";
#if TN_FINDS_DUAL_ROLE
    stream << R"(Dual-role build. Select it as both tools in the stock Windows classic batch.
Classic resampler calls always write a note model. Clear old PCM caches when
switching builds. Classic wavtool calls collect the note models and render
at the end of the batch. The host concatenates the .whd and .dat files.
Manual phrase rendering uses --prepare-note with the classic resampler
arguments, then --phrase with the wavtool arguments and a final LAST_NOTE.
An explicit --wavtool command uses ordinary PCM overlap-add.

  tn_finds [--phrase | --wavtool] output.wav input.wav offset duration [p1 p2 p3 v1 v2 v3 v4 overlap p4 p5 v5] [LAST_NOTE]

)";
#endif
    stream << R"(Native render options (value follows option as a separate argument):
  --pitch-st N         Relative pitch shift in semitones, default 0.
  --target-hz N        Target pitch in Hz, 0 keeps source contour, default 0.
  --formant-st N       Shift both voiced envelopes by N semitones, default 0.
                       Positive values raise the harmonic and breath envelopes.
  --mask-return N      Breath return from the five-cycle average, 0..1, default .75.
                       Classic Tnr divided by 100. --harmonic-amount is an alias.
  --grain-formant      Use envelope-free grain warping instead of spectral gain.
  --phase-lock         Try source-relative harmonic phase on edited voiced grains.
  --f0 TRACK.csv       Source F0, header time_s,f0_hz; absolute input seconds, 0=no F0 hint.
  --length-ms N        Output duration, -1 keeps the duration of the selected input.
  --offset-ms N        Begin at this point in the input, default 0 ms.
  --cutoff-ms N        End trim (negative: take this many ms after offset), default 0.
  --consonant-ms N     Fixed region size after offset, default 0 ms.
  --velocity N        Consonant velocity, 100 leaves fixed region unchanged.
  --volume N          Output gain in percent, default 100.
  --tension N         Mt: change the LF glottal amplitude by remapping source Rd.
  --breathiness N     Mb: scale recorded breath; positive values fade voiced harmonics.
                      +100 leaves only residual in steady voiced regions; no added noise.
  --mouth N           Mo: spectral gain relative to target H1; positive sounds more open.
  --resonance N       Mr: boost or cut the band around 3 kHz.
  --formant-emphasis N  ME: sharpen or soften existing formants.
                      Mt/Mb/Mo/Mr/ME: -100..100, all default to 0.
  --consonant-boost N  Unvoiced boost, default 0.
  --opening N         Voiced brightness control, default 0.
  --dynamics N        Pitch-linked dynamics, default 0.
  --forced-hz N       Same limits as classic W.
  --modulation N      Retain N percent of the source pitch variation around target.
  --tempo N           BPM of the bend grid, default 120.
  --pitchbend TEXT    Classic signed Base64 bend contour.
  --loop              Allow looping of the stretchable region (default linear).
  --force-loop        Prefer safe sustain loops even for short extensions (Me).
  --peak-normalize    Enable peak normalization (default off).
  --pcm16             Output PCM16 instead of the default IEEE float32.
  --report            Print a JSON timing/render report to stdout.

Classic mode takes a note like C#4 or Bb3 and millisecond timing fields.
If there are no flags, pass an empty argument. Default: loop and normalize.
The standalone resampler writes PCM16. Recognized flags (others ignored):
  g[-100..100]        Shift harmonic and voiced-breath envelopes; positive g lowers.
  tN                  Pitch offset in 10-cent units.
  b[0..100]           Boost unvoiced consonants.
  A[0..100]           Pitch-linked dynamics.
  O[-100..100]        Voiced brightness.
  Tnr[0..100]         Amount of the breath mask returned to the residual. Default Tnr75.
                       Tnr0 uses the original average above the protected low band.
  WN                  Negative=unvoiced; 0 or (0,50)=off; 50..1000=forced Hz.
  e / e1 / Me / Me1   Linear stretching / prefer safe sustain loops; last wins.
                      e0 and Me0 keep the default mapping (OpenUtau option 0).
  Mt/Mb/Mo/Mr/ME      Moresampler-style controls listed above; case-sensitive.
                      Independent approximations, not Moresampler's model.
Pitch bends: signed 12-bit Base64 pairs in cents; #n# repeats the previous
value n additional times. The engine assumes zero after the contour.

WAV: RIFF little-endian PCM 8/16/24/32 or IEEE float 32/64, 8000..192000 Hz.
Stereo is averaged to mono, and output retains the input sample rate.
)";
}

void print_report(const tn::Report& r) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(10) << std::boolalpha
        << "{\"analysis_ms\":" << r.analysis_ms << ",\"synthesis_ms\":" << r.synthesis_ms
        << ",\"source_pitch_hz\":" << r.source_pitch_hz << ",\"input_samples\":" << r.input_samples
        << ",\"output_samples\":" << r.output_samples << ",\"islands\":" << r.islands
        << ",\"bypass\":" << r.bypass << ",\"loop_used\":" << r.loop_used << '}';
    std::cout << out.str() << '\n';
}

int run(const Args& args) {
    if (args.empty()) { print_help(std::cerr); return 2; }
    const auto command = args.front().u8string();
    if (command == "--help" || command == "-h" || command == "-?" ||
        (args.size() == 2 && args[1] == "--help" &&
         (command == "--render" || (TN_FINDS_DUAL_ROLE && (command == "--wavtool" || command == "--phrase"))))) {
        print_help(std::cout);
        return 0;
    }
    if (command == "--version" || command == "-V") {
        std::cout << "tn_finds " << TN_FINDS_VERSION << (TN_FINDS_DUAL_ROLE ? "-phrase\n" : "-resampler\n");
        return 0;
    }
    const bool native_render = command == "--render";
#if TN_FINDS_DUAL_ROLE
    using tn::wavtool::JoinMode;
    if (command == "--phrase" || command == "--wavtool") {
        tn::wavtool::append(args, 1, command == "--phrase" ? JoinMode::explicit_phrase : JoinMode::pcm);
        return 0;
    }
    const bool prepare = command == "--prepare-note";
    if (!native_render && !prepare && tn::wavtool::is_classic_call(args)) {
        tn::wavtool::append(args, 0, JoinMode::classic_phrase);
        return 0;
    }
    auto job = native_render ? native(args) : classic(prepare ? Args(args.begin()+1, args.end()) : args);
#else
    if (command == "--phrase" || command == "--wavtool" || command == "--prepare-note")
        bad("select the phrase build for wavtool or note-model commands");
    auto job = native_render ? native(args) : classic(args);
#endif
    tn::wav::require_distinct(job.input, job.output);
    if (!job.f0_csv.empty()) {
        tn::wav::require_distinct(job.f0_csv, job.output);
        job.options.f0 = read_f0(job.f0_csv);
    }
    auto audio = tn::wav::read(job.input);
#if TN_FINDS_DUAL_ROLE
    if (!native_render) {
        const auto model = tn::prepare_note(audio.samples, audio.sample_rate, job.options);
        tn::wav::require_distinct(job.input, job.output);
        tn::wavtool::write_model(job.output, model);
        return 0;
    }
#endif
    tn::Report report;
    const auto rendered = tn::render(audio.samples, audio.sample_rate, job.options, job.report ? &report : nullptr);
    // Recheck after computation, immediately before opening the output.
    tn::wav::require_distinct(job.input, job.output);
    if (!job.f0_csv.empty()) tn::wav::require_distinct(job.f0_csv, job.output);
    tn::wav::write(job.output, rendered, audio.sample_rate, job.pcm16);
    if (job.report) print_report(report);
    return 0;
}

#ifdef _WIN32
using ArgChar = wchar_t;
#else
using ArgChar = char;
#endif
} // namespace

#ifdef _WIN32
int wmain(int argc, wchar_t** argv) { // MinGW: link the CLI target with -municode.
#else
int main(int argc, char** argv) {
#endif
    try {
        Args args;
        if (argc > 1) args.reserve(static_cast<std::size_t>(argc - 1));
        for (int i = 1; i < argc; ++i) args.emplace_back(static_cast<const ArgChar*>(argv[i]));
        return run(args);
    } catch (const std::exception& e) {
        std::cerr << "tn_finds: " << e.what() << '\n';
        return 1;
    }
}
