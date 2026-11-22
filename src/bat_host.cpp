// SPDX-License-Identifier: GPL-3.0-only
#include "bat_host.hpp"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#endif

namespace tn::wavtool {
#ifdef _WIN32
namespace {
namespace fs = std::filesystem;

fs::path environment_path(const wchar_t* name) {
    const DWORD count = GetEnvironmentVariableW(name, nullptr, 0);
    if (count <= 1 || count > 32768) return {};
    std::wstring text(count, L'\0');
    const DWORD length = GetEnvironmentVariableW(name, text.data(), count);
    if (length == 0 || length >= count) return {};
    text.resize(length);
    if (text.size() >= 2 && text.front() == L'"' && text.back() == L'"')
        text = text.substr(1, text.size()-2);
    return text;
}

bool same_path(const fs::path& a, const fs::path& b) {
    if (a.empty() || b.empty()) return false;
    std::error_code ec1, ec2;
    const auto x = fs::weakly_canonical(a, ec1), y = fs::weakly_canonical(b, ec2);
    return !ec1 && !ec2 && CompareStringOrdinal(x.c_str(), -1, y.c_str(), -1, TRUE) == CSTR_EQUAL;
}

bool same_exe(const fs::path& tool) {
    if (tool.empty()) return false;
    std::wstring buffer(32768, L'\0');
    const auto n = GetModuleFileNameW(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
    if (n == 0 || n >= buffer.size()) return false;
    buffer.resize(n);
    const fs::path self(buffer);
    if (same_path(self, tool)) return true;
    // OpenUtau keeps separate copies in its two tool directories.
    std::ifstream a(self, std::ios::binary | std::ios::ate), b(tool, std::ios::binary | std::ios::ate);
    const auto size = static_cast<std::streamoff>(a.tellg());
    if (!a || !b || size <= 0 || size != b.tellg() || size > 32*1024*1024) return false;
    a.seekg(0);
    b.seekg(0);
    std::array<char, 65536> x{}, y{};
    for (auto rest = size; rest > 0;) {
        const auto length = std::min(rest, static_cast<std::streamoff>(x.size()));
        if (!a.read(x.data(), length) || !b.read(y.data(), length) ||
            !std::equal(x.begin(), x.begin()+length, y.begin())) return false;
        rest -= length;
    }
    return true;
}

bool starts(std::string_view line, std::string_view prefix) { return line.substr(0,prefix.size()) == prefix; }
bool blank(char c) { return c == ' ' || c == '\t' || c == '\r'; }

std::string command(std::string_view line) {
    std::string out;
    for (char c : line) {
        if (blank(c) || (c == '@' && out.empty())) continue;
        out.push_back(c >= 'A' && c <= 'Z' ? char(c-'A'+'a') : c);
    }
    return out;
}

std::vector<std::string> commands(const fs::path& path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    const auto bytes = static_cast<std::streamoff>(in.tellg());
    if (!in || bytes <= 0 || bytes > 4*1024*1024) return {};
    in.seekg(0);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        if (lines.empty() && starts(line, "\xef\xbb\xbf")) line.erase(0,3);
        if (line.size() > 8192 || line.find('\0') != line.npos) return {};
        // Refuse compound command lines. Opaque quoted paths and flags can
        // contain metacharacters. Standard setup may redirect errors to nul.
        bool quoted = false;
        for (char c : line) {
            if (c == '"') quoted = !quoted;
            if (!quoted && (c == '&' || c == '|' || c == '^')) return {};
        }
        auto s = command(line);
        if (!s.empty() && !starts(s,"rem")) lines.push_back(std::move(s));
    }
    return lines;
}

bool standard_helper() {
    const auto helper = environment_path(L"helper");
    if (helper.empty() || helper.filename() != L"temp_helper.bat") return false;
    const auto lines = commands(helper);
    int resamp = 0, wave = 0;
    for (const auto& line : lines) {
        if (line == "\"%tool%\"\"%output%\"%temp%%stp%%3%env%") ++wave;
        else if (line == "\"%resamp%\"%1%temp%%2%vel%%flag%%5%6%7%8%params%") ++resamp;
        else if (line == ":a" || line == "ifexist%temp%gotoa" ||
                 line == "ifexist\"%cachedir%\\%9_*.wav\"del\"%cachedir%\\%9_*.wav\"" ||
                 line == "settemp=\"temp$$$.wav\"") continue;
        else return false;
    }
    return resamp == 1 && wave == 1 && lines.back() == "\"%tool%\"\"%output%\"%temp%%stp%%3%env%";
}
} // namespace

int classic_calls(const fs::path& output) {
    if (!same_path(output, environment_path(L"output")) || environment_path(L"tool").empty()) return 0;
    const auto lines = commands("temp.bat");
    constexpr std::string_view footer[] = {
        "ifnotexist\"%output%.whd\"gotoe", "ifnotexist\"%output%.dat\"gotoe",
        "copy/y\"%output%.whd\"/b+\"%output%.dat\"/b\"%output%\"",
        "del\"%output%.whd\"", "del\"%output%.dat\"", ":e"
    };
    constexpr auto n = std::size(footer);
    if (lines.size() <= n || !std::equal(std::begin(footer), std::end(footer), lines.end()-n)) return 0;
    int calls = 0;
    bool resets = false, helper = false;
    for (std::size_t i = 0; i+n < lines.size(); ++i) {
        const auto& s = lines[i];
        if ((s == "del\"%output%\"2>nul" || s == "del\"%output%\"") && !resets && calls == 0) resets = true;
        else if (starts(s,"call%helper%") || starts(s,"call\"%helper%\"")) { ++calls; helper = true; }
        else if (starts(s,"\"%tool%\"\"%output%\"")) ++calls;
        else if (starts(s,"set") || starts(s,"echo") || starts(s,"mkdir\"%cachedir%\"")) continue;
        else return 0;
    }
    if (!resets || calls > 20000 || (helper && !standard_helper())) return 0;
    return calls;
}

bool is_dual_tool_cache(const fs::path& note) {
    return same_path(note, environment_path(L"temp")) && same_exe(environment_path(L"resamp"));
}
#else
int classic_calls(const std::filesystem::path&) { return 0; }
bool is_dual_tool_cache(const std::filesystem::path&) { return false; }
#endif
} // namespace tn::wavtool
