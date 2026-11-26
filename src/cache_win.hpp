// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <bcrypt.h>

#include <array>
#include <atomic>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace tn::detail::cache {
namespace fs = std::filesystem;
using Digest = std::array<unsigned char, 32>;
struct Bytes { const void* data; std::size_t size; };
template<class T> Bytes bytes(const T& value) { return {&value, sizeof(value)}; }
template<class T> Bytes bytes(const std::vector<T>& values) { return {values.data(), values.size()*sizeof(T)}; }

class Sha256 {
    BCRYPT_ALG_HANDLE provider_ = nullptr;
public:
    Sha256() { BCryptOpenAlgorithmProvider(&provider_, BCRYPT_SHA256_ALGORITHM, nullptr, 0); }
    ~Sha256() { if (provider_) BCryptCloseAlgorithmProvider(provider_, 0); }
    Sha256(const Sha256&) = delete;
    Sha256& operator=(const Sha256&) = delete;

    bool hash(Digest& out, std::initializer_list<Bytes> parts) const {
        if (!provider_) return false;
        BCRYPT_HASH_HANDLE handle = nullptr;
        if (BCryptCreateHash(provider_, &handle, nullptr, 0, nullptr, 0, 0) < 0) return false;
        bool ok = true;
        for (const auto& part : parts) {
            if (!part.size) continue;
            if (part.size > std::numeric_limits<ULONG>::max() ||
                BCryptHashData(handle, static_cast<PUCHAR>(const_cast<void*>(part.data)),
                               static_cast<ULONG>(part.size), 0) < 0) { ok = false; break; }
        }
        ok = ok && BCryptFinishHash(handle, out.data(), static_cast<ULONG>(out.size()), 0) >= 0;
        BCryptDestroyHash(handle);
        return ok;
    }
};

inline std::wstring hex(const Digest& key) {
    constexpr wchar_t digits[] = L"0123456789abcdef";
    std::wstring name;
    name.reserve(2*key.size());
    for (const auto byte : key) {
        name += digits[byte >> 4];
        name += digits[byte & 15];
    }
    return name;
}

inline std::wstring environment(const wchar_t* variable) {
    wchar_t local[512];
    const DWORD n = GetEnvironmentVariableW(variable, local, 512);
    if (!n) return {};
    if (n < 512) return std::wstring(local, n);
    std::wstring value(n, L'\0');
    const DWORD length = GetEnvironmentVariableW(variable, value.data(), n);
    if (!length || length >= n) return {};
    value.resize(length);
    return value;
}

inline fs::path directory(const wchar_t* variable, const wchar_t* kind) {
    auto root = environment(variable);
    if (root == L"off") return {};
    if (root.empty()) root = environment(L"LOCALAPPDATA");
    const fs::path path(root);
    return path.is_absolute() ? path / L"tn_finds" / kind : fs::path{};
}

// Local Windows processes sharing a key need only one cold computation.
// The kernel releases a mutex if the process dies before publishing.
class EntryLock {
    HANDLE handle_;
public:
    bool owned = false;
    explicit EntryLock(const std::wstring& name) : handle_(CreateMutexW(nullptr, FALSE, name.c_str())) {
        if (handle_) {
            const auto status = WaitForSingleObject(handle_, 30'000);
            owned = status == WAIT_OBJECT_0 || status == WAIT_ABANDONED;
        }
    }
    ~EntryLock() {
        if (owned) ReleaseMutex(handle_);
        if (handle_) CloseHandle(handle_);
    }
    EntryLock(const EntryLock&) = delete;
    EntryLock& operator=(const EntryLock&) = delete;
};

inline void publish(const fs::path& path, std::initializer_list<Bytes> parts) {
    static std::atomic<unsigned int> serial{0};
    auto tmp = path;
    tmp += L"." + std::to_wstring(GetCurrentProcessId()) + L"." + std::to_wstring(serial++) + L".tmp";
    std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
    if (!out) return;
    for (const auto& part : parts)
        if (part.size) out.write(static_cast<const char*>(part.data), static_cast<std::streamsize>(part.size));
    out.close();
    if (!out || !MoveFileExW(tmp.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING)) {
        std::error_code ignored;
        fs::remove(tmp, ignored);
    }
}
} // namespace tn::detail::cache
