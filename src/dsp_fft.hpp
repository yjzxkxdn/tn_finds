// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <cmath>
#include <complex>
#include <cstdint>
#include <utility>
#include <vector>

#if defined(__SSE2__) || defined(_M_X64)
#include <emmintrin.h>
#endif

namespace tn::detail {
using Complex = std::complex<float>;

// All callers provide finite, bounded DSP samples and a power-of-two length.
// Transform is not normalized.
class Fft {
    int n_;
    std::vector<std::pair<int, int>> swaps_;
    std::vector<Complex> roots_;

    template<bool inverse>
    void transform(std::vector<Complex>& a) const {
        for (const auto& pair : swaps_) std::swap(a[pair.first], a[pair.second]);
        for (int half = 1; half < n_; half *= 2) {
            const auto* roots = roots_.data() + half - 1;
            for (int i = 0; i < n_; i += 2*half) {
                auto* left = a.data() + i;
                auto* right = left + half;
                int j = 0;
#if defined(__SSE2__) || defined(_M_X64)
                // Two independent complex butterflies on baseline x64/SSE2.
                const auto signs = _mm_castsi128_ps(_mm_set_epi32(0, INT32_MIN, 0, INT32_MIN));
                for (; j+1 < half; j += 2) {
                    const auto w = _mm_loadu_ps(reinterpret_cast<const float*>(roots+j));
                    auto* l = reinterpret_cast<float*>(left+j);
                    auto* r = reinterpret_cast<float*>(right+j);
                    const auto u = _mm_loadu_ps(l), b = _mm_loadu_ps(r);
                    const auto wr = _mm_shuffle_ps(w, w, _MM_SHUFFLE(2,2,0,0));
                    const auto wi = _mm_shuffle_ps(w, w, _MM_SHUFFLE(3,3,1,1));
                    const auto swapped = _mm_shuffle_ps(b, b, _MM_SHUFFLE(2,3,0,1));
                    const auto imag = _mm_mul_ps(swapped, wi);
                    auto term = _mm_xor_ps(imag, signs);
                    if constexpr (inverse) term = _mm_xor_ps(term, _mm_set1_ps(-0.f));
                    const auto v = _mm_add_ps(_mm_mul_ps(b, wr), term);
                    _mm_storeu_ps(l, _mm_add_ps(u, v));
                    _mm_storeu_ps(r, _mm_sub_ps(u, v));
                }
#endif
                for (; j < half; ++j) {
                    const float wr = roots[j].real();
                    const float wi = inverse ? -roots[j].imag() : roots[j].imag();
                    const Complex u = left[j], b = right[j];
                    // Use the finite complex product directly. This keeps the
                    // butterfly arithmetic and avoids libgcc's NaN recovery.
                    const Complex v(b.real()*wr-b.imag()*wi, b.real()*wi+b.imag()*wr);
                    left[j] = u + v;
                    right[j] = u - v;
                }
            }
        }
    }

public:
    explicit Fft(int n) : n_(n), roots_(n-1) {
        std::vector<int> reverse(n);
        swaps_.reserve(n/2);
        for (int i = 1; i < n; ++i) {
            reverse[i] = (reverse[i >> 1] >> 1) | ((i & 1) * (n >> 1));
            if (i < reverse[i]) swaps_.emplace_back(i, reverse[i]);
        }
        if (n < 2) return;
        auto* full = roots_.data() + n/2 - 1;
        for (int i = 0; i < n/2; ++i) {
            const double angle = -6.283185307179586476925286766559 * i / n;
            full[i] = Complex(static_cast<float>(std::cos(angle)), static_cast<float>(std::sin(angle)));
        }
        // Lay out the twiddles in butterfly order so both operands are
        // contiguous in every stage.
        for (int half = 1; half < n/2; half *= 2) {
            const int step = n/(2*half);
            for (int j = 0; j < half; ++j) roots_[half-1+j] = full[j*step];
        }
    }

    void run(std::vector<Complex>& a, bool inverse = false) const {
        if (inverse) transform<true>(a);
        else transform<false>(a);
    }
};

// Invert a real signal's nonnegative spectrum with one half-size FFT.
// As with Fft, gain=1 gives the unnormalized inverse (N times the signal).
class RealInverseFft {
    Fft fft_;
    std::vector<Complex> packed_, rotation_;

public:
    explicit RealInverseFft(int n) : fft_(n/2), packed_(n/2), rotation_(n/2) {
        for (int k = 0; k < n/2; ++k) {
            const double angle = -6.283185307179586476925286766559*k/n;
            rotation_[k] = {static_cast<float>(std::sin(angle)), static_cast<float>(std::cos(angle))};
        }
    }

    void run(const std::vector<Complex>& spectrum, std::vector<float>& audio, float gain = 1.f) {
        const auto half = packed_.size();
        for (std::size_t k = 0; k < half; ++k) {
            // Sum and difference give the even and odd time-index spectra.
            // Pack them as the real and imaginary parts of one inverse FFT.
            const auto a = spectrum[k], b = std::conj(spectrum[half-k]);
            const auto sum = a+b, diff = a-b, w = rotation_[k];
            packed_[k] = {sum.real()+(w.real()*diff.real()-w.imag()*diff.imag()),
                          sum.imag()+(w.real()*diff.imag()+w.imag()*diff.real())};
        }
        fft_.run(packed_, true);
        audio.resize(2*half);
        for (std::size_t k = 0; k < half; ++k) {
            audio[2*k] = gain*packed_[k].real();
            audio[2*k+1] = gain*packed_[k].imag();
        }
    }
};
} // namespace tn::detail
