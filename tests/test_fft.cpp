// SPDX-License-Identifier: GPL-3.0-only
#include "dsp_fft.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    using tn::detail::Complex;
    for (int n : {2, 32, 256, 2048, 16384}) {
        const int half = n/2;
        std::vector<Complex> spectrum(n);
        spectrum[0] = -.125f;
        spectrum[half] = .25f;
        for (int i = 1; i < half; ++i) {
            spectrum[i] = {static_cast<float>(std::cos(i*1.3)), static_cast<float>(std::sin(i*.59))};
            spectrum[n-i] = std::conj(spectrum[i]);
        }
        std::vector<float> wave;
        tn::detail::RealInverseFft(n).run(spectrum, wave, 1.f/n);
        tn::detail::Fft(n).run(spectrum, true);
        for (int i = 0; i < n; ++i) {
            if (!(std::abs(spectrum[i].real()/n-wave[i]) < 2.e-6f)) {
                std::fprintf(stderr,"Inverse FFT error at N=%d sample %d\n", n, i);
                return 1;
            }
        }
    }
}
