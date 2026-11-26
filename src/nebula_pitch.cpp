// SPDX-License-Identifier: GPL-3.0-only
#include "pitch.hpp"

extern "C" {
#include <nebula.h>
}

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <stdexcept>

namespace tn::pitch {
namespace {
int native_hop(int rate) {
    int step = 1;
    while (step*200 < rate) step *= 2;
    return step;
}

Track analyze_nebula(const std::vector<float>& input, int rate) {
    Track result;
    // Nebula uses a 5 ms grid on the source clock.
    result.step = native_hop(rate);
    const auto [low, high] = std::minmax_element(input.begin(), input.end());
    if (!std::isfinite(*low) || !std::isfinite(*high) || *low < -64.f || *high > 64.f)
        throw std::invalid_argument("Nebula expects finite audio with magnitude <= 64");
    // Keep two 30 ms frames before the voicing and refinement passes.
    const auto minimum = static_cast<std::size_t>((65*rate+999)/1000);
    if (input.size() < minimum || *high-*low < 1.e-10f) {
        result.frames.resize(std::max<std::size_t>(1, input.size()/result.step));
        return result;
    }

    auto waveform = input;
    int frames = 0;
    std::unique_ptr<float, decltype(&std::free)> f0(nullptr, &std::free);
    {
        // The upstream dithering implementation seeds the C runtime RNG.
        static std::mutex mutex;
        const std::lock_guard<std::mutex> lock(mutex);
        f0.reset(f0_from_wave(70.f, 800.f, static_cast<float>(result.step)/rate,
                             1, 0, .1f, .2f, 2.f, 1.5f, 0.f, 0.f, 3, 10, 10,
                             waveform.data(), static_cast<int>(waveform.size()), rate, &frames));
    }
    if (!f0 || frames <= 0 || static_cast<std::size_t>(frames) > 1+input.size()/result.step)
        throw std::runtime_error("Nebula returned an invalid F0 track");
    result.frames.reserve(frames);
    for (int i = 0; i < frames; ++i) {
        const double hz = f0.get()[i];
        if (std::isfinite(hz) && hz >= 30. && hz <= 2400.) result.frames.push_back({hz, 1.});
        else result.frames.push_back({});
    }
    return result;
}
} // namespace

#ifdef TN_FINDS_NEBULA_CACHE_ID
Track cached_nebula(const std::vector<float>&, int rate, int step,
                    Track (*analyze)(const std::vector<float>&, int));
#endif

Track estimate_nebula(const std::vector<float>& input, int rate) {
    if (input.empty() || input.size() > 16'000'000 || rate < 8000 || rate > 192000) return {};
#ifdef TN_FINDS_NEBULA_CACHE_ID
    return cached_nebula(input, rate, native_hop(rate), analyze_nebula);
#else
    return analyze_nebula(input, rate);
#endif
}

} // namespace tn::pitch
