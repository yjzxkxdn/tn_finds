// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include <vector>

namespace tn::pitch {

// Zero Hz is unvoiced. Confidence is a periodicity/track-fit score in [0, 1],
// not a probability or a gate for a Track supplied by the caller.
struct Frame {
    double hz = 0, confidence = 0;
};

// Centers are at source samples i * step (not window starts).
struct Track {
    int step = 1;
    std::vector<Frame> frames;
};

struct Island {
    double begin = 0, end = 0;  // sample coordinates; end exclusive, <= x.size()
    std::vector<double> marks;  // same-phase fractional sample positions, not necessarily GCIs
};

// Mono audio, 8--192 kHz. The analysis searches 55--1300 Hz, in frames about
// 5 ms apart. Invalid rates and empty input return an empty track. An input
// shorter than the supporting span of a low note need not yield voiced F0.
Track estimate(const std::vector<float>& x, int sr);

// Build a continuous trial clock, not a voicing decision. Fill missing or
// weak proposals using local low-band periodicity and neighboring anchors.
// Where the waveform has no credible F0, continue a neutral trial clock.
// Positive supplied F0 is never replaced. Filtering decisions belong to the
// harmonic separator, not to this clock's Hz or confidence.
Track candidate_clock(const std::vector<float>& x, int sr, const Track& anchors, bool supplied);

// Works with estimated or supplied F0 (e.g. CSV). Accepts 30--2400 Hz, further
// limited to sr/3. Gaps and out-of-range/nonfinite Hz are unvoiced. Estimates
// the phase origin from low-pass audio and refines pulses locally within
// conservative bounds; an incompatible forced track is not a guarantee of
// the acoustic F0. Islands with fewer than two marks are omitted.
std::vector<Island> locate(const std::vector<float>& x, int sr, const Track& tr);

// Linear Hz/score within voiced frames (also linear score between unvoiced
// frames). Switch at the midpoint at a voiced/unvoiced boundary, never
// interpolate frequency toward 0. Outside the track clamp to the endpoint.
Frame sample(const Track& tr, double source_sample);

}  // namespace tn::pitch
