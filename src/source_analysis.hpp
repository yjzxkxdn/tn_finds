// SPDX-License-Identifier: GPL-3.0-only
#pragma once

#include "glottal.hpp"
#include "harmonic_mask.hpp"

namespace tn::detail {
// Keep the target pitch, time map and phase state out of the persistent source data.
// The same immutable analysis feeds standalone WAVs and dual-role note models.
struct SourceAnalysis {
    pitch::Track clock;
    pitch::Island pulses;
    HarmonicSplit split;
    GlottalTrack rd;
};

struct SourceAnalysisInput {
    const std::vector<float>& recording;
    const std::vector<float>& segment;
    const std::vector<F0Point>& f0;
    const pitch::Track& mrq_f0;
    int rate, first; // segment begins here on the untrimmed source clock
    double forced_hz, mask_return;
    bool cleanup, need_rd;
};

#ifdef TN_FINDS_ANALYSIS_CACHE_ID
SourceAnalysis cached_source_analysis(const SourceAnalysisInput&, SourceAnalysis (*)(const SourceAnalysisInput&),
                                      bool& hit);
#endif
} // namespace tn::detail
