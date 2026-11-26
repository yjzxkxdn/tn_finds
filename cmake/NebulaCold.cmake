# SPDX-License-Identifier: GPL-3.0-only
# The documented OpenMoresampler revision computes a 30 ms coarse track only
# to log its median. tn_finds discards the log, but the coarse dither advances
# the CRT PRNG before the real (fine grid) analysis. Keep the same preprocessing
# and exactly two rand() calls per dithered sample; omit the unused coarse DSP.
# Use the unmodified upstream code for other revisions. No source tree or model
# is modified or copied into Git by this build-time adaptation.
option(TN_FINDS_NEBULA_FAST_COLD "Omit the log-only coarse pass for the pinned Nebula sources" ON)
if(NOT TN_FINDS_NEBULA_FAST_COLD)
    return()
endif()

if(NOT nebula_pinned)
    message(STATUS "Using unmodified Nebula analysis (different upstream source revision)")
    return()
endif()

file(READ "${nebula_root}/nebula/nebula_wave.c" wave)
file(READ "${nebula_root}/nebula/nebula_f0.c" estimate)
string(REPLACE "\r\n" "\n" wave "${wave}")
string(REPLACE "\r\n" "\n" estimate "${estimate}")
string(FIND "${estimate}" "    int nx = nx_in;" prefix_begin)
string(FIND "${estimate}" "    int *naxis = (int *)calloc" prefix_end)
math(EXPR prefix_size "${prefix_end} - ${prefix_begin}")
string(SUBSTRING "${estimate}" ${prefix_begin} ${prefix_size} prefix)
string(REPLACE "    int Count = (int)(nx * 0.00012500001f / thop);\n" "" prefix "${prefix}")
string(REPLACE "NEB_FS / fs" "8000.0f / fs" prefix "${prefix}")
string(REPLACE "                x[i] = nebula_normal() * (th - rms[i]) + x[i];"
               "                (void)rand();\n                (void)rand();" prefix "${prefix}")
string(REPLACE "// two-pass F0: coarse nebula_est, then refine"
               "static void tn_advance_coarse_dither(float *wave, int nx_in, float fs, float fmin, float dither_level)\n{\n${prefix}    free(x);\n}\n"
               wave "${wave}")
string(FIND "${wave}" "    nebula_result *r0 = nebula_est(" coarse_begin)
string(FIND "${wave}" "    float *xd =" coarse_end)
math(EXPR coarse_size "${coarse_end} - ${coarse_begin}")
string(SUBSTRING "${wave}" ${coarse_begin} ${coarse_size} coarse)
string(REPLACE "${coarse}"
               "    (void)use_peak_likelihood;\n    (void)lmap_smooth;\n    tn_advance_coarse_dither(x, nx, fs_f, fmin, dither_level);\n\n"
               nebula_wave_override "${wave}")
message(STATUS "Nebula: skipping log-only coarse DSP, retaining its dither sequence")
