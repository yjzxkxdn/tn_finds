# SPDX-License-Identifier: GPL-3.0-only
# Adapt private hot loops of the verified upstream integration without
# copying model tables or making local changes to the dependency checkout.
option(TN_FINDS_NEBULA_FAST_KERNELS "Use workspace-reusing Nebula kernels for the pinned sources" ON)
if(NOT TN_FINDS_NEBULA_FAST_KERNELS OR NOT nebula_pinned)
    return()
endif()

function(tn_replace variable old new)
    string(FIND "${${variable}}" "${old}" at)
    if(at LESS 0)
        message(FATAL_ERROR "Unrecognized layout in the pinned Nebula kernel: ${old}")
    endif()
    string(REPLACE "${old}" "${new}" result "${${variable}}")
    set(${variable} "${result}" PARENT_SCOPE)
endfunction()

file(READ "${nebula_root}/nebula/nebula_features.c" features)
string(REPLACE "\r\n" "\n" features "${features}")
string(FIND "${features}" "// complex overlap-add convolution;" conv_begin)
string(FIND "${features}" "// after first decimation:" conv_end)
math(EXPR conv_size "${conv_end} - ${conv_begin}")
string(SUBSTRING "${features}" ${conv_begin} ${conv_size} convolution)
tn_replace(features "${convolution}" "#include \"nebula_convolution.h\"\n\n")
tn_replace(features "const float *h, const float *hs, int nh)" "tn_nebula_conv *h, tn_nebula_conv *hs, int nh)")
tn_replace(features "fftconv_complex(yr, yi, nwork, h, nh, y2r, y2i);" "tn_nebula_conv_run(h, yr, yi, nwork, y2r, y2i);")
tn_replace(features "fftconv_complex(pw, NULL, nwork, hs, nh, y2r, NULL);" "tn_nebula_conv_run(hs, pw, NULL, nwork, y2r, NULL);")
tn_replace(features "        float *mix_r =" [=[        tn_nebula_conv w_conv = tn_nebula_conv_make(w, nw);
        tn_nebula_conv wf_conv = tn_nebula_conv_make(wf, nh);
        tn_nebula_conv ws_conv = tn_nebula_conv_make(ws, nh);
        float *mix_r =]=])
tn_replace(features "fftconv_complex(mix_r, mix_i, nx, w, nw, conv_r, conv_i);" "tn_nebula_conv_run(&w_conv, mix_r, mix_i, nx, conv_r, conv_i);")
tn_replace(features "x2as_decimated(dec_r, dec_i, nwork, wf, ws, nh);" "x2as_decimated(dec_r, dec_i, nwork, &wf_conv, &ws_conv, nh);")
tn_replace(features "        free(w);" [=[        free(w_conv.kernel);
        free(wf_conv.kernel);
        free(ws_conv.kernel);
        free(w);]=])
# Parallel units write disjoint bands. Balance uneven filter widths and
# heterogeneous CPU cores; no thread reduction changes the arithmetic.
tn_replace(features "#pragma omp parallel for" "#pragma omp parallel for schedule(dynamic, 1)")
tn_replace(features "ifdetector *ifd2 = cig_create_ifdetector(fc[b] + fc[b], fres[b]);"
                    "ifdetector *ifd2 = IF2 ? cig_create_ifdetector(fc[b] + fc[b], fres[b]) : NULL;")
# The upstream reduction can round differently on an unaligned input view.
# Keep the same malloc alignment, but reuse one copied window per band.
tn_replace(features [=[        for (t = 0; t < nfrm; ++t) {
            int nhif = ifd->nh;
            float *frame = fetch_frame(x, nx, naxis[t], nhif);]=]
[=[        const int nhif = ifd->nh;
        float *frame = (float*)malloc((size_t)nhif*sizeof(float));
        for (t = 0; t < nfrm; ++t) {
            const int first = naxis[t] - nhif/2;
            if (first >= 0 && first <= nx-nhif) {
                memcpy(frame, x+first, (size_t)nhif*sizeof(float));
            } else {
                for (int i = 0; i < nhif; ++i) {
                    const int j = first+i;
                    frame[i] = j >= 0 && j < nx ? x[j] : 0.f;
                }
            }]=])
tn_replace(features [=[            free(frame);
        }
        cig_delete_ifdetector(ifd);]=] [=[        }
        free(frame);
        cig_delete_ifdetector(ifd);]=])
# The heterodyne oscillators are only consumed by the mix that follows them.
# Keep the upstream interleaved recurrence and its exact float rounding, but
# write the products straight into the mix buffers instead of materializing a
# 2*nx table.
tn_replace(features "// per-band SNR/IF worker: heterodyne, decimate, x2as, interpolate, scatter."
[=[// Same interleaved cosine/sine recurrence as make_oscillator, fused with the
// heterodyne products so no 2*nx oscillator table is ever materialized.
static void tn_heterodyne(const float *x, int nx, float carrier, float *re, float *im)
{
    const double omega = (double)carrier;
    const float c = (float)cos(omega);
    const float sn = (float)sin(omega);
    const float c2 = c + c;
    float prev2r = 1.0f, prev1r = c, prev2i = 0.0f, prev1i = sn;
    int i;
    if (nx > 0) {
        re[0] = 1.0f * x[0];
        im[0] = x[0] * 0.0f;
    }
    if (nx > 1) {
        re[1] = c * x[1];
        im[1] = x[1] * sn;
    }
    for (i = 2; i < nx; ++i) {
        const float r = prev1r * c2 - prev2r;
        const float q = prev1i * c2 - prev2i;
        re[i] = r * x[i];
        im[i] = x[i] * q;
        prev2r = prev1r;
        prev1r = r;
        prev2i = prev1i;
        prev1i = q;
    }
}

// per-band SNR/IF worker: heterodyne, decimate, x2as, interpolate, scatter.]=])
tn_replace(features [=[        float *osc = make_oscillator(nx, (float)carrier);
        int i, t;

        for (i = 0; i < nx; ++i) {
            mix_r[i] = osc[2 * i] * x[i];
            mix_i[i] = x[i] * osc[2 * i + 1];
        }]=]
[=[        int i, t;

        tn_heterodyne(x, nx, carrier, mix_r, mix_i);]=])
tn_replace(features "        free(osc);\n" "")
tn_replace(features [=[            float *osc2 = make_oscillator(nx, (float)(carrier + carrier));
            for (i = 0; i < nx; ++i) {
                mix_r[i] = osc2[2 * i] * x[i];
                mix_i[i] = x[i] * osc2[2 * i + 1];
            }]=]
[=[            tn_heterodyne(x, nx, (float)(carrier + carrier), mix_r, mix_i);]=])
tn_replace(features "            free(osc2);\n" "")
# Only the centred nx-sample window of the convolution is read afterwards. Keep
# that window in place of the padded result so the per-band PAD_mix/PAD_conv
# allocations and the two nx-sample copies disappear; the block grid and the
# per-bin addition order are untouched.
tn_replace(features [=[        conv_r = (float *)calloc((size_t)(nw + nx), sizeof(*conv_r));
        conv_i = (float *)calloc((size_t)(nw + nx), sizeof(*conv_i));
        tn_nebula_conv_run(&w_conv, mix_r, mix_i, nx, conv_r, conv_i);
        if (nx > 0) {
            memcpy(mix_r, conv_r + nw / 2, (size_t)(4 * nx));
            memcpy(mix_i, conv_i + nw / 2, (size_t)(4 * nx));
        }
        free(conv_r);
        free(conv_i);]=]
[=[        conv_r = (float *)calloc((size_t)(nx > 0 ? nx : 1), sizeof(*conv_r));
        conv_i = (float *)calloc((size_t)(nx > 0 ? nx : 1), sizeof(*conv_i));
        tn_nebula_conv_run_shift(&w_conv, mix_r, mix_i, nx, nw / 2, conv_r, conv_i);
        free(mix_r);
        free(mix_i);
        mix_r = conv_r;
        mix_i = conv_i;]=])
tn_replace(features [=[            conv_r = (float *)calloc((size_t)(nw + nx), sizeof(*conv_r));
            conv_i = (float *)calloc((size_t)(nw + nx), sizeof(*conv_i));
            tn_nebula_conv_run(&w_conv, mix_r, mix_i, nx, conv_r, conv_i);
            if (nx > 0) {
                memcpy(mix_r, conv_r + nw / 2, (size_t)(4 * nx));
                memcpy(mix_i, conv_i + nw / 2, (size_t)(4 * nx));
            }
            free(conv_r);
            free(conv_i);]=]
[=[            conv_r = (float *)calloc((size_t)(nx > 0 ? nx : 1), sizeof(*conv_r));
            conv_i = (float *)calloc((size_t)(nx > 0 ? nx : 1), sizeof(*conv_i));
            tn_nebula_conv_run_shift(&w_conv, mix_r, mix_i, nx, nw / 2, conv_r, conv_i);
            free(mix_r);
            free(mix_i);
            mix_r = conv_r;
            mix_i = conv_i;]=])
# Every band allocates its own buffers, hands them to one of the parallel
# workers and frees them again. Keep one block per (thread, slot) instead.
# Values are unchanged because each producer writes the bytes its consumer
# reads; only the decimated buffers rely on a zeroed tail, so those two slots
# are still cleared.
tn_replace(features "// Same interleaved cosine/sine recurrence as make_oscillator, fused with the"
[=[// One reusable scratch block per thread and slot. Contents match a fresh
// calloc: every path writes the bytes it later reads, except the decimation
// buffers whose unused tail must stay zero (slot zero flag).
static void *tn_band_scratch(int slot, size_t bytes, int zero) {
    static __thread void *block[12];
    static __thread size_t capacity[12];
    if (bytes < 1) bytes = 1;
    if (block[slot] == NULL || bytes > capacity[slot]) {
        free(block[slot]);
        block[slot] = malloc(bytes);
        capacity[slot] = block[slot] != NULL ? bytes : 0;
        if (block[slot] == NULL) return calloc(1, bytes);
    }
    if (zero) memset(block[slot], 0, bytes);
    return block[slot];
}

// Same interleaved cosine/sine recurrence as make_oscillator, fused with the]=])
tn_replace(features "        float *axis = (float *)calloc((size_t)ndec, sizeof(*axis));"
                 "        float *axis = (float *)tn_band_scratch(0, (size_t)(ndec > 0 ? ndec : 1) * sizeof(*axis), 0);")
tn_replace(features "        float *mix_r = (float *)calloc((size_t)(nw + nx), sizeof(*mix_r));"
                 "        float *mix_r = (float *)tn_band_scratch(1, (size_t)(nw + nx) * sizeof(*mix_r), 0);")
tn_replace(features "        float *mix_i = (float *)calloc((size_t)(nw + nx), sizeof(*mix_i));"
                 "        float *mix_i = (float *)tn_band_scratch(2, (size_t)(nw + nx) * sizeof(*mix_i), 0);")
tn_replace(features "\n        conv_r = (float *)calloc((size_t)(nx > 0 ? nx : 1), sizeof(*conv_r));"
                 "\n        conv_r = (float *)tn_band_scratch(3, (size_t)(nx > 0 ? nx : 1) * sizeof(*conv_r), 0);")
tn_replace(features "\n        conv_i = (float *)calloc((size_t)(nx > 0 ? nx : 1), sizeof(*conv_i));"
                 "\n        conv_i = (float *)tn_band_scratch(4, (size_t)(nx > 0 ? nx : 1) * sizeof(*conv_i), 0);")
tn_replace(features "\n            conv_r = (float *)calloc((size_t)(nx > 0 ? nx : 1), sizeof(*conv_r));"
                 "\n            conv_r = (float *)tn_band_scratch(5, (size_t)(nx > 0 ? nx : 1) * sizeof(*conv_r), 0);")
tn_replace(features "\n            conv_i = (float *)calloc((size_t)(nx > 0 ? nx : 1), sizeof(*conv_i));"
                 "\n            conv_i = (float *)tn_band_scratch(6, (size_t)(nx > 0 ? nx : 1) * sizeof(*conv_i), 0);")
tn_replace(features "        dec_r = (float *)calloc((size_t)nwork, sizeof(*dec_r));"
                 "        dec_r = (float *)tn_band_scratch(7, (size_t)nwork * sizeof(*dec_r), 1);")
tn_replace(features "        dec_i = (float *)calloc((size_t)nwork, sizeof(*dec_i));"
                 "        dec_i = (float *)tn_band_scratch(8, (size_t)nwork * sizeof(*dec_i), 1);")
tn_replace(features "        free(mix_r);\n        free(mix_i);\n" "")
tn_replace(features "            free(mix_r);\n            free(mix_i);\n" "")
tn_replace(features "        free(dec_r);\n        free(dec_i);\n" "")
tn_replace(features "        free(axis);\n" "")
# estimate_ifsnr keeps every decim-th sample of the windowed heterodyne, so the
# full-rate overlap-add convolution computes up to decim-1 out of decim bins for
# nothing. Form the needed samples directly instead: the window is reversed once
# per band and each retained bin is a short contiguous dot product.
tn_replace(features "// per-band SNR/IF worker: heterodyne, decimate, x2as, interpolate, scatter."
[=[// Only every decim-th sample of the windowed heterodyne is consumed. Build the
// windowed signal once, then evaluate those bins directly; the arithmetic is
// the same linear convolution, without the discarded full-rate outputs.
static void tn_heterodyne_decimated(const float *x, int nx, float carrier,
                                    const float *win, int nw, int decim, int ndec,
                                    float *out_r, float *out_i) {
    const int head = nw + 2;
    const size_t span = (size_t)nx + 2u * (size_t)head + 2u;
    float *re = (float *)tn_band_scratch(9, span * sizeof(float), 0);
    float *im = (float *)tn_band_scratch(10, span * sizeof(float), 0);
    if (re == NULL || im == NULL) return;
    memset(re, 0, (size_t)head * sizeof(float));
    memset(im, 0, (size_t)head * sizeof(float));
    memset(re + head + nx, 0, (size_t)head * sizeof(float));
    memset(im + head + nx, 0, (size_t)head * sizeof(float));
    tn_heterodyne(x, nx, carrier, re + head, im + head);
    const int start = head + nw / 2 - nw + 1;
    for (int i = 0; i < ndec; ++i) {
        const float *p = re + start + i * decim;
        const float *q = im + start + i * decim;
        // The windowed sum can cancel heavily, but an explicit SIMD reduction
        // lets the compiler keep its own partial sums and use the full vector
        // width; this loop is the dominant cost of the estimate.
        float sr = 0.f, si = 0.f;
        #pragma omp simd reduction(+:sr,si)
        for (int k = 0; k < nw; ++k) {
            sr += win[k] * p[k];
            si += win[k] * q[k];
        }
        out_r[i] = sr;
        out_i[i] = si;
    }
}

// per-band SNR/IF worker: heterodyne, decimate, x2as, interpolate, scatter.]=])
tn_replace(features [=[        tn_nebula_conv w_conv = tn_nebula_conv_make(w, nw);
        tn_nebula_conv wf_conv = tn_nebula_conv_make(wf, nh);]=]
[=[        tn_nebula_conv wf_conv = tn_nebula_conv_make(wf, nh);]=])
tn_replace(features [=[        float *mix_r = (float *)tn_band_scratch(1, (size_t)(nw + nx) * sizeof(*mix_r), 0);
        float *mix_i = (float *)tn_band_scratch(2, (size_t)(nw + nx) * sizeof(*mix_i), 0);
        float *conv_r, *conv_i;
        float *dec_r, *dec_i;]=]
[=[        float *dec_r, *dec_i;]=])
tn_replace(features [=[        tn_heterodyne(x, nx, carrier, mix_r, mix_i);
        conv_r = (float *)tn_band_scratch(3, (size_t)(nx > 0 ? nx : 1) * sizeof(*conv_r), 0);
        conv_i = (float *)tn_band_scratch(4, (size_t)(nx > 0 ? nx : 1) * sizeof(*conv_i), 0);
        tn_nebula_conv_run_shift(&w_conv, mix_r, mix_i, nx, nw / 2, conv_r, conv_i);
        mix_r = conv_r;
        mix_i = conv_i;
        int nwork = ndec + 256;
        dec_r = (float *)tn_band_scratch(7, (size_t)nwork * sizeof(*dec_r), 1);
        dec_i = (float *)tn_band_scratch(8, (size_t)nwork * sizeof(*dec_i), 1);
        for (i = 0; i < ndec; ++i) {
            dec_r[i] = mix_r[i * decim];
            dec_i[i] = mix_i[i * decim];
        }]=]
[=[        int nwork = ndec + 256;
        dec_r = (float *)tn_band_scratch(7, (size_t)nwork * sizeof(*dec_r), 1);
        dec_i = (float *)tn_band_scratch(8, (size_t)nwork * sizeof(*dec_i), 1);
        // Reverse the analysis window once so each retained bin reads forward.
        for (int r = 0; r < nw / 2; ++r) {
            const float tmp = w[r];
            w[r] = w[nw - 1 - r];
            w[nw - 1 - r] = tmp;
        }
        tn_heterodyne_decimated(x, nx, carrier, w, nw, decim, ndec, dec_r, dec_i);]=])
tn_replace(features [=[            tn_heterodyne(x, nx, (float)(carrier + carrier), mix_r, mix_i);
            conv_r = (float *)tn_band_scratch(5, (size_t)(nx > 0 ? nx : 1) * sizeof(*conv_r), 0);
            conv_i = (float *)tn_band_scratch(6, (size_t)(nx > 0 ? nx : 1) * sizeof(*conv_i), 0);
            tn_nebula_conv_run_shift(&w_conv, mix_r, mix_i, nx, nw / 2, conv_r, conv_i);
            mix_r = conv_r;
            mix_i = conv_i;
            for (i = 0; i < ndec; ++i) {
                dec_r[i] = mix_r[i * decim];
                dec_i[i] = mix_i[i * decim];
            }]=]
[=[            tn_heterodyne_decimated(x, nx, (float)(carrier + carrier), w, nw, decim, ndec, dec_r, dec_i);]=])
tn_replace(features "        free(w_conv.kernel);\n" "")
set(nebula_features_override "${features}")

file(READ "${nebula_root}/nebula/nebula_model.c" model)
string(REPLACE "\r\n" "\n" model "${model}")
tn_replace(model "static float *neb_condgmm_density" "static void neb_condgmm_density")
tn_replace(model "const float *axis, int width)" "const float *axis, int width, float *density)")
tn_replace(model [=[    float *weights = (float *)calloc((size_t)nmix, sizeof(float));
    float *means = (float *)calloc((size_t)nmix, sizeof(float));
    float *sigma = (float *)calloc((size_t)nmix, sizeof(float));
    float *logw = (float *)calloc((size_t)nmix, sizeof(float));
    float *y = (float *)calloc((size_t)(2 * g->nord), sizeof(float));]=]
[=[    float weights[NEB_NMIX], means[NEB_NMIX], sigma[NEB_NMIX], logw[NEB_NMIX];
    float y[NEB_NFEAT];]=])
tn_replace(model "    float *density = (float *)calloc((size_t)width, sizeof(float));"
                 "    memset(density, 0, (size_t)width*sizeof(float));")
tn_replace(model [=[    free(weights);
    free(means);
    free(sigma);
    free(logw);
    free(y);
    return density;
]=] "")
tn_replace(model "            float *density = neb_condgmm_density(bands[b], feat, axis, width);"
                 "            float density[NEB_NF0];\n            neb_condgmm_density(bands[b], feat, axis, width, density);")
tn_replace(model "            free(density);\n" "")
# The original density already returns exactly zero for q > 50. On the
# monotonic F0 axis, narrow the scan to a conservative q <= 52 interval;
# the original q <= 50 gate and per-bin mixture addition remain in force.
tn_replace(model [=[        float invsigma = 1.0 / sigma[i];
        for (int k = 0; k < width; ++k) {]=]
[=[        float invsigma = 1.0 / sigma[i];
        const double radius = sqrt(52.0 * sigma[i]);
        const double lo_hz = means[i]-radius, hi_hz = means[i]+radius;
        int first = 0, last = width;
        while (first < last) {
            const int middle = first+(last-first)/2;
            if (axis[middle] < lo_hz) first = middle+1;
            else last = middle;
        }
        int past = width;
        while (last < past) {
            const int middle = last+(past-last)/2;
            if (axis[middle] <= hi_hz) last = middle+1;
            else past = middle;
        }
        // Align both ends to the same group of bins as the upstream loop.
        first &= ~3;
        past = (past+3)&~3;
        if (past > width) past = width;
        for (int k = first; k < past; ++k) {]=])
tn_replace(model "#pragma omp parallel for" "#pragma omp parallel for schedule(dynamic, 4)")
# The 5x5 triangular solves run 16 times per frame and band. Keep the library's
# exact forward/back substitution order, but inline it so y stays in registers.
tn_replace(model "        cig_lusolve(g->inv11[i], y, nfeat);"
[=[        {
            float *lu = g->inv11[i];
            int j, k;
            for (j = 0; j < nfeat; ++j)
                for (k = j + 1; k < nfeat; ++k)
                    y[k] -= lu[k + j * nfeat] * y[j];
            for (j = nfeat - 1; j >= 0; --j) {
                y[j] /= lu[j + j * nfeat];
                for (k = 0; k < j; ++k)
                    y[k] -= lu[k + j * nfeat] * y[j];
            }
        }]=])
# Rewrite the fixed-5-feature conditional density so the bin loop can be
# vectorized: the same fasterpow2 bit trick expressed as plain float ops, and
# a `simd` hint. Every bin is still computed and accumulated independently, so
# the values are unchanged.
tn_replace(model "static void neb_condgmm_density"
[=[static inline float tn_fasterexp(float p)
{
    const float q = 1.442695040f * p;
    const float clipp = (q < -126.0f) ? -126.0f : q;
    const uint32_t bits = (uint32_t)(8388608.0f * (clipp + 126.94269504f));
    float f;
    memcpy(&f, &bits, sizeof(f));
    return f;
}

static void neb_condgmm_density]=])
tn_replace(model "const float *axis, int width, float *density)"
                 "const float *restrict axis, int width, float *restrict density)")
tn_replace(model "            float value = q <= 50.0 ? fasterexp(-0.5 * q) : 0.0f;"
                 "            float value = q <= 50.0 ? tn_fasterexp(-0.5 * q) : 0.0f;")
tn_replace(model "        for (int k = first; k < past; ++k) {"
                 "        #pragma omp simd\n        for (int k = first; k < past; ++k) {")
# Per-mixture quantities that do not depend on the frame: the log weight, the
# inverse and square-root of the conditional sigma, and the scan radius. They
# are computed once when the band is loaded instead of for every frame.
tn_replace(model "neb_gmm *neb_gmm_load_band(int band)"
[=[static float tn_sigma_inv[NEB_NBAND][NEB_NMIX];
static float tn_sigma_sqrt[NEB_NBAND][NEB_NMIX];
static double tn_sigma_radius[NEB_NBAND][NEB_NMIX];
static double tn_logw0[NEB_NBAND][NEB_NMIX];

neb_gmm *neb_gmm_load_band(int band)]=])
tn_replace(model [=[        g->condsigma[m] = cs;]=]
[=[        g->condsigma[m] = cs;
        {
            const float sg = cs < 4.0f ? 4.0f : cs;
            tn_sigma_inv[band][m] = 1.0f / sg;
            tn_sigma_sqrt[band][m] = (float)sqrt(6.283185307179586 * (double)sg);
            tn_sigma_radius[band][m] = sqrt(52.0 * (double)sg);
            tn_logw0[band][m] = (double)fastlog(g->w[m])
                              - (double)(g->nord - 1) * 1.837877066409345;
        }]=])
tn_replace(model "const float *restrict axis, int width, float *restrict density)"
                 "const float *restrict axis, int width, float *restrict density, int band)")
tn_replace(model "            neb_condgmm_density(bands[b], feat, axis, width, density);"
                 "            neb_condgmm_density(bands[b], feat, axis, width, density, b);")
tn_replace(model "        float scale = weights[i] / sqrt(6.283185307179586 * sigma[i]);"
                 "        float scale = weights[i] / tn_sigma_sqrt[band][i];")
tn_replace(model "        float invsigma = 1.0 / sigma[i];"
                 "        float invsigma = tn_sigma_inv[band][i];")
tn_replace(model "        const double radius = sqrt(52.0 * sigma[i]);"
                 "        const double radius = tn_sigma_radius[band][i];")
tn_replace(model [=[        float value = fastlog(g->w[i]) - normalizer
                    + (quad + g->det11[i]) * -0.5;]=]
[=[        float value = (float)(tn_logw0[band][i]
                             + (double)((quad + g->det11[i]) * -0.5));]=])
set(nebula_model_override "${model}")

# Each Lanczos output sample is independent. Leave the per-sample tap
# reduction intact and split only the outer loop on substantial recordings.
file(READ "${nebula_root}/external/ciglet/ciglet.c" ciglet)
string(REPLACE "\r\n" "\n" ciglet "${ciglet}")
tn_replace(ciglet [=[  for(int i = 0; i < *ny; i ++) {
    FP_TYPE ix = (FP_TYPE)i / ratio;]=]
[=[  #pragma omp parallel for if(*ny >= 8192) schedule(static)
  for(int i = 0; i < *ny; i ++) {
    FP_TYPE ix = (FP_TYPE)i / ratio;]=])
# The Lanczos kernel divided by pi twice and by x twice; one division of the
# product removes three dependent divisions per tap. Same closed form.
tn_replace(ciglet "  return a * sin_2(M_PI * x) * sin_2(M_PI * x / a) / M_PI / M_PI / x / x;"
                 "  return a * sin_2(M_PI * x) * sin_2(M_PI * x / a) / (M_PI * M_PI * x * x);")
# The Lanczos weights depend only on the fractional read position. Tabulate
# them once per load instead of running two sin approximations and a division
# for each of the eight taps of every output sample.
tn_replace(ciglet [=[  #pragma omp parallel for if(*ny >= 8192) schedule(static)
  for(int i = 0; i < *ny; i ++) {
    FP_TYPE ix = (FP_TYPE)i / ratio;
    int ix_base = floor(ix);
    FP_TYPE a = ix - ix_base;
    for(int j = -3; j <= 4; j ++)
      if(ix_base + j > 0 && ix_base + j < nx)
        y[i] += interp_kernel(a - j, 4) * xlp[ix_base + j];
  }]=]
[=[  #pragma omp parallel for if(*ny >= 8192) schedule(static)
  for(int i = 0; i < *ny; i ++) {
    FP_TYPE ix = (FP_TYPE)i / ratio;
    int ix_base = floor(ix);
    FP_TYPE a = ix - ix_base;
    const float *wt = tn_lanczos_weights(a);
    for(int j = -3; j <= 4; j ++)
      if(ix_base + j > 0 && ix_base + j < nx)
        y[i] += wt[j + 3] * xlp[ix_base + j];
  }]=])
tn_replace(ciglet [=[FP_TYPE* cig_rresample(FP_TYPE* x, int nx, FP_TYPE ratio, int* ny) {]=]
[=[#define TN_LANCZOS_PHASES 4096
static float tn_lanczos_table[TN_LANCZOS_PHASES + 1][8];
static int tn_lanczos_ready;

// Built on first resample rather than at load: the tool starts a process per
// note, and the warm path never resamples, so an eager table would only cost
// launch time.
static void tn_init_lanczos(void)
{
  if(tn_lanczos_ready) return;
  for(int p = 0; p <= TN_LANCZOS_PHASES; p ++) {
    const FP_TYPE a = (FP_TYPE)p / (FP_TYPE)TN_LANCZOS_PHASES;
    for(int t = 0; t < 8; t ++)
      tn_lanczos_table[p][t] = interp_kernel(a - (t - 3), 4);
  }
  tn_lanczos_ready = 1;
}

static const float *tn_lanczos_weights(FP_TYPE a)
{
  int p = (int)(a * (FP_TYPE)TN_LANCZOS_PHASES + (FP_TYPE)0.5);
  if(p < 0) p = 0;
  if(p > TN_LANCZOS_PHASES) p = TN_LANCZOS_PHASES;
  return tn_lanczos_table[p];
}

FP_TYPE* cig_rresample(FP_TYPE* x, int nx, FP_TYPE ratio, int* ny) {
  tn_init_lanczos();]=])
set(ciglet_override "${ciglet}")

message(STATUS "Nebula: reusing overlap-add workspace, removing GMM allocations and parallelizing long resampling")
