// SPDX-License-Identifier: GPL-3.0-only
// Used inside the pinned Nebula IF/SNR translation unit. Keep Ooura's FFT,
// the same overlap-add block grid, and each output bin's addition order.
#ifndef TN_NEBULA_CONVOLUTION_H
#define TN_NEBULA_CONVOLUTION_H

void cdft(int n, int isgn, float* data);

typedef struct {
    int nfft, nh;
    float* kernel;
    float* block;
} tn_nebula_conv;

static tn_nebula_conv tn_nebula_conv_make(const float* h, int nh) {
    int nfft = 4;
    while (nfft/4 < nh) nfft *= 2;
    float* data = (float*)calloc((size_t)4*nfft, sizeof(float));
    for (int i = 0; i < nh; ++i) data[2*i] = h[i];
    cdft(2*nfft, -1, data);
    return (tn_nebula_conv){nfft, nh, data, data+2*nfft};
}

static void tn_nebula_conv_run(tn_nebula_conv* p, const float* xr, const float* xi,
                               int nx, float* yr, float* yi) {
    const int nfft = p->nfft, step = nfft-p->nh;
    const int outlen = nx+p->nh-1;
    float* restrict a = p->block;
    const float* restrict h = p->kernel;
    memset(yr, 0, (size_t)outlen*sizeof(float));
    if (yi) memset(yi, 0, (size_t)outlen*sizeof(float));
    const float inv_nfft = 1.f/nfft;
    for (int offset = 0; offset <= nx; offset += step) {
        const int take = nx-offset < step ? nx-offset : step;
        if (xi) {
            for (int i = 0; i < take; ++i) {
                a[2*i] = xr[offset+i];
                a[2*i+1] = xi[offset+i];
            }
        } else {
            for (int i = 0; i < take; ++i) {
                a[2*i] = xr[offset+i];
                a[2*i+1] = 0.f;
            }
        }
        memset(a+2*take, 0, (size_t)2*(nfft-take)*sizeof(float));
        cdft(2*nfft, -1, a);
        for (int i = 0; i < 2*nfft; i += 2) {
            const float r = a[i], j = a[i+1];
            a[i] = r*h[i]-j*h[i+1];
            a[i+1] = h[i+1]*r+h[i]*j;
        }
        cdft(2*nfft, 1, a);
        const int count = outlen-offset < nfft ? outlen-offset : nfft;
        for (int i = 0; i < count; ++i) yr[offset+i] += a[2*i]*inv_nfft;
        if (yi)
            for (int i = 0; i < count; ++i) yi[offset+i] += a[2*i+1]*inv_nfft;
    }
}

// Same block grid and per-bin addition order, but keep only the window the
// caller actually reads. The caller no longer allocates a padded result or
// copies the centred region out of it.
static void tn_nebula_conv_run_shift(tn_nebula_conv* p, const float* xr, const float* xi,
                                     int nx, int shift, float* yr, float* yi) {
    const int nfft = p->nfft, step = nfft-p->nh;
    const int outlen = nx+p->nh-1;
    float* restrict a = p->block;
    const float* restrict h = p->kernel;
    const float inv_nfft = 1.f/nfft;
    if (nx > 0) {
        memset(yr, 0, (size_t)nx*sizeof(float));
        if (yi) memset(yi, 0, (size_t)nx*sizeof(float));
    }
    for (int offset = 0; offset <= nx; offset += step) {
        const int take = nx-offset < step ? nx-offset : step;
        if (xi) {
            for (int i = 0; i < take; ++i) {
                a[2*i] = xr[offset+i];
                a[2*i+1] = xi[offset+i];
            }
        } else {
            for (int i = 0; i < take; ++i) {
                a[2*i] = xr[offset+i];
                a[2*i+1] = 0.f;
            }
        }
        memset(a+2*take, 0, (size_t)2*(nfft-take)*sizeof(float));
        cdft(2*nfft, -1, a);
        for (int i = 0; i < 2*nfft; i += 2) {
            const float r = a[i], j = a[i+1];
            a[i] = r*h[i]-j*h[i+1];
            a[i+1] = h[i+1]*r+h[i]*j;
        }
        cdft(2*nfft, 1, a);
        const int count = outlen-offset < nfft ? outlen-offset : nfft;
        for (int i = 0; i < count; ++i) {
            const int at = offset+i-shift;
            if (at >= 0 && at < nx) yr[at] += a[2*i]*inv_nfft;
        }
        if (yi)
            for (int i = 0; i < count; ++i) {
                const int at = offset+i-shift;
                if (at >= 0 && at < nx) yi[at] += a[2*i+1]*inv_nfft;
            }
    }
}

#endif
