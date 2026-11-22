# SPDX-License-Identifier: GPL-3.0-only
"""Shared harmonic/noise fixtures for the flag and mask regressions."""

from __future__ import annotations

import numpy as np

HIGH_BAND = (5800., 10100.)


def fade(t, start, end, ramp):
    u = np.clip(np.minimum(t-start, end-t) / ramp, 0., 1.)
    return np.sin(np.pi*u/2.)**2


def rms(x):
    return float(np.sqrt(np.mean(np.square(x, dtype=np.float64))))


def colored_noise(t, sr, lo, hi, amplitude, seed):
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(len(t), 1/sr)
    spec = np.fft.rfft(rng.standard_normal(len(t)))
    spec[(f < lo) | (f > hi)] = 0.
    z = np.fft.irfft(spec, len(t))
    return z*(amplitude/rms(z))


def fixture(sr=44100, f0=149., *, weak_h1=False, jitter_us=0., high=True,
            high_rms=.031, noise_rms=.020, seed=5617):
    """Aligned harmonic and noise stems with high partials inside the noise band."""
    t = np.arange(round(1.35*sr))/sr
    phase = 2*np.pi*f0*t + .37
    if jitter_us > 0.:
        cycle = t*f0
        i = cycle.astype(int)
        a = cycle-i
        shifts = np.random.default_rng(seed+3).uniform(-1., 1., int(i[-1])+2) * jitter_us * 1e-6
        w = a*a*(3-2*a)
        phase += 2*np.pi*f0*((1-w)*shifts[i] + w*shifts[i+1])
        phase += f0*.004/4.1*np.sin(2*np.pi*4.1*t)
    amps = (.154*(.025 if weak_h1 else 1.), .083, .060, .034, .020, .014)
    h = sum(a*np.cos(k*phase + .18*k**1.36) for k, a in enumerate(amps, 1))
    h *= fade(t, .15, 1.21, .030)
    if high:
        order = np.arange(int(np.ceil(6100/f0)), int(9800/f0)+1)
        upper = sum((.55+.45*np.cos(.32*k))*np.cos(k*phase + .053*k**1.62) for k in order)
        upper *= high_rms/rms(upper)
        h += upper*fade(t, .30, 1.06, .060)
    noise = colored_noise(t, sr, 5650, 10350, noise_rms, seed) * fade(t, .025, 1.31, .030)
    noise += colored_noise(t, sr, 800, 3200, .003, seed+1) * fade(t, .025, 1.31, .045)
    return h.astype(np.float32), noise.astype(np.float32)


def band_coefficients(x, sr, span=(.43, .93), band=HIGH_BAND):
    a, b = (round(t*sr) for t in span)
    x = np.asarray(x[a:b], dtype=np.float64)
    f = np.fft.rfftfreq(len(x), 1/sr)
    return np.fft.rfft(x*np.hanning(len(x)))[(f >= band[0]) & (f < band[1])]
