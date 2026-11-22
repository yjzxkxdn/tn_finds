# SPDX-License-Identifier: GPL-3.0-only
"""Audio + command-line regressions; requires NumPy and SoundFile (tests only)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
import soundfile as sf

SR = 48000
ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / 'tn_finds.exe'


def vowel(hz: float, seconds: float, sr: int = SR) -> np.ndarray:
    """Harmonic fixture with a specified vowel-shaped spectrum and aligned phases."""
    h = np.arange(hz, sr / 2, hz)
    a = (.95 * np.exp(-.5 * ((h - 650) / 230) ** 2)
         + .66 * np.exp(-.5 * ((h - 2050) / 280) ** 2)
         + .22 * np.exp(-.5 * ((h - 3200) / 320) ** 2)
         + .008 * np.exp(-h / 2900))
    t = np.arange(round(seconds * sr)) / sr
    y = np.zeros_like(t)
    for f, amp in zip(h, a):
        y += amp * np.cos(2 * np.pi * f * t)
    return (.45 * y / np.max(np.abs(y))).astype(np.float32)


def compound(sr: int = SR) -> np.ndarray:
    """120 ms /s/-like noise, 0.98 s vowel with aspiration, 100 ms release."""
    r = np.random.default_rng(856)
    n = round(1.2 * sr)
    raw = r.normal(size=n + 1).astype(np.float32)
    x = .018 * (raw[1:] - raw[:-1])
    a, b = round(.12 * sr), round(1.1 * sr)
    body = vowel(200., (b - a) / sr, sr)
    fade = round(.015 * sr)
    body[:fade] *= np.sin(np.linspace(0., np.pi / 2, fade)) ** 2
    body[-fade:] *= np.sin(np.linspace(np.pi / 2, 0., fade)) ** 2
    x[a:b] = .18 * x[a:b] + body
    x[b:] *= np.linspace(.7, 0, n - b)
    return x.astype(np.float32)


def boundary_fixture(phase: float = 0., sr: int = SR):
    """Frication → smooth 160-Hz voiced onset → release; no injected thump."""
    t = np.arange(round(.92*sr)) / sr
    u = t - .15
    env = np.clip(np.minimum(u, .70 - t) / .008, 0., 1.)
    angle = 2*np.pi*160*u + phase
    body = .20*np.cos(angle) + .06*np.cos(2*angle + .4) + .02*np.cos(3*angle - .6)
    noise = .01*np.diff(np.random.default_rng(1242).normal(size=len(t)+1))
    return np.where((t >= .15) & (t < .70), (.5-.5*np.cos(np.pi*env))*body, noise).astype(np.float32)


def subvoice_lp(x: np.ndarray, sr: int):
    """Offline diagnostic only: 40 Hz passband, 65 Hz stopband, no delay."""
    p = np.pad(x.astype(np.float64), sr//4)
    f = np.fft.rfftfreq(len(p), 1/sr)
    g = .5 + .5*np.cos(np.pi*np.clip((f-40)/25, 0, 1))
    y = np.fft.irfft(np.fft.rfft(p)*g, len(p))
    return y[sr//4:sr//4+len(x)]


def stripe_fixture(phase=0., sr=SR):
    """Band-limited UV–V–UV stimulus with smooth input joins and known F0."""
    n = round(.92 * sr)
    t = np.arange(n) / sr
    p = 2*np.pi*160*(t-.15) + phase
    body = .2*np.cos(p) + .075*np.cos(2*p+.4) + .04*np.cos(3*p-.8)
    env = .5 - .5*np.cos(np.pi*np.clip(np.minimum(t-.15, .70-t)/.008, 0, 1))
    rng = np.random.default_rng(58)
    spec = np.fft.rfft(rng.normal(size=n))
    f = np.fft.rfftfreq(n, 1/sr)
    spec[(f < 900) | (f > 1700)] = 0
    breath = np.fft.irfft(spec, n)
    breath *= .007 / np.sqrt(np.mean(breath**2))
    return (env*body + (1-env)*breath).astype(np.float32)


def stripe_band(x, sr=SR):
    """Offline short-impulse diagnostic: smooth 3–12 kHz band, no delay."""
    pad = sr//16
    p = np.pad(x.astype(np.float64), pad)
    f = np.fft.rfftfreq(len(p), 1/sr)
    up = .5 - .5*np.cos(np.pi*np.clip((f-2700)/300, 0, 1))
    down = .5 + .5*np.cos(np.pi*np.clip((f-11700)/300, 0, 1))
    return np.fft.irfft(np.fft.rfft(p)*up*down, len(p))[pad:pad+len(x)]


MIXED_WINDOWS_S = ((.33, .55), (.65, .95), (1.05, 1.27))


def mixed_fixture(sr=SR, *, f0=149., phase=.37, h1_scale=1., noise_rms=.055, seed=1749):
    """Return separate low harmonic and high noise stems, both with smooth envelopes."""
    t = np.arange(round(1.60 * sr)) / sr
    angle = 2 * np.pi * f0 * t + phase
    amplitudes = (.145 * h1_scale, .078, .056, .033, .018, .011)
    harmonic = sum(a * np.cos(k * angle + .29 * k**1.3) for k, a in enumerate(amplitudes, 1))
    v = np.clip(np.minimum(t - .22, 1.38 - t) / .024, 0., 1.)
    harmonic *= np.sin(.5 * np.pi * v)**2

    spectrum = np.fft.rfft(np.random.default_rng(seed).standard_normal(len(t)))
    frequencies = np.fft.rfftfreq(len(t), 1. / sr)
    spectrum[(frequencies < 5000) | (frequencies > 9000)] = 0.
    noise = np.fft.irfft(spectrum, len(t))
    noise *= noise_rms / np.sqrt(np.mean(noise**2))
    quieter = np.clip(np.minimum(t - .57, 1.03 - t) / .055, 0., 1.)
    noise *= 1. - .65 * np.sin(.5 * np.pi * quieter)**2
    return harmonic.astype(np.float32), noise.astype(np.float32)


def mixed_f0_points(*, f0=149., gap=None, all_zero=False):
    """Absolute source seconds for mixed_fixture, optionally with missing F0 labels."""
    if all_zero:
        return [(0., 0.), (1.6, 0.)]
    points = [(0., 0.), (.2199, 0.), (.22, f0)]
    if gap is not None:
        a, b = gap
        points.extend([(a - .0001, f0), (a, 0.), (b - .0001, 0.), (b, f0)])
    return points + [(1.3799, f0), (1.38, 0.), (1.6, 0.)]


def harmonic_grid(x, source_hz, target_hz, sr=SR, *, max_hz=1200.):
    """Joint least-squares fit to the distinguishable lines of two harmonic grids."""
    t = np.arange(len(x), dtype=np.float64) / sr
    separation = sr / len(x)
    source = source_hz * np.arange(1, int(max_hz / source_hz) + 1)
    target_orders = np.arange(1, int(max_hz / target_hz) + 1)
    target = target_hz * target_orders
    distinguishable = np.abs(source[:, None] - target) >= separation
    use_source, use_target = distinguishable.all(axis=1), distinguishable.all(axis=0)
    frequencies = np.r_[source[use_source], target[use_target]]
    if not use_source.any() or not use_target.any():
        raise ValueError('window must distinguish at least one line in each grid')
    phase = 2 * np.pi * t[:, None] * frequencies
    w = np.sqrt(np.hanning(len(x)))
    basis = np.column_stack((np.ones(len(t)), np.cos(phase), np.sin(phase)))
    coeff = np.linalg.lstsq(basis * w[:, None], x * w, rcond=None)[0]
    size = len(frequencies)
    line_power = (coeff[1:1+size]**2 + coeff[1+size:]**2) / 2.
    source_count = int(np.count_nonzero(use_source))
    source_rms = float(np.sqrt(np.sum(line_power[:source_count])))
    target_power = line_power[source_count:]
    target_rms = float(np.sqrt(np.sum(target_power)))
    odd_rms = float(np.sqrt(np.sum(target_power[target_orders[use_target] % 2 == 1])))
    ratio = 20. * np.log10(max(1e-12, target_rms) / max(1e-12, source_rms))
    return {'source_rms': source_rms, 'target_rms': target_rms, 'target_odd_rms': odd_rms,
            'target_over_source_db': float(ratio)}


def paired_noise_metrics(y_plus, y_minus, noise, sr=SR, *, window=None, band=None):
    """Compare the sign-sensitive high-band response with the corresponding dry noise."""
    region = slice(None) if window is None else slice(round(window[0] * sr), round(window[1] * sr))
    def bandpass(x):
        if band is None:
            return stripe_band(x, sr)
        low, high = band
        pad = sr//16
        p = np.pad(x.astype(np.float64), pad)
        f = np.fft.rfftfreq(len(p), 1/sr)
        rise = .5 - .5*np.cos(np.pi*np.clip((f-low+200.)/200., 0., 1.))
        fall = .5 + .5*np.cos(np.pi*np.clip((f-high)/200., 0., 1.))
        return np.fft.irfft(np.fft.rfft(p)*rise*fall, len(p))[pad:pad+len(x)]

    observed = bandpass((y_plus.astype(np.float64) - y_minus) / 2.)[region]
    expected = bandpass(noise)[region]
    expected_norm = float(np.linalg.norm(expected))
    observed_norm = float(np.linalg.norm(observed))
    dot = float(np.dot(expected, observed))
    scale = dot / max(expected_norm**2, 1e-24)
    residual = float(np.linalg.norm(observed - scale * expected))
    db = lambda v: float(20. * np.log10(max(v, 1e-12) / max(expected_norm, 1e-12)))
    return {'gain_db': db(observed_norm), 'correlation': dot / max(expected_norm * observed_norm, 1e-24),
            'residual_db': db(residual)}


def estimate_hz(x: np.ndarray, sr: int, minimum: float = 55, maximum: float = 1300) -> tuple[float, float]:
    """Independent normalized autocorrelation, picking the first strong peak."""
    x = x.astype(np.float64) - np.mean(x)
    if np.max(np.abs(x)) < 1e-10:
        return 0., 0.
    n = len(x)
    l, r = max(2, int(sr / maximum)), min(n // 2 - 1, int(np.ceil(sr / minimum)))
    p = np.arange(r + 2)
    fft_n = 1 << int(2 * n - 1).bit_length()
    cross = np.fft.irfft(np.abs(np.fft.rfft(x, fft_n)) ** 2)[:r + 2]
    ss = np.r_[0., np.cumsum(x * x)]
    corr = cross / np.sqrt(np.maximum(1e-30, ss[n - p] * (ss[n] - ss[p])))
    peaks, = np.nonzero((corr[l:r+1] >= corr[l-1:r]) & (corr[l:r+1] > corr[l+1:r+2]))
    if not len(peaks):
        return 0., 0.
    peaks += l
    top = float(np.max(corr[peaks]))
    hits = peaks[corr[peaks] >= max(.52, top - .045)]
    if not len(hits):
        return 0., top
    k = int(hits[0])
    d = .5 * (corr[k-1] - corr[k+1]) / min(-1e-10, corr[k-1] - 2*corr[k] + corr[k+1])
    return float(sr / (k + np.clip(d, -.5, .5))), float(corr[k])


def f2_peak(x: np.ndarray, sr: int, hz: float = 200.) -> float:
    t = np.arange(len(x)) / sr
    x = x * np.hanning(len(x))
    harmonics = np.arange(hz, 3700., hz)
    a = np.array([abs(np.dot(x, np.exp(-2j * np.pi * f * t))) for f in harmonics])
    a = np.log(np.maximum(1e-15, a))
    inside = np.where((harmonics >= 1400) & (harmonics <= 2850))[0]
    j = int(inside[np.argmax(a[inside])])
    d = .5 * (a[j-1] - a[j+1]) / min(-1e-10, a[j-1] - 2*a[j] + a[j+1])
    return float(harmonics[j] + hz * np.clip(d, -.5, .5))


def phase_fixture():
    # Nonzero glottal and resonant-filter phase. An all-in-phase sum is a
    # weak control for the effect of relative phase on waveform shape.
    src, dst = 100, 85
    h = np.arange(1, 62)
    f = src * h
    poles = [(540., 95.), (1440., 180.), (2520., 260.)]
    h_resp = np.prod([1 + 1.6 * fc**2 / (fc**2 - f**2 + 2j * bw * f) for fc, bw in poles], axis=0)
    a = h_resp / h * np.exp(-1j * .35 * h / (1 + .06*h*h))
    x = np.zeros(int(1.6*SR))
    t = np.arange(len(x), dtype=np.float64) / SR
    for g, c in zip(f, a):
        x += np.real(c * np.exp(2j * np.pi * g * t))
    return (x * (.44 / np.max(np.abs(x)))).astype(np.float32), src, dst


def relative_shape_metrics(source, measured, src, dst, control=None):
    # Integer grids on the exact one-second interval. A time translation
    # adds k*delta to harmonic k; subtracting k*phi_1 cancels that slope.
    def coeff(y, f):
        return np.fft.rfft(y[SR//4:SR//4+SR])[f*np.arange(1, 49)]
    a, b = coeff(source, src), coeff(measured, dst)
    weight = np.abs(a)**2
    hj = np.arange(1, len(a)+1)
    phase = np.angle(b * np.conj(a))
    relative = np.angle(np.exp(1j*(phase - phase[0]*hj)))
    useful = (hj >= 2) & (np.abs(a) >= .035*np.abs(a).max()) & (np.abs(b) >= .02*np.abs(b).max())
    importance = np.clip(weight[useful] / np.max(weight), .01, .6)
    out = {'partial_count': int(np.count_nonzero(useful)),
           'relative_phase_rms_degrees': float(np.rad2deg(np.sqrt(np.average(relative[useful]**2, weights=importance))))}
    if control is not None:
        c = coeff(control, dst)
        d = 20*np.log10(np.maximum(np.abs(b[useful]), 1e-12) / np.maximum(np.abs(c[useful]), 1e-12))
        out['amplitude_rms_change_vs_phase_off_db'] = float(np.sqrt(np.average(d*d, weights=importance)))
    return out


def bend_code(cents: int) -> str:
    b64 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
    c = cents % 4096
    return b64[c // 64] + b64[c % 64]


class TestAudio(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='tn-finds-')
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        self.inp = self.path / '原声 かな.wav'
        self.out = self.path / '结果 ü.wav'

    def cmd(self, *args, ok=True):
        p = subprocess.run([str(EXE), *map(str, args)], capture_output=True, text=True,
                           encoding='utf-8', errors='replace', timeout=20)
        if ok:
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        else:
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        return p

    def native(self, *args):
        p = self.cmd('--render', self.inp, self.out, *args, '--report')
        y, sr = sf.read(self.out, dtype='float32')
        return y, sr, json.loads(p.stdout.strip().splitlines()[-1])

    def legacy(self, note, flags, length, bend='', velocity=100, consonant=120, offset=0, cutoff=0):
        self.cmd(self.inp, self.out, note, velocity, flags, offset, length, consonant, cutoff, 100, 0, '!120', bend)
        self.assertEqual(sf.info(self.out).subtype, 'PCM_16')
        y, sr = sf.read(self.out, dtype='float32')
        self.assertEqual(sr, SR)
        self.assertEqual(len(y), round(length * SR / 1000))
        return y

    def assert_tuned(self, x, hz, cents=3.):
        measured, cor = estimate_hz(x, SR)
        self.assertGreater(cor, .75, (measured, hz))
        self.assertGreater(measured, 0.)
        self.assertLess(abs(1200 * np.log2(measured / hz)), cents)

    def test_passthrough_unicode_and_negative_cutoff(self):
        x = compound()
        sf.write(self.inp, x, SR, subtype='FLOAT')
        y, sr, report = self.native()
        self.assertEqual(sr, SR)
        self.assertTrue(report['bypass'])
        np.testing.assert_array_equal(x, y)
        y, _, _ = self.native('--offset-ms', 25, '--cutoff-ms', -110)
        np.testing.assert_array_equal(x[1200:6480], y)
        y, _, _ = self.native('--offset-ms', 25, '--cutoff-ms', 110)
        np.testing.assert_array_equal(x[1200:-5280], y)

    def test_formant_only_leaves_clear_consonants_and_f0(self):
        x = compound()
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'ふぉーまんと.csv'
        track.write_text('time_s,f0_hz\n0,0\n.116,0\n.120,200\n1.097,200\n1.101,0\n1.2,0\n', encoding='utf-8')
        yr, _, report = self.native('--f0', track, '--formant-st', 3)
        yd, _, _ = self.native('--f0', track, '--formant-st', -3)
        self.assertFalse(report['bypass'])
        self.assertEqual(report['islands'], 1)
        for y in (yr, yd):
            self.assertEqual(len(y), len(x))
            np.testing.assert_array_equal(x[:int(.105*SR)], y[:int(.105*SR)])
            np.testing.assert_array_equal(x[int(1.12*SR):], y[int(1.12*SR):])
            self.assert_tuned(y[int(.45*SR):int(.8*SR)], 200.)
        m = slice(int(.4*SR), int(.6*SR))
        self.assertGreater(f2_peak(yr[m], SR), f2_peak(x[m], SR) + 170)
        self.assertLess(f2_peak(yd[m], SR), f2_peak(x[m], SR) - 140)

    def test_formant_shifts_the_mb_controlled_breath_in_both_directions(self):
        body, _ = mixed_fixture(f0=160., noise_rms=0.)
        hz = np.fft.rfftfreq(len(body), 1. / SR)
        # A broadband residual with a known broad ridge and nonzero skirts.
        profile = (.04 + .96 * np.exp(-.5*((hz-4200.)/650.)**2)) * np.clip((hz-1700.)/700., 0., 1.)
        z = np.random.default_rng(848).standard_normal(len(body))
        noise = np.fft.irfft(np.fft.rfft(z) * profile, len(body))
        t = np.arange(len(body)) / SR
        fade = np.sin(.5*np.pi*np.clip(np.minimum(t, len(body)/SR-t)/.04, 0., 1.))**2
        x = body + (.03 / np.std(noise) * noise * fade).astype(np.float32)
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'known.csv'
        points = mixed_f0_points(f0=160.)
        track.write_text('time_s,f0_hz\n' + ''.join(f'{t},{f}\n' for t, f in points), encoding='utf-8')
        for opts, samples, f0 in (((), len(x), 160.),
                                  (('--target-hz', 190., '--length-ms', 1980, '--consonant-ms', 180),
                                   round(1.98*SR), 190.)):
            with self.subTest(samples=samples):
                middle = slice(round(.52*SR), round(1.1*SR))
                band = np.fft.rfftfreq(middle.stop-middle.start, 1. / SR)
                win = np.hanning(len(band)*2-2)
                balances = []
                for step in (-3., 0., 3.):
                    args = ('--f0', track, '--formant-st', step, *opts)
                    low, _, _ = self.native(*args, '--breathiness', -100)
                    high, _, _ = self.native(*args)
                    self.assertEqual(len(low), samples)
                    self.assertEqual(len(high), samples)
                    np.testing.assert_array_equal(low[:round(.14*SR)], x[:round(.14*SR)])
                    np.testing.assert_array_equal(high[:round(.14*SR)], x[:round(.14*SR)])
                    # Use the low harmonic band for the F0 check. The breath
                    # ridge can otherwise bias the sample-lag interpolation.
                    harmonic_band = np.fft.rfft(high[middle])
                    harmonic_band[band > 1200.] = 0.
                    self.assert_tuned(np.fft.irfft(harmonic_band, len(win)), f0, cents=5.)
                    power = np.abs(np.fft.rfft((high-low)[middle] * win))**2
                    lower = np.sum(power[(band >= 3000.) & (band < 3900.)])
                    upper = np.sum(power[(band >= 4500.) & (band < 5600.)])
                    balances.append(upper / lower)
                self.assertLess(balances[0], .60 * balances[1], balances)
                self.assertGreater(balances[2], 1.7 * balances[1], balances)

    def test_edge_filter_does_not_add_a_large_low_band_burst(self):
        # A real target_hz forces synthesis even if it equals the source F0.
        # Test several starting phases so a near-zero first sample cannot hide
        # a bad startup state. Known F0 makes this independent of estimation.
        track = self.path / 'known.csv'
        track.write_text('time_s,f0_hz\n0,0\n.1499,0\n.15,160\n.6999,160\n.70,0\n.92,0\n', encoding='utf-8')
        attacks = slice(round(.12*SR), round(.235*SR))
        for phase in (0., .65, 1.4, 2.2):
            with self.subTest(phase=phase):
                x = boundary_fixture(phase)
                sf.write(self.inp, x, SR, subtype='FLOAT')
                y, _, r = self.native('--f0', track, '--target-hz', 120,
                                      '--consonant-ms', 150, '--length-ms', 1570)
                self.assertEqual(r['islands'], 1)
                self.assertEqual(len(y), round(1.57*SR))
                # Packed complex FFT channels can differ by float roundoff.
                np.testing.assert_allclose(x[:int(.13*SR)], y[:int(.13*SR)], rtol=0., atol=1e-8)
                thump = float(np.max(np.abs(subvoice_lp(y, SR)[attacks])))
                # The two-pitch crossfade still has some subvoice energy.
                # The old endpoint seed added a 0.045..0.075 full-scale burst
                # on this fixture; do not require the edit to be a null.
                self.assertLess(thump, .028, f'low-band transient {thump:.6f}')
                self.assert_tuned(y[SR//2:SR*3//4], 120.)

    def test_pitch_cycle_boundary_does_not_add_broadband_spikes(self):
        track = self.path / 'known.csv'
        track.write_text('time_s,f0_hz\n0,0\n.1499,0\n.15,160\n.6999,160\n.70,0\n.92,0\n', encoding='utf-8')
        for dst in (120, 240):
            for phase in (0., .65, 1.4, 2.2):
                with self.subTest(dst=dst, phase=phase):
                    x = stripe_fixture(phase)
                    sf.write(self.inp, x, SR, subtype='FLOAT')
                    y, _, r = self.native('--f0', track, '--target-hz', dst,
                                          '--consonant-ms', 150, '--length-ms', 1570)
                    self.assertEqual(r['islands'], 1)
                    np.testing.assert_array_equal(x[:int(.112*SR)], y[:int(.112*SR)])
                    high = stripe_band(y)
                    for center in (.15, .15 + (.70-.15) * (1.57-.15) / (.92-.15)):
                        z = high[round((center-.038)*SR):round((center+.038)*SR)]
                        peak = float(np.max(np.abs(z)))
                        self.assertLess(peak, .0015, f'{center=:.4f} : high-band spike {peak:.7f}')
                    self.assert_tuned(y[SR//2:SR*7//8], dst)

    def test_mixed_voicing_transposes_harmonics_and_preserves_noise(self):
        h, n = mixed_fixture()
        track = self.path / 'known.csv'
        np.savetxt(track, mixed_f0_points(), delimiter=',', header='time_s,f0_hz', comments='')
        for target in (113., 197.):
            with self.subTest(target=target):
                # The five-cycle endpoint intentionally retains some high
                # noise. Test preservation with the full breath mask enabled.
                opts = ('--f0', track, '--target-hz', target, '--mask-return', 1)
                signals = []
                for x in (h+n, h-n, h):
                    sf.write(self.inp, x, SR, subtype='FLOAT')
                    y, rate, report = self.native(*opts)
                    self.assertEqual((rate, len(y)), (SR, len(x)))
                    self.assertFalse(report['bypass'])
                    signals.append(y.astype(np.float64))
                plus, minus, harmonic_only = signals
                common = (plus+minus) / 2.
                for window in MIXED_WINDOWS_S:
                    with self.subTest(window=window):
                        a, b = (round(t*SR) for t in window)
                        base = harmonic_grid(harmonic_only[a:b], 149., target)
                        mixed = harmonic_grid(common[a:b], 149., target)
                        self.assertGreater(mixed['target_over_source_db'], 15., mixed)
                        self.assertGreater(mixed['target_rms'], .025, mixed)
                        self.assertGreater(mixed['target_rms'], .708*base['target_rms'], (mixed, base))
                        self.assert_tuned(common[a:b], target, cents=5.)
                        high = paired_noise_metrics(plus, minus, n, window=window)
                        self.assertLess(abs(high['gain_db']), .5, high)
                        self.assertGreater(high['correlation'], .99, high)
                        self.assertLess(high['residual_db'], -25., high)

    def test_mixed_voicing_recovers_from_missing_csv_pitch(self):
        h, n = mixed_fixture()
        sf.write(self.inp, h+n, SR, subtype='FLOAT')
        track = self.path / 'missing.csv'
        gap = (.40, .48)
        cases = [('vowel gap', mixed_f0_points(gap=gap), ((.411, .469),)),
                 ('all zero', mixed_f0_points(all_zero=True), MIXED_WINDOWS_S)]
        for name, points, windows in cases:
            with self.subTest(name=name):
                np.savetxt(track, points, delimiter=',', header='time_s,f0_hz', comments='')
                y, _, _ = self.native('--f0', track, '--target-hz', 197., '--mask-return', 1)
                self.assertEqual(len(y), len(h))
                for window in windows:
                    with self.subTest(window=window):
                        a, b = (round(t*SR) for t in window)
                        # Keep well-separated low lines for the shorter window.
                        grid = harmonic_grid(y[a:b], 149., 197., max_hz=500.)
                        self.assertGreater(grid['target_over_source_db'], 12., grid)
                        self.assertGreater(grid['target_rms'], .025, grid)
                        # y and -y extract y directly since this check uses
                        # only one sign of the noise input.
                        high = paired_noise_metrics(y, -y, n, window=window)
                        self.assertLess(abs(high['gain_db']), .5, high)
                        self.assertGreater(high['correlation'], .99, high)
                        self.assertLess(high['residual_db'], -22., high)

    def test_supplied_voicing_does_not_synthesize_breath(self):
        _, n = mixed_fixture()
        t = np.fft.rfftfreq(len(n), 1/SR)
        spec = np.fft.rfft(np.random.default_rng(618).standard_normal(len(n)))
        spec[(t < 350) | (t > 1500)] = 0.
        breath = np.fft.irfft(spec, len(n))
        x = n + (.012 * breath / np.sqrt(np.mean(breath**2))).astype(np.float32)
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'false_positive.csv'
        track.write_text('time_s,f0_hz\n0,149\n1.6,149\n', encoding='utf-8')
        y, _, _ = self.native('--f0', track, '--target-hz', 197.)
        self.assertEqual(len(y), len(x))
        for window in MIXED_WINDOWS_S:
            with self.subTest(window=window):
                a, b = (round(t*SR) for t in window)
                change = y[a:b].astype(np.float64) - x[a:b]
                self.assertLess(np.linalg.norm(change) / np.linalg.norm(x[a:b]), .018)
                grid = harmonic_grid(change, 149., 197.)
                self.assertLess(grid['target_rms'], .00025, grid)
                high = paired_noise_metrics(y, -y, x, window=window)
                self.assertLess(abs(high['gain_db']), .5, high)
                self.assertLess(high['residual_db'], -25., high)

    def test_mixed_voicing_transposes_a_weak_fundamental(self):
        h, n = mixed_fixture(f0=83., h1_scale=.025)
        sf.write(self.inp, h+n, SR, subtype='FLOAT')
        track = self.path / 'weak.csv'
        np.savetxt(track, mixed_f0_points(f0=83.), delimiter=',', header='time_s,f0_hz', comments='')
        y, _, _ = self.native('--f0', track, '--target-hz', 109.)
        for a, b in MIXED_WINDOWS_S:
            with self.subTest(window=(a, b)):
                grid = harmonic_grid(y[round(a*SR):round(b*SR)], 83., 109., max_hz=650.)
                self.assertGreater(grid['target_over_source_db'], 12., grid)
                self.assertGreater(grid['target_rms'], .025, grid)
                # H1 is weak by design. Odd target lines still distinguish
                # a genuine pitch shift from synthesis at twice the target.
                self.assertGreater(grid['target_odd_rms'], .010, grid)

    def test_narrow_resonance_is_not_a_voice_by_itself(self):
        # A 6 kHz carrier can correlate strongly at the supplied 150 Hz lag.
        rng = np.random.default_rng(611)
        t = np.arange(SR) / SR
        span = np.fft.rfftfreq(len(t), 1/SR)
        spec = np.fft.rfft(rng.normal(size=len(t)))
        spec[(span < 5985) | (span > 6015)] = 0.
        x = np.fft.irfft(spec, len(t))
        x = (.04 * x / np.sqrt(np.mean(x**2))).astype(np.float32)
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'wrong.csv'
        track.write_text('time_s,f0_hz\n0,150\n1,150\n', encoding='utf-8')
        for amount in (0., .35, 1.):
            with self.subTest(mask_return=amount):
                y, _, report = self.native('--f0', track, '--target-hz', 197., '--mask-return', amount)
                self.assertEqual(report['islands'], 0)
                np.testing.assert_array_equal(x, y)

    def test_mask_keeps_a_supplied_fundamental_above_the_crossover(self):
        f0, target = 2310., 1510.
        t = np.arange(SR) / SR
        fade = np.clip(np.minimum(t, 1.-t)/.05, 0., 1.)
        x = .13*fade*np.sin(2.*np.pi*f0*t)
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'upper_f0.csv'
        track.write_text(f'time_s,f0_hz\n0,{f0}\n1,{f0}\n', encoding='utf-8')
        y, _, _ = self.native('--f0', track, '--target-hz', target, '--mask-return', 1)
        segment = slice(SR//4, 3*SR//4)
        def amplitude(hz):
            return 2.*abs(np.dot(y[segment], np.exp(-2j*np.pi*hz*t[segment])))/len(y[segment])
        self.assertGreater(amplitude(target), .06)
        self.assertLess(amplitude(f0), .007)

    def test_octave_up_does_not_introduce_dc(self):
        t = np.arange(SR) / SR
        env = np.clip(np.minimum(t-.16, .84-t)/.022, 0., 1.)
        x = (.20 * np.cos(2*np.pi*160*t+.4) * np.sin(np.pi/2*env)**2).astype(np.float32)
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'octave.csv'
        track.write_text('time_s,f0_hz\n0,160\n1,160\n', encoding='utf-8')
        y, _, _ = self.native('--f0', track, '--target-hz', 320.)
        mid = y[SR//4:SR*3//4]
        self.assert_tuned(mid, 320.)
        self.assertGreater(np.std(mid), .025)
        self.assertLess(abs(np.mean(mid)), .001)
        self.assertLess(np.max(np.abs(subvoice_lp(y, SR)[round(.11*SR):round(.89*SR)])), .026)

    def test_pitch_change_rejects_the_original_harmonic_grid(self):
        sf.write(self.inp, vowel(200., 1.2), SR, subtype='FLOAT')
        track = self.path / 'known.csv'
        track.write_text('time_s,f0_hz\n0,200\n1.2,200\n', encoding='utf-8')
        y, _, _ = self.native('--f0', track, '--target-hz', 220)
        z = y[SR//4:3*SR//4].astype(np.float64)
        hz = np.fft.rfftfreq(len(z), 1/SR)
        a = np.abs(np.fft.rfft(z))**2
        old_only = (hz % 200 == 0) & (hz % 220 != 0)
        self.assertLess(np.sum(a[old_only]) / np.sum(a[1:-1]), 10**(-46/10))
        self.assert_tuned(z, 220.)

    def test_transposing_a_single_partial_keeps_the_target_fundamental(self):
        sr = SR
        t = np.arange(2*sr) / sr
        fade = np.sin(.5*np.pi*np.clip(np.minimum(t, 2.-t)/.08, 0., 1.))**2
        x = .25 * fade * np.cos(2*np.pi*240*t + .37)
        sf.write(self.inp, x, sr, subtype='FLOAT')
        track = self.path / 'single_partial.csv'
        track.write_text('time_s,f0_hz\n0,240\n2,240\n', encoding='ascii')
        for target in (100, 120, 160, 360, 480):
            with self.subTest(target=target):
                y, rate, _ = self.native('--f0', track, '--target-hz', target)
                self.assertEqual((len(y), rate), (len(x), sr))
                mid = y[sr//2:3*sr//2]
                powers = np.abs(np.fft.rfft(mid))**2
                fundamental = powers[target]
                powers[target] = 0.
                self.assertGreater(2*np.sqrt(fundamental)/len(mid), .22)
                self.assertLess(powers.sum()/fundamental, 10**(-38/10))

    def test_resynthesis_does_not_comb_filter_the_envelope(self):
        t = np.arange(2*SR) / SR
        freq = 240*np.arange(1, 18)
        amplitude = np.clip((4100-freq)/1000., 0., 1.)
        phases = np.random.default_rng(66).uniform(-np.pi, np.pi, len(freq))
        x = sum(a*np.cos(2*np.pi*f*t + p) for f, a, p in zip(freq, amplitude, phases))
        fade = np.sin(.5*np.pi*np.clip(np.minimum(t, 2.-t)/.08, 0., 1.))**2
        sf.write(self.inp, .45*fade*x/np.max(np.abs(x)), SR, subtype='FLOAT')
        track = self.path / 'comb.csv'
        track.write_text('time_s,f0_hz\n0,240\n2,240\n', encoding='ascii')
        for target in (151, 173, 343):
            with self.subTest(target=target):
                y, _, _ = self.native('--f0', track, '--target-hz', target)
                lines = target*np.arange(1, 3000//target+1)
                measured = np.abs(np.fft.rfft(y[SR//2:3*SR//2])[lines])
                self.assertGreater(np.median(measured), 10.)
                error_db = 20*np.log10(np.maximum(measured, 1e-9)/np.median(measured))
                self.assertLess(np.sqrt(np.mean(error_db**2)), 1.2, error_db)
                self.assertGreater(np.min(error_db), -3., error_db)

    def test_stretched_upshift_preserves_envelope_without_cycle_rate_sidebands(self):
        src, target, duration, length = 240, 360, .48, 3.
        t = np.arange(round(duration*SR)) / SR
        phi = 2*np.pi*src*t
        mod_hz = 28.
        x = (.095*np.cos(phi+.31) + .08*np.cos(2*phi-.61)
             + .09*(1+.28*np.cos(2*np.pi*mod_hz*t))*np.cos(3*phi+.74)
             + .06*np.cos(4*phi-1.3) + .04*np.cos(5*phi+.4) + .035*np.cos(6*phi-1.7))
        fade = np.sin(.5*np.pi*np.clip(np.minimum(t, duration-t)/.04, 0., 1.))**2
        sf.write(self.inp, x*fade, SR, subtype='FLOAT')
        y, rate, _ = self.native('--forced-hz', src, '--target-hz', target,
                                  '--length-ms', round(1000*length), '--breathiness', -100, '--mask-return', 0)
        self.assertEqual((len(y), rate), (round(length*SR), SR))

        # Source H3 and target H2 share 720 Hz. Source AM should appear at
        # +/-4.48 Hz from that carrier after the stretch. Selecting the
        # closest single cycle creates extra +/-33.92 and +/-42.88 Hz images.
        core = y[3*SR//4:9*SR//4]
        win = np.hanning(len(core))
        freq = np.fft.rfftfreq(8*SR, 1/SR)
        power = np.abs(np.fft.rfft(core*win, n=8*SR))**2
        gap = np.abs(freq-2*target)
        carrier = power[gap<1.2].sum()
        self.assertGreater(2*np.sqrt(power[gap==0].sum())/win.sum(), .04)
        am = power[np.abs(gap-mod_hz*duration/length)<2.].sum()
        spurs = power[(gap>20.) & (gap<145.)].sum()
        self.assertGreater(am/carrier, 10**(-23./10), 'lost original amplitude movement')
        self.assertLess(spurs/carrier, 10**(-40./10),
                        f'cycle-rate sidebands: {10*np.log10(max(spurs, 1e-30)/carrier):.1f} dBc')

    def test_upshift_wave_table_does_not_fold_audio_band_images(self):
        t = np.arange(2*SR) / SR
        source, target = 240, 353
        orders = np.arange(1, 76)
        f = source*orders
        a = (.30*np.exp(-.5*((f-12600)/2800)**2)
             + .16*np.exp(-.5*((f-850)/700)**2) + .05*np.exp(-f/7000))
        phases = np.random.default_rng(66).uniform(-np.pi, np.pi, len(orders))
        x = sum(amp*np.cos(2*np.pi*hz*t+p) for amp, hz, p in zip(a, f, phases))
        fade = np.sin(.5*np.pi*np.clip(np.minimum(t, 2-t)/.08, 0., 1.))**2
        sf.write(self.inp, .48*fade*x/np.max(np.abs(x)), SR, subtype='FLOAT')
        track = self.path / 'bright.csv'
        track.write_text('time_s,f0_hz\n0,240\n2,240\n', encoding='ascii')
        y, _, _ = self.native('--f0', track, '--target-hz', target, '--breathiness', -100)
        power = np.abs(np.fft.rfft(y[SR//2:3*SR//2]))**2
        strongest = np.max(power[target::target])
        # The 256-point cubic table leaked (256+k)*target, folded down
        # around 4.6–5.7 kHz. These bins are separate from both grids.
        image_bins = (256+np.arange(29, 33))*target-2*SR
        self.assertLess(np.sum(power[image_bins])/strongest, 10**(-63./10))

    def assert_smooth_partial(self, y, hz):
        t = np.arange(len(y)) / SR
        # Demodulate one partial without averaging over the complete note;
        # a brief 180-degree flip can disappear from a long FFT measurement.
        mixed = 2*y*np.exp(-2j*np.pi*hz*t)
        baseband_hz = np.fft.fftfreq(len(y), 1/SR)
        gaussian = np.exp(-.5*(2*np.pi*.003*baseband_hz)**2)
        partial = np.fft.ifft(np.fft.fft(mixed)*gaussian)[SR//2:3*SR//2]
        step = round(.006*SR)
        delta = np.abs(np.angle(partial[step:] * np.conj(partial[:-step])))
        amplitude = np.abs(partial)
        self.assertGreater(np.median(amplitude), .045)
        self.assertGreater(np.min(amplitude)/np.median(amplitude), .80)
        self.assertLess(np.max(delta), .22, f'{hz} Hz phase step {np.rad2deg(np.max(delta)):.1f} degrees')

    def test_phase_arc_crossing_does_not_flip_the_target_partial(self):
        t = np.arange(2*SR) / SR
        phase = 2*np.pi*240*t
        fade = np.sin(.5*np.pi*np.clip(np.minimum(t, 2-t)/.08, 0., 1.))**2
        # The difference of source H2 and H3 keeps crossing the +/-pi cut.
        # Both shifts place a target partial halfway between those harmonics.
        relative = np.pi + .08*np.sin(2*np.pi*10*(t-1.003))
        x = .1*fade*(np.cos(phase) + np.cos(2*phase) + np.cos(3*phase+relative))
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'crossing.csv'
        track.write_text('time_s,f0_hz\n0,240\n2,240\n', encoding='ascii')
        for target in (200, 300):
            with self.subTest(target=target):
                y, _, _ = self.native('--f0', track, '--target-hz', target, '--breathiness', -100)
                self.assert_smooth_partial(y, 600)

    def test_weak_fundamental_does_not_flip_stronger_harmonics(self):
        t = np.arange(2*SR) / SR
        phase = 2*np.pi*240*t
        fade = np.sin(.5*np.pi*np.clip(np.minimum(t, 2-t)/.08, 0., 1.))**2
        weak_h1 = .16*(1.-np.exp(-((t-1.)/.12)**2))
        x = fade*(weak_h1*np.cos(phase+.84) + .16*np.cos(2*phase+.7) + .16*np.cos(3*phase-.8))
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'weak-h1.csv'
        track.write_text('time_s,f0_hz\n0,240\n2,240\n', encoding='ascii')
        y, _, _ = self.native('--f0', track, '--target-hz', 150, '--breathiness', -100)
        self.assert_smooth_partial(y, 600)

    def test_nearly_sinusoidal_source_does_not_spawn_strong_overtones(self):
        t = np.arange(2*SR) / SR
        phase = 2*np.pi*240*t
        fade = np.sin(.5*np.pi*np.clip(np.minimum(t, 2-t)/.08, 0., 1.))**2
        # H2 slowly exceeds the old pure-tone threshold. That switch used
        # to copy source H1 onto multiple low-order target harmonics.
        x = fade*(.25*np.cos(phase) + (.0035355+.0025*np.tanh(6*(t-1.))) * np.cos(2*phase))
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'nearly-sinusoidal.csv'
        track.write_text('time_s,f0_hz\n0,240\n2,240\n', encoding='ascii')
        y, _, _ = self.native('--f0', track, '--target-hz', 80, '--breathiness', -100)
        window = np.hanning(SR//4)
        for start in range(SR//2, 5*SR//4, SR//16):
            amplitudes = 2*np.abs(np.fft.rfft(y[start:start+len(window)]*window))/window.sum()
            self.assertGreater(amplitudes[80//4], .20)
            self.assertLess(max(amplitudes[160//4], amplitudes[240//4]), .025)

    def test_phase_formant_preserves_steady_harmonic_phase(self):
        x = vowel(200., 1.2)
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'known.csv'
        track.write_text('time_s,f0_hz\n0,200\n1.2,200\n', encoding='utf-8')
        y, _, _ = self.native('--f0', track, '--formant-st', 3)
        mid = slice(SR//4, 3*SR//4)
        a = np.fft.rfft(x[mid])[100:2000:100]
        b = np.fft.rfft(y[mid])[100:2000:100]
        energy = abs(a)**2
        delta = np.angle(b * np.conj(a))
        self.assertLess(np.sqrt(np.average(delta**2, weights=energy)), .045)
        self.assert_tuned(y[mid], 200.)

    def test_phase_lock_corrects_relative_phase_without_changing_level(self):
        x, src, dst = phase_fixture()
        sf.write(self.inp, x, SR, subtype='FLOAT')
        track = self.path / 'known.csv'
        track.write_text(f'time_s,f0_hz\n0,{src}\n1.6,{src}\n', encoding='utf-8')
        off, _, _ = self.native('--target-hz', dst, '--f0', track)
        on, _, _ = self.native('--target-hz', dst, '--f0', track, '--phase-lock')
        before = relative_shape_metrics(x, off, src, dst)
        after = relative_shape_metrics(x, on, src, dst)
        self.assertEqual(before['partial_count'], after['partial_count'])
        self.assertLess(after['relative_phase_rms_degrees'], before['relative_phase_rms_degrees'] * .55)

        # The default renderer can now change magnitude. Check the legacy
        # phase-lock path against its log-envelope interpolation directly.
        orders = np.arange(1, 49)
        mid = slice(SR//4, SR//4+SR)
        source_amp = np.abs(np.fft.rfft(x[mid])[orders*src])
        locked_amp = np.abs(np.fft.rfft(on[mid])[orders*dst])
        envelope = np.interp(orders*dst, orders*src, np.log(source_amp))
        relative_db = (20./np.log(10.))*(np.log(locked_amp)-envelope)
        reliable = (orders > 1) & (source_amp > .035*source_amp.max()) & (locked_amp > .02*locked_amp.max())
        self.assertLess(np.sqrt(np.mean((relative_db[reliable]-relative_db[0])**2)), .15)
        self.assert_tuned(on[SR//4:SR//4+SR], dst, cents=4.)

    def test_utau_pitch_bend_flags_and_length(self):
        sf.write(self.inp, compound(), SR, subtype='PCM_16')
        bends = bend_code(-100) + '#126#' + bend_code(100) + '#234#'
        y = self.legacy('A3', 'eO2g-7.5', 1510, bends)
        self.assert_tuned(y[int(.22*SR):int(.50*SR)], 220 * 2**(-1/12))
        self.assert_tuned(y[int(.86*SR):int(1.20*SR)], 220 * 2**(1/12))
        y = self.legacy('Bb3', 'et-10', 850)
        self.assert_tuned(y[int(.30*SR):int(.70*SR)], 220.)
        self.assertLessEqual(float(np.max(np.abs(y))), .501)

    def test_fast_bend_matches_the_integrated_pitch_curve(self):
        t = np.arange(2*SR)/SR
        sf.write(self.inp, .25*np.cos(2*np.pi*162*t+.32), SR, subtype='FLOAT')
        bpm = 180
        hop = SR*60./(96*bpm)
        positions = np.arange(int(np.ceil(len(t)/hop))+1)*hop
        bends = np.rint(560*np.sin(2*np.pi*17*positions/SR)).astype(int)
        self.cmd(self.inp, self.out, 'A2', 100, 'W162eMb-100', 0, 2000, 0, 0,
                 100, 0, f'!{bpm}', ''.join(bend_code(x) for x in bends))
        y, _ = sf.read(self.out)
        knots = np.arange(0, len(t)+1, round(.001*SR))
        hz = 110.*np.exp2(np.interp(knots, positions, bends)/1200.)
        hz = np.interp(np.arange(len(t)+1), knots, hz)
        phase = 2*np.pi/SR*np.r_[0., np.cumsum((hz[:-1]+hz[1:])*.5)][:len(t)]
        mid = slice(SR//4, -SR//4)
        basis = np.c_[np.cos(phase[mid]), np.sin(phase[mid])]
        ideal = basis@np.linalg.lstsq(basis, y[mid], rcond=None)[0]
        error = np.linalg.norm(y[mid]-ideal)/np.linalg.norm(y[mid])
        self.assertLess(error, 10**(-70/20))

    def test_pitch_and_stretch_preserve_fixed_consonant_and_use_sustain(self):
        x = compound()
        sf.write(self.inp, x, SR, subtype='FLOAT')
        y, sr, report = self.native('--pitch-st', 5, '--length-ms', 2050,
                                    '--consonant-ms', 100, '--loop')
        self.assertEqual((sr, len(y)), (SR, 98400))
        self.assertTrue(report['loop_used'])
        np.testing.assert_allclose(x[:int(.1*SR)], y[:int(.1*SR)], rtol=0., atol=1e-8)
        self.assert_tuned(y[int(.9*SR):int(1.5*SR)], 200 * 2**(5/12), cents=6.)
        self.assertTrue(np.isfinite(y).all())

    def test_unvoiced_only_and_silent_render(self):
        x = np.random.default_rng(1).normal(0, .02, size=SR // 2).astype(np.float32)
        sf.write(self.inp, x, SR, subtype='FLOAT')
        y, _, report = self.native('--pitch-st', 5, '--formant-st', -5, '--breathiness', 100)
        self.assertEqual(report['islands'], 0)
        np.testing.assert_array_equal(x, y)
        sf.write(self.inp, np.zeros(1000, dtype=np.float32), SR, subtype='PCM_16')
        self.assertEqual(np.count_nonzero(self.legacy('C4', 'Mb100W-1', 62, consonant=0)), 0)

    def test_all_fixed_input_pads_instead_of_repeating_last_sample(self):
        x = (.15 + .1 * np.cos(np.arange(4800) / 13)).astype(np.float32)
        sf.write(self.inp, x, SR, subtype='FLOAT')
        y, _, _ = self.native('--length-ms', 240, '--consonant-ms', 100, '--loop')
        self.assertEqual(len(y), 11520)
        np.testing.assert_array_equal(x[:4600], y[:4600])
        self.assertEqual(np.count_nonzero(y[4800:]), 0)

    def test_negative_mb_is_linear_and_preserves_clear_noise(self):
        x = compound()
        sf.write(self.inp, x, SR, subtype='FLOAT')
        low, _, _ = self.native('--breathiness', -100)
        neutral, _, _ = self.native('--breathiness', 0)
        mid, _, _ = self.native('--breathiness', -50)
        np.testing.assert_array_equal(x, neutral)
        np.testing.assert_allclose(mid, (low + neutral) / 2, rtol=0., atol=1e-7)
        for y in (low, mid):
            np.testing.assert_array_equal(x[:int(.10*SR)], y[:int(.10*SR)])
            np.testing.assert_array_equal(x[int(1.12*SR):], y[int(1.12*SR):])
            self.assert_tuned(y[SR//3:SR*3//4], 200., cents=6.)
        self.assertGreater(np.linalg.norm(low - neutral), .05)

    def test_classic_mb_and_formant_match_native_controls(self):
        sf.write(self.inp, compound(), SR, subtype='FLOAT')
        g = -25
        formant_st = -g * (12. * np.log2(10.) / 200.)
        for mb in (-100, 0, 100):
            with self.subTest(mb=mb):
                y = self.legacy('A3', f'eMb{mb}g{g}', 1200)
                ref, _, _ = self.native('--target-hz', 220, '--consonant-ms', 120,
                                         '--breathiness', mb, '--formant-st', formant_st,
                                        '--peak-normalize', '--pcm16')
                np.testing.assert_array_equal(y, ref)

    def test_stereo_pcm24_downmix_and_odd_chunk(self):
        n = 1024
        x = np.cos(np.arange(n) / 10) * .25
        sf.write(self.inp, np.column_stack((x, x * .5)), 24000, subtype='PCM_24')
        expected, sr = sf.read(self.inp, dtype='float32')
        initial = self.inp.read_bytes()
        chunk = b'JUNK\x03\x00\x00\x00abc\x00'
        data = initial[:4] + (len(initial) + len(chunk) - 8).to_bytes(4, 'little') + initial[8:12] + chunk + initial[12:]
        self.inp.write_bytes(data)
        y, got_sr, _ = self.native()
        self.assertEqual(sr, got_sr)
        np.testing.assert_array_equal(expected.mean(axis=1), y)

    def test_bad_parameters_fail_without_creating_output(self):
        sf.write(self.inp, compound(), SR, subtype='FLOAT')
        self.cmd('--render', self.inp, self.out, '--length-ms', 'nan', ok=False)
        self.cmd('--render', self.inp, self.out, '--offset-ms', 2000, ok=False)
        self.cmd('--render', self.inp, self.inp, ok=False)
        self.cmd(self.inp, self.out, 'no-note', 100, '', 0, 400, 120, 0, 100, 0, ok=False)
        self.cmd(self.inp, self.out, 'A3', 100, '', 0, 400, 120, 0, 100, 0, '!120', 'AA#999999999#', ok=False)
        bad = self.path / 'broken.wav'
        bad.write_bytes(b'RIFF\x7f\xff\xff\xffWAVEfmt ')
        self.cmd('--render', bad, self.out, ok=False)
        track = self.path / 'bad.csv'
        track.write_text('time_s,f0_hz\n0.4,220\n0.4,225\n', encoding='utf-8')
        self.cmd('--render', self.inp, self.out, '--formant-st', 3, '--f0', track, ok=False)
        self.assertFalse(self.out.exists())
        self.assertAlmostEqual(sf.info(self.inp).duration, 1.2)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--exe', type=Path, default=EXE)
    args, unittest_args = ap.parse_known_args()
    EXE = args.exe.resolve()
    unittest.main(argv=[__file__, *unittest_args], verbosity=2)
