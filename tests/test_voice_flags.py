# SPDX-License-Identifier: GPL-3.0-only
"""Moresampler-style flag regressions (synthetic audio, NumPy and SoundFile)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
import soundfile as sf

SR, F0, TARGET, MS = 24000, 125, 220, 950
N = round(MS*SR/1000)
SUSTAIN = slice(round(.34*SR), round(.74*SR))
CONTROLS = ('--tension', '--breathiness', '--mouth', '--resonance', '--formant-emphasis')
EXE = BASELINE = None


def fixture():
    t = np.arange(N)/SR
    f = np.arange(F0, 9500, F0)
    amplitude = (.90*np.exp(-.5*((f-640.)/125.)**2)
                 + .53*np.exp(-.5*((f-2250.)/240.)**2)
                 + .31*np.exp(-.5*((f-3180.)/335.)**2)
                 + .047*np.exp(-f/6000))
    amplitude[0] += .5  # A credible source F0 despite the vowel-shaped spectrum.
    rng = np.random.default_rng(1379)
    harmonic = sum(a*np.cos(2*np.pi*k*t+phase) for k, a, phase in zip(f, amplitude, rng.uniform(-np.pi, np.pi, len(f))))
    gate = np.sin(.5*np.pi*np.clip(np.minimum(t-.15, .87-t)/.021, 0., 1.))**2
    harmonic *= gate*.57/np.max(np.abs(harmonic))
    noise = np.fft.rfft(rng.normal(size=N))
    bins = np.fft.rfftfreq(N, 1/SR)
    noise[(bins < 5300.) | (bins > 8800.)] = 0.
    noise = np.fft.irfft(noise, N)
    noise *= .019/np.sqrt(np.mean(noise**2))
    clear = .8*noise*(1.-gate)
    return (harmonic+clear).astype('float32'), (harmonic+clear+.65*gate*noise).astype('float32')


def lf_reference(rd):
    """Numerically sample one glottal derivative period instead of evaluating its spectrum."""
    ta = (-1.+4.8*rd)/100.
    rk = .224+.118*rd
    rg = .25*rk/(.11*rd/(.5+1.2*rk)-ta)
    te = (1.+rk)/(2.*rg)
    wg = 2.*np.pi*rg
    timing = np.arange(65536)/65536.
    t = timing[timing < te]
    u = timing[timing >= te]-te

    def solve(lo, hi, positive):
        for _ in range(42):
            mid = (lo+hi)/2.
            if positive(mid):
                lo = mid
            else:
                hi = mid
        return (lo+hi)/2.

    epsilon = solve(0., 1./ta, lambda e: -np.expm1(-e*(1.-te)) > e*ta)
    tail = -(np.exp(-epsilon*u)-np.exp(-epsilon*(1.-te)))/(epsilon*ta)
    shape = -np.sin(wg*t)/np.sin(wg*te)
    alpha = solve(0., 256., lambda a: np.sum(shape*np.exp(a*(t-te))) > -tail.sum())
    return np.fft.rfft(np.r_[shape*np.exp(alpha*(t-te)), tail])


def harmonics(x, f0=F0, high=9000):
    x = x[SUSTAIN]
    f = np.arange(f0, high+1, f0)
    w = np.hanning(len(x))
    z = np.fft.rfft(x*w)[np.rint(f*len(x)/SR).astype(int)]/w.sum()
    return f, z


def db_ratio(a, b):
    return float(10.*np.log10(max(a, 1e-30)/max(b, 1e-30)))


def harmonic_share(x, f0):
    y = x[SUSTAIN].astype('float64')
    p = (np.abs(np.fft.rfft(y))/len(y))**2
    orders = np.arange(f0, SR/2, f0)
    at = np.rint(orders*len(y)/SR).astype(int)
    return float(2*p[at].sum()/max(np.mean(y*y), 1e-30))


class TestVoiceFlags(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if EXE is None:
            raise unittest.SkipTest('select a renderer with --exe')
        temp = tempfile.TemporaryDirectory(prefix='tnf-voice-flags-')
        cls.addClassCleanup(temp.cleanup)
        sandbox = Path(temp.name)
        cls.inp, cls.mix, cls.out = [sandbox/name for name in ('vowel.wav', 'mix.wav', 'out.wav')]
        cls.x, cls.with_noise = fixture()
        sf.write(cls.inp, cls.x, SR, subtype='FLOAT')
        sf.write(cls.mix, cls.with_noise, SR, subtype='FLOAT')

    def run_cli(self, *args, exe=None, ok=True):
        self.out.unlink(missing_ok=True)
        p = subprocess.run([str(exe or EXE), *map(str, args)], capture_output=True,
                           text=True, encoding='utf-8', errors='replace', timeout=20)
        if not ok:
            self.assertNotEqual(p.returncode, 0, p.stdout)
            self.assertFalse(self.out.exists())
        else:
            self.assertEqual(p.returncode, 0, p.stderr + p.stdout)
        return p

    def native(self, *args, inp=None, exe=None, pcm=False):
        fmt = ('--peak-normalize', '--pcm16') if pcm else ()
        p = self.run_cli('--render', inp or self.inp, self.out, '--report', *fmt, *args, exe=exe)
        y, sr = sf.read(self.out, dtype='float32')
        self.assertEqual(sr, SR)
        self.assertTrue(np.isfinite(y).all())
        return y, json.loads(p.stdout)

    def voice(self, *args, inp=None):
        return self.native('--forced-hz', F0, '--breathiness', -100, *args, inp=inp)[0]

    def classic(self, flags, *, length=MS, exe=None):
        self.run_cli(self.inp, self.out, 'A3', 100, flags, 0, length, 110, 0, 100, 0, '!120', exe=exe)
        return self.out.read_bytes()

    def test_neutral_flags_stay_byte_identical(self):
        zero = tuple(v for flag in CONTROLS for v in (flag, '-0'))
        for settings in ((), ('--forced-hz', F0, '--target-hz', TARGET, '--length-ms', 1130),
                         ('--forced-hz', F0, '--target-hz', 94, '--formant-st', 2.6),
                         ('--forced-hz', F0, '--target-hz', 155, '--length-ms', 1300, '--loop',
                          '--formant-st', -1.6, '--mask-return', .65, '--pitchbend', 'AAAbAZAMAA#19#')):
            with self.subTest(settings=settings):
                _, report = self.native(*settings, inp=self.mix)
                current = self.out.read_bytes()
                self.assertEqual(report['bypass'], not settings)
                self.native(*settings, *zero, inp=self.mix)
                self.assertEqual(self.out.read_bytes(), current)
                if BASELINE:
                    self.native(*settings, inp=self.mix, exe=BASELINE)
                    self.assertEqual(self.out.read_bytes(), current)
        current = self.classic('eW125t-7')
        self.assertEqual(current, self.classic('eMt0Mb0Mo0Mr0ME0W125t-7'))
        self.assertEqual(current, self.classic('eMd65MC40MG50MD65Ms7Mm0W125t-7'))
        if BASELINE:
            self.assertEqual(current, self.classic('eW125t-7', exe=BASELINE))

    def test_removed_b_flag_is_ignored_but_lowercase_b_still_works(self):
        settings = 'eW125Mb30g-20'
        plain = self.classic(settings)
        for flag in ('B0', 'B50', 'B100'):
            with self.subTest(flag=flag):
                self.assertEqual(plain, self.classic(f'e{flag}W125Mb30g-20'))
                self.assertEqual(plain, self.classic(settings + flag))
        self.assertNotEqual(plain, self.classic(settings + 'b40'))

    def test_removed_native_breath_option_is_rejected(self):
        result = self.run_cli('--render', self.inp, self.out, '--breath', 50, ok=False)
        self.assertIn('unknown option --breath', result.stderr)

    def test_classic_compound_tokens_match_native_and_clamp(self):
        common = ('--forced-hz', F0, '--target-hz', TARGET, '--consonant-ms', 110,
                  '--length-ms', MS, '--consonant-boost', 40, '--mask-return', .45)
        self.native(*common, '--tension', 40, '--breathiness', -27.5, '--mouth', 60,
                    '--resonance', 28, '--formant-emphasis', -61, pcm=True)
        self.assertEqual(self.out.read_bytes(), self.classic('eMt+40Mb-27.5Mo60Mr28ME-61W125b40Tnr45'))
        self.assertEqual(self.classic('eMt100Mb-100Mo100Mr-100ME100W125'),
                         self.classic('eMt109Mb-109Mo109Mr-109ME109W125'))
        for flag in CONTROLS:
            for value in ('100.01', '-100.01', 'nan'):
                with self.subTest(flag=flag, value=value):
                    self.run_cli('--render', self.inp, self.out, flag, value, ok=False)
        self.run_cli('--render', self.inp, self.out, '--mouth', ok=False)

    def test_tension_changes_slope_and_ring_is_local(self):
        f, z = harmonics(self.voice())
        power = abs(z)**2
        energy = lambda p, lo, hi: p[(f >= lo) & (f <= hi)].sum()
        for sign in (-1, 1):
            with self.subTest(sign=sign):
                _, mt = harmonics(self.voice('--tension', sign*100))
                _, mr = harmonics(self.voice('--resonance', sign*100))
                tilt = lambda p: energy(p, 2300, 5000)/energy(p, 200, 900)
                self.assertGreater(sign*db_ratio(tilt(abs(mt)**2), tilt(power)), 4.)
                center = db_ratio(energy(abs(mr)**2, 2800, 3200), energy(power, 2800, 3200))
                low = db_ratio(energy(abs(mr)**2, 200, 1100), energy(power, 200, 1100))
                self.assertGreater(sign*center, 5.)
                self.assertLess(abs(low), .5)
                # Stable real gains must not rotate the strong partials.
                strong = power > power.max()*.012
                self.assertLess(np.max(np.abs(np.angle((mt*np.conj(z))[strong], deg=True))), 5.)
                self.assertLess(np.max(np.abs(np.angle((mr*np.conj(z))[strong], deg=True))), 5.)

    def test_tension_matches_lf_magnitude_ratios(self):
        rd = .70
        source = lf_reference(rd)
        order = np.arange(1, int(.47*SR/F0))
        f = F0*order
        rr = 128./(9.*np.pi**2)
        lr = 8.*.015/(3.*np.pi*340.)
        # Flow derivative at the lips: integrate once, then apply radiation.
        coeff = source[order]/np.sqrt(1.+(2.*np.pi*f*lr/rr)**2)
        y = np.real(np.exp(2j*np.pi*(np.arange(N)/SR)[:, None]*f)@coeff)
        inp = self.inp.with_name('lf-vowel.wav')
        sf.write(inp, .50*y/np.max(np.abs(y)), SR, subtype='FLOAT')
        reference_rd_by_flag = {
            50: .34+.04*(1.32458-1.25649)/(1.32908-1.25649),
            -50: 1.22+.04*(2.32458-2.29938)/(2.32964-2.29938),
        }
        for target in (95, F0, TARGET):
            args = () if target == F0 else ('--target-hz', target)
            _, base = harmonics(self.voice(*args, inp=inp), target, high=1700)
            for flag, changed in reference_rd_by_flag.items():
                with self.subTest(target=target, flag=flag):
                    _, edited = harmonics(self.voice(*args, '--tension', flag, inp=inp), target, high=1700)
                    ratio = edited/base
                    want = np.abs(lf_reference(changed)[1:len(base)+1]/source[1:len(base)+1])
                    observed_db = 20.*np.log10(np.abs(ratio[:6]/ratio[0]))
                    np.testing.assert_allclose(observed_db, 20.*np.log10(want[:6]/want[0]), atol=.8, rtol=0.)
                    self.assertLess(np.max(np.abs(np.angle(ratio[:6], deg=True))), 4.)

    def test_pitch_shift_rebuilds_lf_magnitude_on_target_harmonics(self):
        lf = lf_reference(.8)
        orders = np.arange(1, 81)
        lip = (2.*np.pi*8.*.015/(3.*np.pi*340.))/(128./(9.*np.pi**2))
        f = F0*orders
        source = lf[orders]/np.sqrt(1.+(lip*f)**2)
        t = np.arange(N)/SR
        x = np.real(np.exp(2j*np.pi*t[:, None]*f) @ source)
        inp = self.inp.with_name('lf-shift.wav')
        sf.write(inp, .50*x/np.max(np.abs(x)), SR, subtype='FLOAT')

        for target in (100, 175):
            with self.subTest(target=target):
                y, _ = self.native('--forced-hz', F0, '--target-hz', target,
                                   '--breathiness', -100, inp=inp)
                f, actual = harmonics(y, target, high=8*target)
                expected = np.abs(lf[1:9])/np.sqrt(1.+(lip*f)**2)
                ratio = np.abs(actual)/expected
                np.testing.assert_allclose(20.*np.log10(ratio/ratio[0]), 0., atol=.4, rtol=0.)

    def test_mouth_matches_target_h1_response(self):
        fs, bandwidth = 49226., 153.
        radius = np.exp(-2.*np.pi*bandwidth/fs)
        for target in (95, F0, TARGET):
            args = () if target == F0 else ('--target-hz', target)
            f, z = harmonics(self.voice(*args), target, high=5500)
            strong = (f > target) & (np.abs(z) > np.abs(z).max()*.065)
            self.assertGreater(strong.sum(), 5)
            d = 1.+4.*radius*np.sin(np.pi*(f-target)/fs)**2/(1.-radius)**2
            for value in (-100, -25, 50, 100):
                with self.subTest(target=target, value=value):
                    _, edited = harmonics(self.voice(*args, '--mouth', value), target, high=5500)
                    observed = edited/z
                    delta = 20.*np.log10(np.abs(observed/observed[0]))
                    expected = (20.*value/300.)*np.log10(d)
                    np.testing.assert_allclose(delta[strong], expected[strong], atol=.30, rtol=0.)
                    self.assertLess(np.max(np.abs(np.angle(observed[strong], deg=True))), 3.)

    def test_emphasis_changes_existing_peak_contrast(self):
        f, z = harmonics(self.voice())
        peaks = ((f >= 500) & (f <= 750)) | ((f >= 2125) & (f <= 2525))
        valleys = (f >= 1100) & (f <= 1525)
        contrast = lambda p: p[peaks].mean()/p[valleys].mean()
        power = abs(z)**2
        for sign in (-1, 1):
            with self.subTest(sign=sign):
                me = abs(harmonics(self.voice('--formant-emphasis', sign*100))[1])**2
                self.assertGreater(sign*db_ratio(contrast(me), contrast(power)), 3.)

    def test_envelope_controls_keep_the_up_and_down_shifted_harmonic_grid(self):
        for target in (95, TARGET):
            with self.subTest(target=target):
                args = ('--target-hz', target)
                baseline = self.voice(*args)
                colored = self.voice(*args, '--mouth', 55, '--resonance', -40,
                                     '--formant-emphasis', 40, '--tension', 25)
                self.assertGreater(harmonic_share(colored, target), .9)
                _, original = harmonics(baseline, target, high=4900)
                _, edited = harmonics(colored, target, high=4900)
                strong = abs(original) > .12*max(abs(original))
                phase_deg = np.abs(np.angle((edited*np.conj(original))[strong], deg=True))
                self.assertLess(np.max(phase_deg), 8.)
                self.assertGreater(np.std((colored-baseline)[SUSTAIN]), .004)

    def test_breathiness_uses_the_recorded_residual_only(self):
        # Recover the recorded residual from the difference between neutral
        # and muted breath. Check positive Mb against the same waveform.
        for inp, settings in ((self.inp, ('--target-hz', TARGET)),
                              (self.mix, ('--target-hz', 94, '--formant-st', 2.6, '--mask-return', .65))):
            with self.subTest(input=inp.name):
                def render(mb, exe=None):
                    return self.native('--forced-hz', F0, *settings, '--breathiness', mb, inp=inp, exe=exe)[0]

                baseline = render(0)
                h = render(-100)[SUSTAIN].astype('float64')
                r = baseline[SUSTAIN].astype('float64') - h
                self.assertGreater(np.std(h), .025)
                if inp == self.mix:
                    self.assertGreater(np.std(r), .004)
                for mb in (-100, -60, 25, 55, 100):
                    with self.subTest(mb=mb):
                        y = render(mb)
                        expected = (1. - max(0, mb)/100.)*h + (1. + mb/100.)*r
                        np.testing.assert_allclose(y[SUSTAIN], expected, rtol=0., atol=4.e-7)
                        np.testing.assert_array_equal(y[:round(.07*SR)], baseline[:round(.07*SR)])
                        if BASELINE:
                            current_wav = self.out.read_bytes()
                            render(mb, exe=BASELINE)
                            self.assertEqual(self.out.read_bytes(), current_wav)

    def test_breathiness_whisper_classic_matches_native(self):
        self.native('--forced-hz', F0, '--target-hz', TARGET, '--consonant-ms', 110,
                     '--length-ms', MS, '--breathiness', 100, pcm=True)
        self.assertEqual(self.out.read_bytes(), self.classic('eW125Mb100'))

    def test_clear_audio_is_not_colored(self):
        flags = ('--tension', 100, '--breathiness', 100, '--mouth', 100, '--resonance', 100,
                 '--formant-emphasis', 100)
        colored, report = self.native('--forced-hz', -1, *flags)
        self.assertEqual(report['islands'], 0)
        np.testing.assert_array_equal(colored, self.x)
        voiced, _ = self.native('--forced-hz', F0, *flags)
        np.testing.assert_array_equal(voiced[:round(.07*SR)], self.x[:round(.07*SR)])

    def test_me_prefers_safe_short_sustain_and_last_loop_flag_wins(self):
        common = ('--forced-hz', F0, '--target-hz', TARGET, '--consonant-ms', 110)
        length = MS+60
        regular, report = self.native(*common, '--length-ms', length, '--loop', pcm=True)
        self.assertFalse(report['loop_used'])
        regular_bytes = self.out.read_bytes()
        forced, report = self.native(*common, '--length-ms', length, '--force-loop', pcm=True)
        self.assertTrue(report['loop_used'])
        self.assertEqual(len(forced), round(length*SR/1000))
        self.assertNotEqual(self.out.read_bytes(), regular_bytes)
        forced_bytes = self.out.read_bytes()
        self.assertEqual(self.classic('eMeW125', length=length), forced_bytes)
        self.assertEqual(self.classic('eMeB0W125', length=length), forced_bytes)
        self.assertEqual(self.classic('eMeMs4W125', length=length), forced_bytes)
        self.assertEqual(self.classic('MeW125', length=length), forced_bytes)
        for flags in ('eW125', 'MeeW125', 'eMetal80W125', 'eMe66W125'):
            with self.subTest(flags=flags):
                self.assertEqual(self.classic(flags, length=length), regular_bytes)
        _, report = self.native(*common, '--length-ms', MS+14, '--force-loop')
        self.assertFalse(report['loop_used'])
        _, report = self.native(*common, '--forced-hz', -1, '--length-ms', length, '--force-loop')
        self.assertFalse(report['loop_used'])
        # Without peak normalization, both mappings leave the fixed attack at its original level.
        y, _ = self.native(*common, '--length-ms', length, '--force-loop')
        np.testing.assert_array_equal(y[:round(.07*SR)], self.x[:round(.07*SR)])
        self.assertGreater(np.std((regular-forced)[SUSTAIN]), .0008)

    def test_openutau_option_values_keep_e_unambiguous(self):
        # OpenUtau appends an Options expression's text verbatim, so e0/e1 and
        # Me0/Me1 must remain separate from the numeric flag before them.
        long, short = MS + 200, MS + 60
        looped = self.classic('N', length=long)
        linear = self.classic('e', length=long)
        self.assertNotEqual(looped, linear)
        self.assertEqual(self.classic('e0', length=long), looped)
        self.assertEqual(self.classic('e1', length=long), linear)
        self.assertEqual(self.classic('e50', length=long), looped)
        self.assertEqual(self.classic('Mt0e0', length=long), self.classic('Mt0', length=long))
        self.assertEqual(self.classic('Mt0e1', length=long), self.classic('Mt0e', length=long))
        self.assertNotEqual(self.classic('Mt0e1', length=long), self.classic('Mt0e0', length=long))
        # A numeric flag keeps all its digits; an option must not steal one.
        self.native('--forced-hz', F0, '--target-hz', TARGET, '--length-ms', long,
                    '--consonant-ms', 110, '--tension', 10, pcm=True)
        self.assertEqual(self.out.read_bytes(), self.classic('eMt10W125', length=long))
        self.assertNotEqual(self.classic('eMt1W125', length=long),
                            self.classic('eMt10W125', length=long))
        # Me option: 0 keeps the default mapping, 1 enables force-loop.
        self.assertEqual(self.classic('Me0', length=short), self.classic('N', length=short))
        forced = self.classic('Me1', length=short)
        self.assertEqual(forced, self.classic('Me', length=short))
        self.assertNotEqual(forced, self.classic('N', length=short))
        self.assertEqual(self.classic('Mee0', length=short), forced)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, help='optional prior executable built with the same compiler')
    args, unittest_args = parser.parse_known_args()
    EXE = args.exe.resolve(strict=True)
    BASELINE = args.baseline.resolve(strict=True) if args.baseline else None
    unittest.main(argv=[__file__, *unittest_args], verbosity=2)
