# SPDX-License-Identifier: GPL-3.0-only
"""Render/CLI regression for UTAU Tnr and native --harmonic-amount.

Run with --exe NEW. Only NumPy and SoundFile are needed. Fixtures and
outputs live in a temporary directory.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
import soundfile as sf

from test_audio import harmonic_grid
from fixtures import band_coefficients, fixture

SR, SOURCE_HZ, NOTE_HZ, PITCH_ST, LENGTH_MS = 44100, 149, 220, -.5, 1350
VOICED_SPAN = (.43, .93)
EXE = None


class TestHarmonicFlag(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if EXE is None:
            raise unittest.SkipTest('Run with --exe to select the harmonic-flag renderer')
        temp = tempfile.TemporaryDirectory(prefix='tnf-harmonic-flag-')
        cls.addClassCleanup(temp.cleanup)
        cls.inp = Path(temp.name) / 'mix.wav'
        cls.out = Path(temp.name) / 'out.wav'
        harmonic, noise = fixture(SR, SOURCE_HZ, high_rms=.035, noise_rms=.020)
        cls.x = harmonic + noise
        assert cls.x.shape == (round(SR*LENGTH_MS/1000),)
        sf.write(cls.inp, cls.x, SR, subtype='FLOAT')

    def run_cli(self, *args, exe=None, ok=True):
        self.out.unlink(missing_ok=True)
        proc = subprocess.run([str(exe or EXE), *map(str, args)], capture_output=True,
                              text=True, encoding='utf-8', errors='replace', timeout=20)
        if ok:
            self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
        else:
            self.assertNotEqual(proc.returncode, 0, proc.stdout)
            self.assertFalse(self.out.exists())
        return proc

    def read_render(self, pcm16):
        info = sf.info(self.out)
        self.assertEqual((info.samplerate, info.frames, info.channels), (SR, len(self.x), 1))
        self.assertEqual(info.subtype, 'PCM_16' if pcm16 else 'FLOAT')
        y, _ = sf.read(self.out, dtype='float32')
        self.assertTrue(np.isfinite(y).all())
        return y

    def native(self, *args, pcm16=False, exe=None):
        fmt = ('--peak-normalize', '--pcm16') if pcm16 else ()
        proc = self.run_cli('--render', self.inp, self.out, *args, *fmt, '--report', exe=exe)
        return self.read_render(pcm16), json.loads(proc.stdout.strip().splitlines()[-1])

    def pitched(self, *args, pcm16=False, exe=None):
        y, report = self.native('--forced-hz', SOURCE_HZ, '--target-hz', NOTE_HZ, '--pitch-st', PITCH_ST,
                                '--length-ms', LENGTH_MS, *args, pcm16=pcm16, exe=exe)
        self.assertFalse(report['bypass'])
        self.assertGreater(report['islands'], 0)
        return y

    def classic(self, flags):
        self.run_cli(self.inp, self.out, 'A3', 100, flags, 0, LENGTH_MS, 0, 0, 100, 0, '!120')
        return self.read_render(True)

    def test_default_is_tnr75(self):
        default = self.pitched()
        np.testing.assert_array_equal(self.pitched('--mask-return', .75), default)
        np.testing.assert_array_equal(self.pitched('--harmonic-amount', .15), self.pitched('--mask-return', .15))
        np.testing.assert_array_equal(self.classic('eW149t-5'), self.classic('eTnr75W149t-5'))

    def test_old_h_is_ignored(self):
        for h in ('H0', 'H15', 'H35', 'H100'):
            with self.subTest(h=h):
                np.testing.assert_array_equal(self.classic(f'e{h}W149t-5'), self.classic('eW149t-5'))
        np.testing.assert_array_equal(self.classic('eH100Tnr15W149t-5'), self.classic('eTnr15W149t-5'))

    def test_classic_mapping_clamping_and_audible_range(self):
        flags = (('0', 0.), ('15', .15), ('35', .35), ('67.5', .675), ('100', 1.))
        # Tnr immediately follows bare e and is followed by existing W/t.
        results = {amount: self.classic(f'eTnr{flag}W149t-5') for flag, amount in flags}
        for flag, amount in flags:
            with self.subTest(flag=flag):
                np.testing.assert_array_equal(results[amount],
                                              self.pitched('--harmonic-amount', amount, pcm16=True))
        for flag, endpoint in (('-.75', 0.), ('+100.25', 1.)):
            with self.subTest(clamp=flag):
                np.testing.assert_array_equal(self.classic(f'eTnr{flag}W149t-5'), results[endpoint])

        # A parser that accepts but ignores Tnr would also pass native/classic
        # parity. Require changes in the upper band and preserve the low F0.
        lo, hi = (round(t*SR) for t in VOICED_SPAN)
        target = NOTE_HZ * 2**(PITCH_ST/12.)
        upper = {k: band_coefficients(y, SR, VOICED_SPAN) for k, y in results.items()}
        for amount, y in results.items():
            with self.subTest(rendered_amount=amount):
                low = harmonic_grid(y[lo:hi], SOURCE_HZ, target, SR)
                self.assertGreater(low['target_over_source_db'], 20., low)
                self.assertGreater(low['target_rms'], .08, low)
        for lower, higher in zip((0., .15, .35, .675), (.15, .35, .675, 1.)):
            with self.subTest(changing_amount=(lower, higher)):
                change = np.linalg.norm(upper[higher]-upper[lower])
                level = np.linalg.norm(upper[higher])
                self.assertGreater(change, .005*level)

    def test_maximum_keeps_low_band_and_clear_frication(self):
        ref = self.pitched('--harmonic-amount', .35)
        maximum = self.pitched('--harmonic-amount', 1.)
        # Before the first voiced onset, audio is pure frication even though
        # the fixture supplies a candidate pitch throughout the file.
        stop = round(.08*SR)
        np.testing.assert_array_equal(maximum[:stop], self.x[:stop])
        # Compare before peak normalization: normalization can change the
        # amplitude of the low harmonics even if neither splitter touched them.
        low = band_coefficients(ref, SR, VOICED_SPAN, (65., 1200.))
        change = band_coefficients(maximum-ref, SR, VOICED_SPAN, (65., 1200.))
        self.assertLess(np.linalg.norm(change), .02*np.linalg.norm(low))

    def test_neutral_bypass_and_force_unvoiced(self):
        for amount in (0., 1.):
            with self.subTest(amount=amount):
                y, report = self.native('--harmonic-amount', amount)
                self.assertTrue(report['bypass'])
                np.testing.assert_array_equal(y, self.x)
        y, report = self.native('--forced-hz', -1, '--target-hz', NOTE_HZ, '--harmonic-amount', 1.)
        self.assertFalse(report['bypass'])
        self.assertEqual(report['islands'], 0)
        np.testing.assert_array_equal(y, self.x)
        normalized, _ = self.native(pcm16=True)
        np.testing.assert_array_equal(self.classic('eTnr100W-1t-5'), normalized)

    def test_native_rejects_out_of_range_or_missing_values(self):
        for suffix in (('-.0001',), ('1.0001',), ('nan',), ()):
            with self.subTest(value=suffix):
                self.run_cli('--render', self.inp, self.out, '--harmonic-amount', *suffix, ok=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    args, unittest_args = parser.parse_known_args()
    EXE = args.exe.resolve(strict=True)
    unittest.main(argv=[__file__, *unittest_args], verbosity=2)
