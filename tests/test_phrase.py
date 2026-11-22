# SPDX-License-Identifier: GPL-3.0-only
"""Check the parameter handoff, phase integral, and stock Windows batch footer."""
from __future__ import annotations

import argparse
import io
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import unittest
import wave


SR = 44100
EXE: Path | None = None
LENGTH = 550
OVERLAP = 240
ENV_A = (0, 0, OVERLAP, 100, 100, 100, 0, 0, 0, 0, 100)
ENV_B = (0, OVERLAP, 0, 0, 100, 100, 100, OVERLAP, 0, 0, 100)


def frames(ms):
    return round(ms*SR/1000)


def read_wave(path):
    with wave.open(io.BytesIO(path.read_bytes()), 'rb') as reader:
        assert (reader.getframerate(), reader.getnchannels(), reader.getsampwidth()) == (SR, 1, 2)
        count = reader.getnframes()
        return struct.unpack(f'<{count}h', reader.readframes(count))


def source(path, hz, phase):
    samples = [round(10000*math.sin(2*math.pi*i*hz/SR+phase)) for i in range(frames(650))]
    with wave.open(str(path), 'wb') as out:
        out.setparams((1, 2, SR, 0, 'NONE', ''))
        out.writeframes(struct.pack(f'<{len(samples)}h', *samples))


def model(path, hz, length=LENGTH, amplitude=.32, residual=0., overtone=None, support=None):
    """Exact analytic frame coefficients with canonical cosine H1 phase."""
    n = frames(length)
    hop, pitch_hop = frames(6), frames(1)
    count = 1+(n+hop-1)//hop if hz else 0
    pitch_count = 1+(n+pitch_hop-1)//pitch_hop if hz else 0
    header = struct.pack('<8s4i2I', b'TNFMOD1\0', SR, n, hop if hz else 0,
                         pitch_hop if hz else 0, n, count)
    data = bytearray(header + struct.pack('<f', hz)*pitch_count + struct.pack('<f', residual)*n)
    for j in range(count):
        orders = [complex(amplitude/2)]
        if overtone is not None:
            orders.append(overtone(min(n, j*hop)))
        weight = support(j) if support is not None else 1.
        data.extend(struct.pack('<ffI', hz, weight, len(orders)))
        for coefficient in orders:
            data.extend(struct.pack('<ff', weight*coefficient.real, weight*coefficient.imag))
    path.write_bytes(data)


def peak_powers(samples, frequencies):
    size = len(samples)
    windowed = [x*(.5-.5*math.cos(2*math.pi*i/(size-1))) for i,x in enumerate(samples)]
    return [abs(sum(x*complex(math.cos(w*i), math.sin(w*i)) for i,x in enumerate(windowed)))**2
            for w in (2*math.pi*hz/SR for hz in frequencies)]


class PhraseTests(unittest.TestCase):
    def setUp(self):
        if EXE is None:
            self.skipTest('Select a tn_finds executable with --exe')
        tmp = tempfile.TemporaryDirectory(prefix='tnf phrase ')
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.out = self.root/'phrase.wav'

    def call(self, *args, ok=True):
        result = subprocess.run([str(EXE), *map(str,args)], cwd=self.root, capture_output=True, timeout=20)
        if ok:
            self.assertEqual(result.returncode, 0, result.stderr.decode(errors='replace'))
        else:
            self.assertNotEqual(result.returncode, 0)

    def prepare(self, name, note, pitch):
        self.call('--prepare-note', name, note, pitch, 100, 'e', 0, 650, 0, 0, 60, 0, '!120')
        self.assertEqual((self.root/note).read_bytes()[:8], b'TNFMOD1\0')

    def test_cross_pitch_integrates_one_clock_and_snapshots_reused_cache(self):
        cache = self.root/'reused note.wav'
        model(cache, 220.)
        self.call('--phrase', self.out, cache, 0, '480@120+50', *ENV_A)
        self.assertEqual(read_wave(self.out), ())
        self.assertFalse((self.root/'phrase.wav.dat').exists())
        self.assertTrue((self.root/'phrase.wav.tnphrase').exists())
        hz = 293.6647679
        model(cache, hz)
        self.call('--phrase', self.out, cache, 0, LENGTH, *ENV_B, 'LAST_NOTE')
        self.assertFalse((self.root/'phrase.wav.tnphrase').exists())
        actual = read_wave(self.out)
        self.assertEqual(len(actual), frames(860))
        a, b = frames(310), frames(550)
        phase = 0.
        max_error = 0
        # Skip only the final spectrum-to-silence fade of the last grid interval.
        for i, sample in enumerate(actual[:-frames(4)]):
            ideal = round(32768*.32*math.cos(2*math.pi*phase))
            max_error = max(max_error, abs(sample-ideal))
            u = max(0., min(1., (i+.5-a)/(b-a)))
            phase += (220.+(hz-220.)*u)/SR
            phase %= 1.
        self.assertLessEqual(max_error, 3, f'largest phase-locked PCM error {max_error} LSB')

    def test_overtone_rotates_through_antiphase_inside_one_note(self):
        note = self.root/'phase-turn.wav'
        hz = 250.
        hop = frames(6)
        start = 30*hop
        angle = lambda n: math.pi*max(0., min(1., (n-start)/hop))
        model(note, hz, overtone=lambda n: .13*complex(math.cos(angle(n)), math.sin(angle(n))))
        self.call('--phrase', self.out, note, 0, LENGTH, 'LAST_NOTE')
        actual = read_wave(self.out)
        # Measure the middle of the rotation after subtracting known H1.
        # H2 should move around the unit circle at full amplitude.
        center = start+hop/2
        self.assert_overtone(actual, hz, round(center), round(hop*.29), angle)

    def test_overtone_overlap_keeps_same_arc_across_pi(self):
        left, right = self.root/'left.wav', self.root/'right.wav'
        hz = 250.
        begin, span = frames(LENGTH-OVERLAP), frames(OVERLAP)
        model(left, hz, overtone=lambda n: .13+0j)
        turn = lambda local: math.pi+2.*(local/span-.5)
        model(right, hz, overtone=lambda n: .13*complex(math.cos(turn(n)), math.sin(turn(n))))
        self.call('--phrase', self.out, left, 0, LENGTH, *ENV_A)
        self.call('--phrase', self.out, right, 0, LENGTH, *ENV_B, 'LAST_NOTE')
        actual = read_wave(self.out)
        angle = lambda n: (n-begin)/span*turn(n-begin)
        for offset in (-12, -4, 0, 4, 12):
            with self.subTest(offset=offset):
                self.assert_overtone(actual, hz, begin+span//2+frames(offset), frames(5), angle)

    def assert_overtone(self, signal, hz, center, radius, relative_phase):
        real = quadrature = expected_real = expected_quadrature = 0.
        for n in range(center-radius, center+radius+1):
            phase = math.tau*hz*n/SR
            angle = 2.*phase+relative_phase(n)
            value = signal[n]/32768.-.32*math.cos(phase)
            window = .5+.5*math.cos(math.pi*(n-center)/radius)
            a, b = window*math.cos(angle), window*math.sin(angle)
            real += value*a
            quadrature += value*b
            expected_real += .26*math.cos(angle)*a
            expected_quadrature += .26*math.cos(angle)*b
        self.assertGreater(real/expected_real, .90)
        self.assertLess(abs(quadrature-expected_quadrature)/expected_real, .18)

    def test_zero_overlap_keeps_h1_phase_and_rest_has_no_oscillator(self):
        a, b = self.root/'left.wav', self.root/'right.wav'
        model(a, 237., 163)
        model(b, 237., 163)
        self.call('--phrase', self.out, a, 0, 163)
        self.call('--phrase', self.out, b, 0, 163)
        self.call('--phrase', self.out, 'R.wav', 0, 20, 0, 0, 'LAST_NOTE')
        actual = read_wave(self.out)
        self.assertEqual(len(actual), frames(346))
        seam = frames(163)
        expected = [round(32768*.32*math.cos(2*math.pi*237*i/SR)) for i in range(seam-150,seam+150)]
        self.assertLessEqual(max(abs(a-b) for a,b in zip(actual[seam-150:seam+150],expected)), 2)
        self.assertEqual(set(actual[-frames(20):]), {0})

    def test_resampler_removes_source_antiphase_in_native_overlap(self):
        source(self.root/'a.wav', 220., .3)
        source(self.root/'b.wav', 220., .3+math.pi)
        self.prepare('a.wav', 'a-cache.wav', 'A3')
        self.prepare('b.wav', 'b-cache.wav', 'A3')
        self.call('--phrase', self.out, 'a-cache.wav', 30, LENGTH, *ENV_A)
        self.call('--phrase', self.out, 'b-cache.wav', 40, LENGTH, *ENV_B, 'LAST_NOTE')
        actual = read_wave(self.out)
        def energy(center):
            part = actual[frames(center-10):frames(center+10)]
            return sum(x*x for x in part)/len(part)
        self.assertGreater(energy(430), .75*min(energy(250),energy(610)))
        self.assertLess(energy(430), 1.2*max(energy(250),energy(610)))

    def test_clock_advances_through_an_internal_confidence_gap(self):
        hz = 237.
        clip = self.root/'partial voicing.wav'
        model(clip, hz, support=lambda j: 0. if 28 <= j <= 35 else 1.)
        self.call('--phrase', self.out, clip, 0, LENGTH, 'LAST_NOTE')
        actual = read_wave(self.out)
        self.assertEqual(set(actual[frames(180):frames(204)]), {0})
        indices = range(frames(235), frames(250))
        error = max(abs(actual[n]-round(32768*.32*math.cos(math.tau*hz*n/SR))) for n in indices)
        self.assertLessEqual(error, 2)

    def test_classic_prepares_a_model_and_never_falls_back_to_pcm(self):
        source(self.root/'vowel.wav', 220., .3)
        note = self.root/'cached note.wav'
        self.call('vowel.wav', note, 'A3', 100, 'e', 0, 650, 0, 0, 60, 0, '!120')
        self.assertEqual(note.read_bytes()[:8], b'TNFMOD1\0')
        self.call(self.out, note, 0, LENGTH, ok=False)
        self.call('--wavtool', self.out, note, 0, LENGTH, 'LAST_NOTE', ok=False)
        self.assertFalse(self.out.exists())
        self.assertFalse((self.root/'phrase.wav.dat').exists())
        self.call('--phrase', self.out, note, 0, LENGTH)
        self.assertEqual(read_wave(self.out), ())
        self.assertFalse((self.root/'phrase.wav.dat').exists())
        self.call('--phrase', self.out, note, 0, LENGTH, 'LAST_NOTE')
        self.assertEqual(len(read_wave(self.out)), frames(2*LENGTH))

    def test_bad_note_can_be_retried_with_committed_phrase_preserved(self):
        a, b = self.root/'left.wav', self.root/'right.wav'
        model(a, 220.)
        model(b, 220.)
        self.call('--phrase', self.out, a, 30, LENGTH, *ENV_A)
        state = self.root/'phrase.wav.tnphrase'
        previous = state.read_bytes(), self.out.read_bytes()
        payload = b.read_bytes()
        b.write_bytes(payload[:-1])
        self.call('--phrase', self.out, b, 40, LENGTH, *ENV_B, 'LAST_NOTE', ok=False)
        self.assertEqual((state.read_bytes(), self.out.read_bytes()), previous)
        b.write_bytes(payload)
        self.call('--phrase', self.out, b, 40, LENGTH, *ENV_B, 'LAST_NOTE')
        self.assertEqual(len(read_wave(self.out)), frames(860))

    @unittest.skipUnless(os.name == 'nt', 'Windows batch fixture')
    def test_stock_batch_footer_cached_and_reused_notes(self):
        source(self.root/'vowel.wav', 220., .5)
        tool, sampler = 'tool copy.exe', 'sampler copy.exe'
        shutil.copyfile(EXE, self.root/tool)
        shutil.copyfile(EXE, self.root/sampler)
        common_setup = ['@set output=phrase.wav', '@set helper=temp_helper.bat', f'@set tool={tool}',
                        f'@set resamp={sampler}', '@set cachedir=.', '@del "%output%" 2>nul']
        footer = ['@if not exist "%output%.whd" goto E', '@if not exist "%output%.dat" goto E',
                  'copy /Y "%output%.whd" /B + "%output%.dat" /B "%output%"',
                  'del "%output%.whd"', 'del "%output%.dat"', ':E']
        resample = '@"%resamp%" %1 %temp% %2 %vel% %flag% %5 %6 %7 %8 %params%'
        join = '@"%tool%" "%output%" %temp% %stp% %3 %env%'
        expected = None
        for cache in (True, False):
            with self.subTest(cache=cache):
                helper = ['@if exist %temp% goto A', resample, ':A', join] if cache else [
                    '@set temp="temp$$$.wav"', resample, join]
                (self.root/'temp_helper.bat').write_text('\n'.join(helper)+'\n', encoding='ascii')
                commands = list(common_setup)
                for i,(pitch,offset,env) in enumerate(zip(('A3','D4'),(30,40),(ENV_A,ENV_B))):
                    commands.extend(['@set params=60 0 !120', '@set flag="e"', '@set vel=100',
                                     f'@set stp={offset}', '@set env='+' '.join(map(str,env)),
                                     f'@set temp="note {i}.wav"',
                                     f'@call %helper% "vowel.wav" {pitch} 480@120+50 0 0 650 0 0 {i}'])
                commands.append('@"%tool%" "%output%" "R.wav" 0 20 0 0')
                (self.root/'temp.bat').write_text('\n'.join(commands+footer)+'\n', encoding='ascii')
                def run_batch():
                    res = subprocess.run(['cmd.exe','/d','/c','temp.bat'],cwd=self.root,
                                         capture_output=True,timeout=30)
                    self.assertEqual(res.returncode, 0, res.stdout + res.stderr)
                    self.assertNotIn(b'tn_finds:', res.stderr)
                    self.assertFalse((self.root/'phrase.wav.tnphrase').exists())
                    return read_wave(self.out)
                actual = run_batch()
                self.assertEqual(len(actual), frames(880))
                self.assertEqual(set(actual[-frames(20):]), {0})
                phase_center = actual[frames(398):frames(462)]
                low, mid, high = peak_powers(phase_center, (220., 256.83, 293.66))
                self.assertGreater(mid, 15.*max(low,high))
                if cache:
                    for i in range(2):
                        self.assertEqual((self.root/f'note {i}.wav').read_bytes()[:8], b'TNFMOD1\0')
                    self.assertEqual(run_batch(), actual)
                    expected = actual
                else:
                    self.assertEqual((self.root/'temp$$$.wav').read_bytes()[:8], b'TNFMOD1\0')
                    self.assertEqual(actual, expected)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    options, rest = parser.parse_known_args()
    EXE = options.exe.resolve(strict=True)
    unittest.main(argv=[parser.prog, *rest], verbosity=2)
