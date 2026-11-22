# SPDX-License-Identifier: GPL-3.0-only
"""Black-box PCM WAVTOOL tests. The standalone resampler is a separate build."""
from __future__ import annotations

import argparse
from fractions import Fraction
import io
import math
from pathlib import Path
import random
import struct
import subprocess
import tempfile
import unittest
import wave

SR = 44100
EXE: Path | None = None
RESAMPLER: Path | None = None


def frames(ms):
    # Nearest frame, with exact rational fixture times and half-away ties.
    n = Fraction(str(ms)) * SR / 1000
    return (-1 if n < 0 else 1) * int(abs(n) + Fraction(1, 2))


def noise(n, seed):
    rng = random.Random(seed)
    return tuple(rng.randrange(-12000, 12001) for _ in range(n))


def pack(samples):
    return struct.pack(f'<{len(samples)}h', *samples)


def linear_gain(samples, points):
    """Expected PCM for the strictly ordered envelope segments of a fixture."""
    values = []
    for (t0, v0), (t1, v1) in zip(points, points[1:]):
        a, b = frames(t0), frames(t1)
        assert b > a
        for i in range(a, b):
            level = ((b-i)*v0 + (i-a)*v1)/(100*(b-a))
            values.append(round(samples[i]*level))
    assert len(values) == len(samples)
    return tuple(values)


class TestWavtool(unittest.TestCase):
    def setUp(self):
        if EXE is None:
            self.skipTest('Run with --exe to select the dual-role tn_finds binary')
        temp = tempfile.TemporaryDirectory(prefix='tnf-wavtool-')
        self.addCleanup(temp.cleanup)
        self.work = Path(temp.name) / '混合 かな ü'
        self.work.mkdir()
        self.out = self.work / '合成 output.wav'
        self.head = Path(str(self.out)+'.whd')
        self.data = Path(str(self.out)+'.dat')

    def run_cli(self, *args, exe=None, ok=True):
        result = subprocess.run([str(exe or EXE), *map(str, args)], cwd=self.work,
                                capture_output=True, text=True, encoding='utf-8',
                                errors='replace', timeout=20)
        diagnostic = f'{args!r}\n{result.stdout}{result.stderr}'
        if ok:
            self.assertEqual(result.returncode, 0, diagnostic)
        else:
            self.assertNotEqual(result.returncode, 0, diagnostic)

    def join(self, source, offset, length, *envelope, last=False, output=None, ok=True):
        finish = ('LAST_NOTE',) if last else ()
        self.run_cli('--wavtool', output or self.out, source, offset, length,
                     *envelope, *finish, ok=ok)

    def source(self, name, samples, *, sr=SR, channels=1):
        path = self.work / name
        with wave.open(str(path), 'wb') as writer:
            writer.setparams((channels, 2, sr, 0, 'NONE', 'not compressed'))
            writer.writeframes(pack(samples))
        return path

    def decode(self, payload):
        self.assertEqual(payload[:4], b'RIFF')
        self.assertEqual(payload[8:12], b'WAVE')
        self.assertEqual(struct.unpack_from('<I', payload, 4)[0] + 8, len(payload))
        with wave.open(io.BytesIO(payload), 'rb') as reader:
            self.assertEqual((reader.getframerate(), reader.getnchannels(),
                              reader.getsampwidth(), reader.getcomptype()),
                             (SR, 1, 2, 'NONE'))
            count = reader.getnframes()
            raw = reader.readframes(count + 1)
            self.assertEqual(len(raw), 2*count)
        return struct.unpack(f'<{count}h', raw)

    def split_wave(self, output=None):
        output = output or self.out
        header = Path(str(output)+'.whd').read_bytes()
        data = Path(str(output)+'.dat').read_bytes()
        # The host copies these two files byte for byte. Do not assume a
        # 44-byte header: valid RIFF metadata may precede the data chunk.
        self.assertEqual(header[-8:-4], b'data')
        self.assertEqual(struct.unpack_from('<I', header, len(header)-4)[0], len(data))
        return header+data

    def samples_close(self, actual, expected, tolerance=1):
        self.assertEqual(len(actual), len(expected))
        for i, (a, b) in enumerate(zip(actual, expected)):
            if abs(a-b) > tolerance:
                self.fail(f'PCM frame {i}: got {a}, expected {b} (±{tolerance})')

    def test_explicit_pcm_trim_at_fractional_offset(self):
        x = noise(frames(60), 523)
        src = self.source('素材 原.wav', x)
        start, count = frames('2.5'), frames('23.25')
        self.join(src, '+2.5', '23.25')
        self.assertEqual(self.decode(self.split_wave()), x[start:start+count])
        self.assertEqual(self.decode(self.out.read_bytes()), ())

    def test_negative_offset_prepends_silence_and_short_source_zero_pads(self):
        x = noise(frames(20), 92)
        src = self.source('short.wav', x)
        self.join(src, -10, 50)
        self.assertEqual(self.decode(self.split_wave()), (0,)*frames(10) + x + (0,)*frames(20))

    def test_rate_conversion_and_stereo_mix_keep_the_fixed_output_format(self):
        # Well away from either recording edge, a ramp has a known sample
        # value at each new time. Its channels deliberately differ in slope.
        stereo = tuple(v for i in range(2205) for v in (-12000+8*i, 1200+2*i))
        src = self.source('22050 stereo.wav', stereo, sr=22050, channels=2)
        self.join(src, 20, 40)
        expected = tuple(round(-5400 + 5*(441 + i/2)) for i in range(frames(40)))
        self.samples_close(self.decode(self.split_wave()), expected, tolerance=2)

    def test_full_envelope_uses_linear_segments_and_extra_control_point(self):
        x = (10000,)*frames(150)
        src = self.source('envelope.wav', x)
        # p1, p2, p3, v1..v4, overlap, p4, p5, v5; first overlap is ignored.
        self.join(src, 0, 150, 10, 20, 30, 30, 85, 95, 20, 30, 20, 30, 55)
        expected = linear_gain(x, ((0, 0), (10, 30), (30, 85), (60, 55),
                                   (100, 95), (130, 20), (150, 0)))
        self.samples_close(self.decode(self.split_wave()), expected)

    def test_coincident_endpoints_use_the_values_facing_the_interior(self):
        x = (-10000,)*frames(90)
        src = self.source('endpoint.wav', x)
        self.join(src, 0, 90, 0, 0, 0, 20, 90, 75, 35, 0, 0, 0, 90)
        expected = linear_gain(x, ((0, 90), (90, 75)))
        self.samples_close(self.decode(self.split_wave()), expected)

    def test_overlap_replaces_stale_tail_and_negative_overlap_inserts_gap(self):
        a, b = (12000,)*frames(100), (24000,)*frames(10)
        first = self.source('long.wav', a)
        second = self.source('inside overlap.wav', b)
        self.join(first, 0, 100, 0, 0, 0, 100, 100, 100, 100, 50)
        self.assertEqual(self.decode(self.split_wave()), a)
        self.join(second, 0, 10, 0, 0, 0, 100, 100, 100, 100, 20)
        # The second note finishes before the previous end, not after it.
        # Saturate rather than wrapping when the two PCM contributions sum.
        expected = a[:frames(80)] + (32767,)*frames(10)
        self.assertEqual(self.decode(self.split_wave()), expected)
        third = self.source('after gap.wav', (-7000,)*frames(20))
        self.join(third, 0, 20, 0, 0, 0, 100, 100, 100, 100, -10, last=True)
        final = expected + (0,)*frames(10) + (-7000,)*frames(20)
        self.assertEqual(self.decode(self.split_wave()), final)
        self.assertEqual(self.decode(self.out.read_bytes()), final)

    def test_fractional_tempo_lengths_keep_a_phrase_clock(self):
        expected = ()
        for i in range(7):
            value = 8000 + 500*i
            src = self.source(f'tick {i}.wav', (value,)*frames(10))
            # 1 tick at 125 BPM is 1 ms; the correction makes it 1.1 ms.
            self.join(src, 0, '1@125+0.1')
            start, end = (frames(Fraction(11*k, 10)) for k in (i, i+1))
            expected += (value,)*(end-start)
        self.assertEqual(len(expected), 340)  # Not 7 * round(1.1 ms).
        self.assertEqual(self.decode(self.split_wave()), expected)

    def test_starting_new_host_session_replaces_previous_split_files(self):
        one = self.source('old.wav', noise(frames(90), 85))
        x = noise(frames(20), 201)
        two = self.source('new.wav', x)
        self.join(one, 0, 90)
        old_length = len(self.data.read_bytes())
        self.out.unlink()  # The host deletes only the marker on restart.
        self.join(two, 0, 20, 0, 0, 0, 100, 100, 100, 100, 10)
        self.assertLess(len(self.data.read_bytes()), old_length)
        self.assertEqual(self.decode(self.split_wave()), x)

    def test_classic_rest_forms_append_silence(self):
        x = noise(frames(20), 71)
        src = self.source('active.wav', x)
        self.join(src, 0, 20)
        rest = self.work / 'R.wav'
        self.join(rest, 0, 10)
        self.assertFalse(rest.exists())
        # The compact two-zero rest envelope also suppresses a real sample.
        rest = self.source(rest.name, (8000,)*frames(10))
        self.join(rest, 0, 10, 0, 0, last=True)
        self.assertEqual(self.decode(self.out.read_bytes()), x + (0,)*frames(20))

    def test_failed_calls_do_not_change_input_or_live_split_files(self):
        src = self.source('original.wav', noise(frames(30), 318))
        self.join(src, 0, 30)
        truncated = self.work/'truncated.wav'
        truncated.write_bytes(src.read_bytes()[:-2])
        before = {p: p.read_bytes() for p in (src, truncated, self.out, self.head, self.data)}
        for invalid_src, offset, length in ((src, 'nan', 20), (src, 0, -10),
                                            (src, 0, '480@0+0'), (src, 0, 1_000_000),
                                            (truncated, 0, 20), (self.work/'missing.wav', 0, 20)):
            with self.subTest(offset=offset, length=length, source=invalid_src.name):
                self.join(invalid_src, offset, length, ok=False)
                self.assertEqual({p: p.read_bytes() for p in before}, before)
        self.join(src, 0, 20, output=src, ok=False)
        self.assertEqual(src.read_bytes(), before[src])
        self.assertFalse(Path(str(src)+'.whd').exists())
        self.assertFalse(Path(str(src)+'.dat').exists())

    def test_standalone_resampler_feeds_pcm_wavtool(self):
        if RESAMPLER is None:
            self.skipTest('Select the standalone build with --resampler')
        first = self.source('oto A.wav', noise(frames(120), 369))
        second = self.source('oto B.wav', noise(frames(120), 249))
        cache = [self.work/'cache 1.wav', self.work/'cache 2.wav']
        for src, out, tone in zip((first, second), cache, ('C4', 'D4')):
            # Nonvoiced fixtures keep the join transparent and independent
            # of the phase-aware treatment of sustained notes.
            self.run_cli(src.name, out.name, tone, 100, 'eW-1', 0, 120,
                         0, 0, 50, 0, '!120', exe=RESAMPLER)
        a, b = (self.decode(p.read_bytes()) for p in cache)
        # Emulate the host: 11 envelope numbers, two tempo-corrected notes,
        # and binary assembly without the standalone LAST_NOTE extension.
        self.join(cache[0].name, 10, '96@120-10', 0, 0, 10, 100, 100, 100, 0, 0, 0, 0, 100)
        self.join(cache[1].name, 0, '96@120+10', 0, 10, 0, 0, 100, 100, 100, 10, 0, 0, 100)
        length = frames(90)
        left = linear_gain(a[frames(10):frames(10)+length], ((0, 100), (80, 100), (90, 0)))
        right = linear_gain(b[:frames(110)], ((0, 0), (10, 100), (110, 100)))
        ov = frames(10)
        expected = left[:-ov] + tuple(x+y for x, y in zip(left[-ov:], right[:ov])) + right[ov:]
        self.out.write_bytes(self.split_wave())
        self.samples_close(self.decode(self.out.read_bytes()), expected)

    def test_plain_pcm_overlap_preserves_antiphase_cancellation(self):
        n, ov = frames(160), frames(60)
        start = n-ov

        def tone(i):
            phi = 2*math.pi*i/210
            return round(10000*math.sin(phi + .39) + 2000*math.sin(3*phi - .64))

        a = tuple(tone(i) for i in range(n))
        b = tuple(-tone(start+i) for i in range(n))
        first = self.source('steady A.wav', a)
        second = self.source('antiphase B.wav', b)
        self.join(first, 0, 160, 0, 0, 60, 100, 100, 100, 0, 0)
        self.join(second, 0, 160, 0, 60, 0, 0, 100, 100, 100, 60)
        result = self.decode(self.split_wave())
        self.assertEqual(len(result), 2*n-ov)
        self.assertEqual(result[:start], a[:start])
        self.assertEqual(result[n:], b[ov:])

        # A fixed linear sum has an exact null in the middle. Check a band
        # of samples rather than testing the phase of one arbitrary sample.
        lo, hi = frames(24), frames(36)
        classic = [(1.-2.*i/ov)*a[start+i] for i in range(lo, hi)]
        self.samples_close(result[start+lo:start+hi], tuple(round(x) for x in classic))

    def test_plain_pcm_overlap_preserves_both_endpoint_pitches(self):
        length, overlap = frames(550), frames(240)
        first_hz, last_hz = 220., 293.6647679

        def tone(hz, phase):
            return tuple(round(10000*math.sin(2*math.pi*hz*i/SR+phase)) for i in range(length))

        left = self.source('A3.wav', tone(first_hz, .4))
        right = self.source('D4.wav', tone(last_hz, 1.7))
        self.join(left, 0, 550, 0, 0, 240, 100, 100, 100, 0, 0)
        self.join(right, 0, 550, 0, 240, 0, 0, 100, 100, 100, 240)
        data = self.decode(self.split_wave())
        self.assertEqual(len(data), 2*length-overlap)
        center = length-overlap//2
        half = frames(32)
        signal = data[center-half:center+half]

        def power(hz):
            frequency = 2*math.pi*hz/SR
            samples = [x*(.5-.5*math.cos(2*math.pi*i/(len(signal)-1))) for i, x in enumerate(signal)]
            real = sum(x*math.cos(i*frequency) for i, x in enumerate(samples))
            imag = sum(x*math.sin(i*frequency) for i, x in enumerate(samples))
            return real*real+imag*imag

        midpoint = (first_hz+last_hz)*.5
        self.assertGreater(min(power(first_hz), power(last_hz)), 15.*power(midpoint))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--resampler', type=Path)
    options, rest = parser.parse_known_args()
    EXE = options.exe.resolve(strict=True)
    RESAMPLER = options.resampler.resolve(strict=True) if options.resampler else None
    unittest.main(argv=[parser.prog, *rest], verbosity=2)
