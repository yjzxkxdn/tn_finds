# SPDX-License-Identifier: GPL-3.0-only
"""Default MRQ precedence, source cache keys, and the OpenUtau alias handoff."""
import csv
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import wave


RATE = 44100
N = int(.82*RATE)
HOP = 256


def recording(path, hz=197):
    def sample(i):
        a = i*math.tau*hz/RATE
        gate = min(1., i/850., (N-1-i)/850.)
        return round(gate*(7800.*math.sin(a)+2700.*math.sin(2*a+.46)+1600.*math.sin(3*a)))

    with wave.open(str(path), 'wb') as out:
        out.setparams((1, 2, RATE, N, 'NONE', ''))
        out.writeframes(struct.pack(f'<{N}h', *(sample(i) for i in range(N))))


def entry(name, values, *, hop=HOP, rate=RATE, tail=True):
    name = name.encode('utf-16-le')
    payload = struct.pack('<iii', len(values), rate, hop)+struct.pack(f'<{len(values)}f', *values)
    if tail:
        payload += struct.pack('<ii', 1, 0)  # old hand-edited track: use it anyway
    return struct.pack('<i', len(name)//2)+name+struct.pack('<i', len(payload))+payload


def archive(path, *records, version=2):
    path.write_bytes(struct.pack('<4sii', b'mrq ', version, len(records))+b''.join(records))


def alias(original):
    # OpenUtau emits XXH32(UTF8(oto.File), seed=0) as the final component.
    # Minimal source-independent test implementation of that host hash.
    data = str(original).encode('utf-8')
    mask = (1 << 32)-1
    p1, p2, p3, p4, p5 = 2654435761, 2246822519, 3266489917, 668265263, 374761393

    def rot(n, bits):
        n &= mask
        return ((n << bits) | (n >> (32-bits))) & mask

    at = 0
    h = p5
    if len(data) >= 16:
        lanes = [p1+p2, p2, 0, -p1]
        while at+16 <= len(data):
            for j, word in enumerate(struct.unpack_from('<4I', data, at)):
                lanes[j] = rot(lanes[j]+word*p2, 13)*p1 & mask
            at += 16
        h = sum(rot(v, s) for v, s in zip(lanes, (1, 7, 12, 18)))
    h = (h+len(data)) & mask
    while at+4 <= len(data):
        h = rot(h+struct.unpack_from('<I', data, at)[0]*p3,17)*p4 & mask
        at += 4
    for b in data[at:]:
        h = rot(h+b*p5,11)*p1 & mask
    h = (h^(h >> 15))*p2 & mask
    h = (h^(h >> 13))*p3 & mask
    return f'src-01234567-89abcdef-{h^(h >> 16):08x}.wav'


def exercise(root, resampler, phrase):
    bank = root/'声库 保存'
    bank.mkdir()
    original = bank/'女声_A.wav'
    mrq = bank/'desc.mrq'
    f0_root = root/'engine-cache'
    analysis_root = root/'analysis-cache'
    bindings_root = root/'bind-cache'
    env = {**os.environ, 'TN_FINDS_NEBULA_CACHE': str(f0_root), 'TN_FINDS_ANALYSIS_CACHE': str(analysis_root),
           'TN_FINDS_MRQ_BINDINGS': str(bindings_root), 'OMP_NUM_THREADS': '2'}
    env.pop('TN_FINDS_MRQ', None)
    recording(original)

    def run(exe, *args, ok=True, settings=None):
        proc = subprocess.run([str(exe), *map(str, args)], cwd=root, env={**env, **(settings or {})},
                              capture_output=True, timeout=30)
        assert (proc.returncode == 0) == ok, (args, proc.returncode, proc.stdout, proc.stderr)
        return proc.stdout

    def export(source=original, *flags, exe=resampler, **kw):
        path = root/'export.csv'
        run(exe, '--estimate-f0', source, path, *flags, **kw)
        with path.open(newline='', encoding='ascii') as data:
            return [(float(row['time_s']), float(row['f0_hz'])) for row in csv.DictReader(data)]

    opts = ('--offset-ms', 120, '--cutoff-ms', -148, '--consonant-ms', 144, '--length-ms', 700,
            '--target-hz', 249, '--fixed-residual-mix', 15, '--residual-mix', 100)

    def render(*flags, source=original, settings=None, options=opts):
        path = root/'sound.wav'
        report = json.loads(run(resampler, '--render', source, path, *options, '--report', *flags, settings=settings))
        return path.read_bytes(), report

    # The direct MinGW script builds without the persistent analysis cache; it
    # still has to render identically, just without a second-call hit.
    cached = False

    def cache_available():
        probe = root/'cache-probe'
        render(settings={'TN_FINDS_ANALYSIS_CACHE': str(probe)})
        return probe.exists() and any(probe.rglob('*'))

    def note(name, pitch, source=original):
        path = root/name
        run(phrase, '--prepare-note', source, path, pitch, 100, 'eTnc15Tns100', 120, 700, 144, -148, 100, 0)
        assert path.read_bytes().startswith(b'TNFMOD1\0')
        return path

    # Match the whole UTF-16 filename, not the first entry or a tombstone;
    # timestamps and the hand-edited mark do not affect precedence.
    wrong = entry(original.name.upper(), [320.]*(N//HOP))
    tombstone = entry('\0'*len(original.name), [0.]*(N//HOP))
    values = [0.]*5+[194.5]*(N//HOP-9)+[0.]*4
    archive(mrq, wrong, tombstone, entry(original.name, values))
    initial = mrq.read_bytes(), mrq.stat().st_mtime_ns
    expected = [(i*HOP/RATE, v) for i, v in enumerate(values)]
    assert export() == expected
    assert export(exe=phrase) == expected
    audio, report = render()
    assert report['f0_backend'] == 'mrq' and 192. < report['source_pitch_hz'] < 197., report
    assert (mrq.read_bytes(), mrq.stat().st_mtime_ns) == initial
    assert not f0_root.exists()  # neither render nor export invoked the auto estimator
    assert render(settings={'TN_FINDS_ANALYSIS_CACHE': 'off'})[0] == audio
    if sys.platform == 'win32':
        cached = cache_available()
        rendered, report = render()
        assert rendered == audio and report['analysis_cache_hit'] == cached, report

    first = note('mrq-note-a.model', 'B3')
    second = note('mrq-note-b.model', 'D4')
    phrase_out = root/'whole.wav'
    run(phrase, '--phrase', phrase_out, first, 0, 650, 0, 0, 90, 100, 100, 100, 0)
    run(phrase, '--phrase', phrase_out, second, 0, 680, 0, 60, 0, 0, 100, 100, 100, 70, 'LAST_NOTE')
    with wave.open(str(phrase_out), 'rb') as out:
        assert out.getframerate() == RATE and out.getnframes() > RATE
    assert not f0_root.exists()

    # The same length and restored timestamp do not keep old DSP in cache.
    changed_values = [192.]*(N//HOP)
    archive(mrq, wrong, tombstone, entry(original.name, changed_values))
    os.utime(mrq, ns=(initial[1], initial[1]))
    new_audio, report = render()
    assert new_audio != audio and report['f0_backend'] == 'mrq' and not report['analysis_cache_hit']
    assert new_audio == render(settings={'TN_FINDS_ANALYSIS_CACHE': 'off'})[0]
    archive(mrq, wrong, tombstone, entry(original.name, [0.]*(N//HOP)))
    uv = export()
    assert len(uv) == N//HOP and all(hz == 0. for _, hz in uv)
    _, report = render()
    assert report['f0_backend'] == 'mrq' and report['source_pitch_hz'] == 0. and report['islands'] == 0, report
    zero = note('mrq-uv.model', 'D4')
    assert struct.unpack_from('<I', zero.read_bytes(), 28)[0] == 0
    assert not f0_root.exists()

    # A zero span inside a periodic recording also stays out of the target
    # harmonic stem, even though the alignment clock bridges the gap.
    mixed_uv = [197. if not .30 < i*HOP/RATE < .52 else 0. for i in range(N//HOP)]
    archive(mrq, entry(original.name, mixed_uv))
    middle = root/'mrq-middle-uv.model'
    run(phrase, '--prepare-note', original, middle, 'D4', 100, 'e', 0, 820, 0, 0, 100, 0)
    model = middle.read_bytes()
    _, active, hop, pitch_hop, samples, frames = struct.unpack_from('<4i2I', model, 8)
    at = 32+4*(samples+1+(active+pitch_hop-1)//pitch_hop)
    support = []
    for _ in range(frames):
        _, weight, orders = struct.unpack_from('<ffI', model, at)
        at += 12+8*orders
        support.append(weight)
    assert max(support[:int(.25*RATE/hop)]) > .1 and max(support[int(.6*RATE/hop):]) > .1
    assert max(support[int(.36*RATE/hop):int(.46*RATE/hop)]) == 0.
    assert not f0_root.exists()

    # Non-Nebula hop is preserved. Mixed recorded and unvoiced spans select
    # the supplied source-time section, not the new output note or OTO origin.
    hop = 1536
    native = [186.]*(N//hop//2)+[199.]*(N//hop-N//hop//2)
    archive(mrq, entry(original.name, native, hop=hop, tail=False), version=1)
    assert export() == [(i*hop/RATE, v) for i, v in enumerate(native)]
    high_offset = ('--offset-ms', 560, '--cutoff-ms', -210, '--consonant-ms', 50,
                   '--length-ms', 340, '--target-hz', 279)
    _, report = render(options=high_offset)
    assert report['f0_backend'] == 'mrq' and 198. < report['source_pitch_hz'] < 200., report
    if sys.platform == 'win32' and cached:
        assert render(options=high_offset)[1]['analysis_cache_hit']

    explicit = root/'known.csv'
    explicit.write_text('time_s,f0_hz\n0,198\n1,198\n', encoding='ascii')
    assert render('--f0', explicit)[1]['f0_backend'] == 'csv'
    assert render('--forced-hz', 200)[1]['f0_backend'] == 'forced'
    assert render('--forced-hz', -1)[1]['f0_backend'] == 'forced'
    assert not f0_root.exists()

    estimated = export(original, '--no-mrq')
    assert export(settings={'TN_FINDS_MRQ': 'off'}) == estimated and estimated != export()
    auto, fallback = render('--no-mrq')
    assert fallback['f0_backend'] in ('nebula', 'standard')
    for malformed in [entry('unrelated.wav', values), entry(original.name, values, rate=48000),
                      entry(original.name, [190.]), entry(original.name, [float('nan')]*len(values)),
                      entry(original.name, values)[:-2]]:
        archive(mrq, malformed)
        assert export() == estimated
        a, r = render()
        assert a == auto and r['f0_backend'] == fallback['f0_backend']
    mrq.unlink()
    assert export() == estimated
    archive(mrq, entry(original.name, values))
    assert export() == expected  # use it even though auto F0 was already cached
    preserved = mrq.read_bytes()
    run(resampler, '--estimate-f0', original, mrq, ok=False)
    assert mrq.read_bytes() == preserved

    if sys.platform == 'win32':
        temp = root/'OpenUtau Cache'
        temp.mkdir()
        copied = temp/alias(original)
        shutil.copyfile(original, copied)
        assert export(copied) == estimated  # no binding: original path not in the host call
        host_mrq = temp/'desc.mrq'
        alternate = [198.]*(N//HOP)
        archive(host_mrq, entry(copied.name, alternate))
        assert [hz for _, hz in export(copied)] == alternate
        assert b'Bound 1 MRQ source' in run(phrase, '--bind-mrq', bank)
        assert export(copied) == expected  # registered original beats stale host MRQ
        assert note('same-bound.model', 'D4', copied).read_bytes() == second.read_bytes()
        archive(mrq, entry(original.name, alternate))
        assert [hz for _, hz in export(copied)] == alternate
        archive(mrq, entry(original.name, values))
        host_mrq.unlink()
        recording(copied, 223)
        assert export(copied) != expected
        shutil.copyfile(original, copied)
        assert export(copied) == expected
        assert export(copied, settings={'TN_FINDS_MRQ_BINDINGS': 'off'}) == estimated
    # Producer side: emit the descriptor a Moresampler-style bank ships, then
    # read it back through the normal path.
    def walk(blob):
        version, count = struct.unpack_from('<ii', blob, 4)
        at, rows = 12, []
        for _ in range(count):
            units = struct.unpack_from('<i', blob, at)[0]
            name = blob[at+4:at+4+2*units].decode('utf-16-le')
            payload = struct.unpack_from('<i', blob, at+4+2*units)[0]
            trunk = at+8+2*units
            nf0, fs, hop = struct.unpack_from('<iii', blob, trunk)
            rows.append(dict(name=name, trunk=trunk, payload=payload, nf0=nf0, fs=fs, hop=hop))
            at = trunk+payload
        assert at == len(blob), (at, len(blob))
        return version, rows

    def longer(path, factor=2):
        with wave.open(str(path), 'wb') as out:
            out.setparams((1, 2, RATE, N*factor, 'NONE', ''))
            out.writeframes(struct.pack(f'<{N*factor}h', *(int(7000*math.sin(i*math.tau*233/RATE))
                                                            for i in range(N*factor))))

    written = root/'written'
    written.mkdir()
    source = written/'voice.wav'
    recording(source, 213)
    produce = subprocess.run([str(resampler), '--write-mrq', str(source)], cwd=root, env=env,
                             capture_output=True, timeout=60)
    if produce.returncode != 0:
        # Builds without the Nebula backend must refuse: their trajectory has a
        # different grid and voicing, so a desc.mrq written from it would claim
        # to be a Moresampler-compatible measure of the recording.
        assert b'Nebula' in produce.stderr and not (written/'desc.mrq').exists(), produce.stderr
    else:
        desc = written/'desc.mrq'
        blob = desc.read_bytes()
        version, rows = walk(blob)
        assert version == 2 and len(rows) == 1 and rows[0]['name'] == source.name
        # The hop belongs to whichever backend measured the recording, so only
        # its range is fixed here; the payload framing is not backend dependent.
        assert rows[0]['payload'] == 20+4*rows[0]['nf0']
        assert rows[0]['fs'] == RATE and 0 < rows[0]['hop'] <= RATE
        assert struct.unpack_from('<i', blob, rows[0]['trunk']+12+4*rows[0]['nf0']+4)[0] == 524
        f0 = struct.unpack_from(f"<{rows[0]['nf0']}f", blob, rows[0]['trunk']+12)
        exported = [hz for _, hz in export(source)]
        assert len(exported) == len(f0) and max(abs(a-b) for a, b in zip(exported, f0)) < 1e-3

        # Same trajectory again fits the old payload: rewrite in place, no growth.
        assert b'Wrote' in run(resampler, '--write-mrq', source)
        assert desc.stat().st_size == len(blob)
        assert len(walk(desc.read_bytes())[1]) == 1

        # A different source appends rather than overwriting.
        other = written/'second.wav'
        recording(other, 251)
        run(resampler, '--write-mrq', other)
        version, rows = walk(desc.read_bytes())
        assert version == 2 and [row['name'] for row in rows] == [source.name, other.name]

        # A longer recording no longer fits, so the old record is tombstoned and
        # a fresh one appended - the same rule Moresampler applies.
        longer(written/'voice.wav')
        run(resampler, '--write-mrq', written/'voice.wav')
        version, rows = walk(desc.read_bytes())
        assert version == 2 and len(rows) == 3
        assert rows[0]['name'] == '\0'*len(source.name)  # tombstoned, name zeroed
        assert [rows[1]['name'], rows[2]['name']] == [other.name, source.name]

        # And the reader really consumes that file: a patched frame shows up.
        patched = bytearray(desc.read_bytes())
        struct.pack_into('<f', patched, rows[2]['trunk']+12+4*7, 321.5)
        desc.write_bytes(bytes(patched))
        assert abs([hz for _, hz in export(written/'voice.wav')][7] - 321.5) < 1e-6

    assert sorted(p.name for p in bank.iterdir()) == sorted((original.name, mrq.name))


if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='tn MRQ 测试 ') as tmp:
        exercise(Path(tmp).resolve(), *(Path(v).resolve(strict=True) for v in sys.argv[1:3]))
