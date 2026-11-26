# SPDX-License-Identifier: GPL-3.0-only
"""Exercise exact OTO-local cache reuse, misses, and concurrent publication."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile

import numpy as np
import soundfile as sf


def exercise(resampler, phrase, root):
    bank = root/'声音库'
    bank.mkdir()
    src = bank/'原音.wav'
    rate = 44100
    t = np.arange(int(.81*rate))/rate
    rng = np.random.default_rng(19)
    noise = rng.standard_normal(t.size)

    def record(hz, sr=rate):
        cycles = hz*t+.03*np.sin(2*np.pi*4*t)
        tone = .22*np.sin(2*np.pi*cycles) + .05*np.sin(4*np.pi*cycles+.31)
        sf.write(src, (.75+.25*np.cos(2*np.pi*t))*tone+.011*noise, sr, 'PCM_16')

    record(197.)
    cache = root/'分析'
    entries = cache/'tn_finds'/'oto-analysis-v1'
    base = {'--offset-ms': 120, '--cutoff-ms': -148, '--consonant-ms': 138,
            '--length-ms': 730, '--target-hz': 249, '--residual-mix': 100, '--fixed-residual-mix': 21}

    def run(exe, where, *args, threads=2):
        env = {**os.environ, 'TN_FINDS_NEBULA_CACHE': 'off', 'TN_FINDS_ANALYSIS_CACHE': str(where),
               'OMP_NUM_THREADS': str(threads)}
        return subprocess.run([str(exe), *map(str, args)], cwd=root, env=env,
                              capture_output=True, check=True, timeout=60).stdout

    def render(name, where=cache, source=src, flags=(), **options):
        dst = root/(name+'.wav')
        opts = {**base, **options}
        report = json.loads(run(resampler, where, '--render', source, dst, '--report',
                                *(value for item in opts.items() for value in item), *flags))
        return dst.read_bytes(), report

    # The direct MinGW script links no persistent cache; there is nothing to
    # reuse, so report that instead of failing the reuse assertions.
    probe = root/'probe'
    render('probe', probe)
    if not any(probe.rglob('*.ana')):
        print('analysis cache is not compiled into this build; skipped')
        return

    expected, off = render('off', 'off')
    assert off['source_pitch_hz'] > 0 and not off['analysis_cache_hit'], off
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(lambda i: render(f'parallel-{i}'), range(4)))
    assert all(audio == expected for audio, _ in results)
    assert sum(not req['analysis_cache_hit'] for _, req in results) == 1, results[0][1]
    [entry] = entries.iterdir()
    assert entry.suffix == '.ana'
    stamp = entry.stat().st_mtime_ns
    renamed = bank/'同一段录音.wav'
    shutil.copyfile(src, renamed)
    changed_output = {'--target-hz': 310, '--length-ms': 820, '--fixed-residual-mix': 80,
                      '--volume': 81, '--formant-st': -1.5, '--tension': 27}
    edited, hit = render('modified', source=renamed, **changed_output)
    assert edited != expected and hit['analysis_cache_hit'], hit
    assert edited == render('modified-uncached', 'off', **changed_output)[0]
    assert entry.stat().st_mtime_ns == stamp and list(entries.iterdir()) == [entry]

    def model(note, name, where):
        dst = root/(name+'.model')
        run(phrase, where, '--prepare-note', src, dst, note, 100, 'eTnc21Tns100', 120, 730, 138, -148, 100, 0)
        return dst

    def assemble(name, notes):
        dst = root/(name+'.wav')
        run(phrase, 'off', '--phrase', dst, notes[0], 0, 560, 0, 0, 90, 100, 100, 100, 0)
        run(phrase, 'off', '--phrase', dst, notes[1], 0, 590, 0, 60, 0, 0, 100, 100, 100, 70, 'LAST_NOTE')
        return dst.read_bytes()

    warm = [model(pitch,name,cache) for pitch,name in [('B3','ref-model'), ('D4','edit-model')]]
    fresh = [model(pitch,name,'off') for pitch,name in [('B3','ref-fresh'), ('D4','edit-fresh')]]
    assert [p.read_bytes() for p in warm] == [p.read_bytes() for p in fresh]
    assert warm[0].read_bytes() != warm[1].read_bytes()
    assert list(entries.iterdir()) == [entry] and entry.stat().st_mtime_ns == stamp
    edit_phrase = assemble('phrase-edited', warm)
    assert edit_phrase == assemble('phrase-edited-off', fresh)
    assert edit_phrase != assemble('phrase-original', [warm[0], warm[0]])

    # Same (offset, cutoff, fixed) but new output duration can change which
    # source-side context is available to the short-tail guard.
    cases = [('new-guard', {'--length-ms': 154}, ()), ('new-oto', {'--cutoff-ms': -160}, ()),
             ('new-front', {'--offset-ms': 121}, ()),
             ('new-mask', {'--mask-return': .7500000001}, ()), ('no-rd', {}, ('--phase-lock',)),
             ('no-cleanup', {'--offset-ms': 0, '--cutoff-ms': 0, '--length-ms': 810,
                             '--target-hz': 0, '--residual-mix': 100}, ()),
             ('forced', {'--forced-hz': 201}, ())]
    for name, opts, flags in cases:
        audio, report = render(name, flags=flags, **opts)
        assert not report['analysis_cache_hit'], (name, report)
        assert audio == render(name+'-off', 'off', flags=flags, **opts)[0]
        again, report = render(name+'-second', flags=flags, **opts)
        assert again == audio and report['analysis_cache_hit'], (name, report)
    (root/'external.csv').write_text('time_s,f0_hz\n0,197\n1,197\n', encoding='ascii')
    csvpath = root/'external.csv'
    flags = ('--f0', csvpath)
    ext, report = render('csv', flags=flags)
    assert not report['analysis_cache_hit']
    assert ext == render('csv-second', flags=flags)[0]
    csvpath.write_text('time_s,f0_hz\n0,199\n1,199\n', encoding='ascii')
    assert not render('csv-edited', flags=flags)[1]['analysis_cache_hit']

    good = entry.read_bytes()
    # Change a payload byte, truncate, then inject an unchecked dimension
    # with an otherwise intact entry. Reads must fail safely and republish.
    header = bytearray(good)
    struct.pack_into('<I', header, 76, 0xffffffff)
    for i, bad in enumerate([good[:-1]+bytes([good[-1]^1]), good[:120], bytes(header)]):
        entry.write_bytes(bad)
        audio, report = render(f'repair-{i}')
        assert audio == expected and not report['analysis_cache_hit']
        assert entry.read_bytes() == good

    # Standard/supplied analysis only depends on the slice. Nebula may use
    # the rest of the original recording when resolving source F0.
    exterior = root/'另外一段录音.wav'
    pcm, _ = sf.read(src, dtype='int16')
    pcm[-2200:] = -pcm[-2200:]
    sf.write(exterior, pcm, rate, 'PCM_16')
    audio, report = render('outside-oto', source=exterior)
    assert report['analysis_cache_hit'] == (report['f0_backend'] != 'nebula')
    assert audio == render('outside-oto-off', 'off', exterior)[0]

    meta = src.stat()
    record(231.)
    os.utime(src, ns=(meta.st_atime_ns, meta.st_mtime_ns))
    assert src.stat().st_size == meta.st_size
    res, report = render('same-mtime')
    assert not report['analysis_cache_hit'] and res != expected
    assert res == render('same-mtime-off', 'off')[0]
    record(231., 48000)
    res, report = render('different-rate')
    assert not report['analysis_cache_hit'] and res == render('rate-off', 'off')[0]

    unusable = root/'not-a-directory'
    unusable.write_bytes(b'This is not a directory.')
    got, report = render('unavailable', unusable)
    assert got == res and not report['analysis_cache_hit']
    unused = root/'unused'
    run(resampler, unused, '--render', src, root/'bypass.wav')
    assert not render('explicit-uv', unused, **{'--forced-hz': -1})[1]['analysis_cache_hit']
    assert not unused.exists()
    assert sorted(p.name for p in bank.iterdir()) == sorted((src.name, renamed.name))
    assert all(p.suffix == '.ana' for p in entries.iterdir())
    assert entry.read_bytes() == good


if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='tn OTO 分析 ') as tmp:
        exercise(*(Path(x).resolve(strict=True) for x in sys.argv[1:3]), Path(tmp))
