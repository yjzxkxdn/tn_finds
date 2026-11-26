# SPDX-License-Identifier: GPL-3.0-only
"""Cross-process Nebula reuse, content invalidation, and cache failure fallback."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import soundfile as sf


def call(exe, root, cache, *args):
    env = {**os.environ, 'TN_FINDS_NEBULA_CACHE': str(cache), 'OMP_NUM_THREADS': '2'}
    subprocess.run([str(exe), *map(str, args)], cwd=root, env=env,
                   check=True, capture_output=True, timeout=60)


def main():
    exe = Path(sys.argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix='tn Nebula 测试 ') as tmp:
        root = Path(tmp)
        source = root/'原音.wav'
        sr = 44100
        t = np.arange(sr)/sr

        def record(hz, rate=sr):
            x = .2*np.sin(2*np.pi*hz*t)+.07*np.sin(4*np.pi*hz*t+.37)
            sf.write(source, x, rate, 'PCM_16')

        with_cache = root/'cache'
        cache_dir = with_cache/'tn_finds'/'nebula-f0-v1'

        def estimate(cache, name, wav=source):
            call(exe, root, cache, '--estimate-f0', wav, name)
            return (root/name).read_bytes()

        record(167)
        original = estimate('off', 'cold.csv')
        with ThreadPoolExecutor(max_workers=4) as pool:
            tracks = list(pool.map(lambda i: estimate(with_cache, f'parallel-{i}.csv'), range(4)))
        assert all(track == original for track in tracks)
        entries = list(cache_dir.iterdir())
        assert len(entries) == 1 and entries[0].suffix == '.f0', entries
        entry = entries[0]
        mark = entry.stat().st_mtime_ns
        renamed = root/'same decoded audio.wav'
        shutil.copyfile(source, renamed)
        assert estimate(with_cache, 'renamed.csv', renamed) == original
        assert list(cache_dir.iterdir()) == entries and entry.stat().st_mtime_ns == mark

        # A changed payload, including a plausible positive F0, fails integrity
        # validation. A fresh estimate must replace that partial/corrupt result.
        good = entry.read_bytes()
        entry.write_bytes(good[:-1] + bytes([good[-1] ^ 0x01]))
        assert estimate(with_cache, 'repair.csv') == original
        assert entry.read_bytes() == good
        entry.write_bytes(good[:90])
        assert estimate(with_cache, 'truncated.csv') == original
        assert entry.read_bytes() == good

        # Keep the path, size and mtime identical while changing real audio.
        meta = source.stat()
        record(238)
        os.utime(source, ns=(meta.st_atime_ns, meta.st_mtime_ns))
        assert source.stat().st_size == meta.st_size
        edited = estimate(with_cache, 'edited.csv')
        assert edited != original and edited == estimate('off', 'edited-off.csv')
        assert len(list(cache_dir.iterdir())) == 2
        # The same decoded samples on a different timeline also miss.
        record(238, 48000)
        assert estimate(with_cache, 'rate.csv') == estimate('off', 'rate-off.csv')
        assert len(list(cache_dir.iterdir())) == 3

        # An unavailable cache never blocks the user from rendering.
        unusable = root/'not-a-directory'
        unusable.write_bytes(b'no cache here')
        assert estimate(unusable, 'unavailable.csv') == (root/'rate-off.csv').read_bytes()
        # Bypass and explicitly supplied source F0 must not generate a cache.
        unused = root/'unused'
        call(exe, root, unused, '--render', source, 'bypass.wav')
        call(exe, root, unused, '--render', source, 'forced.wav', '--target-hz', 220, '--forced-hz', 288)
        call(exe, root, unused, '--render', source, 'csv.wav', '--target-hz', 220, '--f0', 'rate-off.csv')
        assert not unused.exists()

        opts = ('--offset-ms', 280, '--cutoff-ms', -85, '--consonant-ms', 80,
                '--length-ms', 630, '--target-hz', 220, '--residual-mix', 100)
        call(exe, root, 'off', '--render', source, 'direct.wav', *opts)
        call(exe, root, with_cache, '--render', source, 'cached.wav', *opts)
        assert hashlib.sha256((root/'direct.wav').read_bytes()).digest() == hashlib.sha256((root/'cached.wav').read_bytes()).digest()
        assert sorted(p.suffix for p in cache_dir.iterdir()) == ['.f0']*3


if __name__ == '__main__':
    main()
