# SPDX-License-Identifier: GPL-3.0-only
"""Check Nebula's source grid, short inputs and OTO time origin."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import soundfile as sf


def run(exe, *args):
    return subprocess.run([str(exe), *map(str, args)], check=True, capture_output=True, text=True, timeout=30)


def analyze(exe, source, output, signal, rate):
    sf.write(source, signal, rate, 'FLOAT')
    run(exe, '--estimate-f0', source, output)
    time, hz = np.loadtxt(output, delimiter=',', skiprows=1, ndmin=2).T
    hop = 1 << int(np.ceil(np.log2(rate/200)))
    np.testing.assert_allclose(time, np.arange(len(hz))*hop/rate, atol=1e-14, rtol=0)
    assert np.isfinite(hz).all() and (hz >= 0).all()
    return time, hz


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exe', type=Path)
    exe = parser.parse_args().exe.resolve(strict=True)
    rate = 44100
    with tempfile.TemporaryDirectory(prefix='tn-f0-regression-') as tmp:
        source, csv, output = (Path(tmp)/name for name in ('source.wav', 'source.csv', 'output.wav'))
        time = np.arange(int(rate*.72))/rate
        x = .12*np.sin(2*np.pi*196*time)+.06*np.sin(4*np.pi*196*time+.2)
        t, hz = analyze(exe, source, csv, x, rate)
        keep = (.12 < t) & (t < .61)
        assert np.count_nonzero(keep) > 40
        np.testing.assert_allclose(hz[keep], 196., atol=1., rtol=0.)
        for signal in (np.ones(len(x))*.2, x[:int(.02*rate)]):
            _, hz = analyze(exe, source, csv, signal, rate)
            assert not np.any(hz)

        time = np.arange(rate)/rate
        freq = np.where(time < .40, 139., 280.)
        theta = 2*np.pi*np.cumsum(freq)/rate
        x = .17*np.sin(theta)+.07*np.sin(2*theta+.6)
        t, hz = analyze(exe, source, csv, x, rate)
        for left, right, f0 in ((.12, .29, 139.), (.60, .85, 280.)):
            np.testing.assert_allclose(hz[(left < t) & (t < right)], f0, atol=1., rtol=0.)

        # A 50 ms OTO segment below the whole-file guard must still use the
        # full file's F0 and the correct source time origin.
        result = run(exe, '--render', source, output, '--target-hz', '370', '--offset-ms', '700',
                     '--cutoff-ms', '-50', '--length-ms', '580', '--residual-mix', '100', '--report')
        report = json.loads(result.stdout)
        assert abs(report['source_pitch_hz']-280.) < 2., report
        assert report['f0_backend'] == 'nebula', report
        assert sf.info(output).frames == round(.580*rate)


if __name__ == '__main__':
    main()
