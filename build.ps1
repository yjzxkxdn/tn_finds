# SPDX-License-Identifier: GPL-3.0-only
[CmdletBinding()]
param(
    [string]$Cxx = 'g++',
    [ValidateSet('Both','Resampler','Phrase')]
    [string]$Mode = 'Both',
    [ValidatePattern('^tn_finds(?:_[a-zA-Z0-9][a-zA-Z0-9.-]*)?\.exe\z')]
    [string]$OutputName
)
$ErrorActionPreference = 'Stop'
if ($OutputName -and $Mode -eq 'Both') {
    if ($PSBoundParameters.ContainsKey('Mode')) { throw '-OutputName requires a single build mode.' }
    $Mode = 'Phrase'
}
$compiler = (Get-Command -Name $Cxx -CommandType Application).Source
$target = & $compiler -dumpmachine
if ($LASTEXITCODE -ne 0 -or $target -notmatch '^x86_64-.*mingw') {
    throw 'Use an x86-64 MinGW g++ (-Cxx), or use CMake for another toolchain.'
}
Push-Location -LiteralPath $PSScriptRoot
try {
    $sources = @('src/main.cpp', 'src/wav.cpp', 'src/engine.cpp', 'src/pitch.cpp',
        'src/formant.cpp', 'src/glottal.cpp', 'src/harmonic_phase.cpp', 'src/harmonic_mask.cpp',
        'src/harmonic_ridges.cpp', 'src/harmonic_synth.cpp')
    $builds = if ($Mode -eq 'Both') { 'Resampler', 'Phrase' } else { $Mode }
    foreach ($build in $builds) {
        $name = if ($OutputName) { $OutputName } else { 'tn_finds_v0.4.7-' + $build.ToLowerInvariant() + '.exe' }
        $extra = if ($build -eq 'Phrase') {
            '-DTN_FINDS_DUAL_ROLE=1', 'src/wavtool.cpp', 'src/bat_host.cpp', 'src/phrase.cpp', 'src/phrase_io.cpp'
        } else { @() }
        & $compiler -std=c++17 -O3 -Wall -Wextra -Wpedantic -municode -static @sources @extra -o $name
        if ($LASTEXITCODE -ne 0) { throw "g++ failed with exit code $LASTEXITCODE." }
        & "./$name" --version
        if ($LASTEXITCODE -ne 0) { throw 'The built executable did not run.' }
    }
} finally {
    Pop-Location
}
