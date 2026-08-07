# PluginDNA

**PluginDNA is an experimental audio-plugin measurement and DSP research environment for understanding what processors actually do to a signal.**

Rather than relying on descriptions such as *warm*, *glue*, *analogue* or *console colour*, PluginDNA runs repeatable measurements across plugin parameters and input levels and records the resulting behaviour as structured evidence.

The long-term aim is to make it possible to examine individual processors, complete signal chains and console-style architectures in terms of measurable DSP behaviour.

## What PluginDNA measures

PluginDNA currently includes analysis of:

* RMS, peak and crest-factor behaviour
* Static transfer curves
* Harmonic distortion and H2–H10 structure
* Harmonic behaviour across frequency
* Even/odd harmonic balance
* Harmonic decay and stability
* Linear frequency response
* Phase and group-delay behaviour
* Residual signal characteristics
* Timing, ringing, overshoot and settling
* Intermodulation / signal interaction
* Stereo integrity and crosstalk
* Summing non-additivity
* Boundary behaviour
* Nyquist-aware alias stress testing

Measurements can be performed across configurable plugin parameter positions and input levels.

Input-level grids can currently be run at 1 dB, 3 dB or 6 dB resolution.

## Plugin DNA

PluginDNA does more than store raw measurements.

The measurement engine generates progressively distilled layers of evidence:

**Raw Evidence**
Full analyser CSV outputs containing the underlying measurements.

**Operating Point DNA**
A compact description of behaviour for each tested combination of plugin setting and input level.

**Feature DNA**
Higher-level characteristics derived from the complete measurement set.

This allows the same plugin to be examined both scientifically and practically without discarding the underlying observations.

## Evidence exports

PluginDNA has two export modes.

### Compact Evidence

A JSON evidence file containing:

* test configuration
* parameter definitions
* operating-point behaviour
* harmonic characteristics
* transfer characteristics
* residual characteristics
* alias behaviour
* stereo and summing behaviour
* measured extrema and thresholds
* FeatureDNA and parameter influence information

### Full Evidence

A lossless research package containing the compact evidence file plus the complete analyser CSV datasets.

The principle is simple:

> **Preserve the measurements. Interpretations can improve later.**

This means future analysis methods can be applied to previously measured plugins without necessarily running the plugin again.

## Current direction

PluginDNA is currently being developed toward analysis of complete signal paths.

Planned work includes:

* expandable serial plugin chains
* measurement probes between processing stages
* Stage DNA
* per-stage contribution analysis
* capture of plugin GUI settings for measurement
* visualisation of harmonic, transfer, residual and alias behaviour
* signal-path / console evolution views

One of the primary research uses is investigating channel, group-bus and mix-bus processing as complete systems rather than evaluating processors only in isolation.

## Supported plugin formats

Current macOS build:

* VST3
* Audio Unit

VST2 support is intentionally not included.

## Building

PluginDNA is currently a development build.

Requirements include:

* macOS
* CMake
* C++17-compatible compiler
* JUCE

JUCE is fetched by the CMake project.

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j1
```

Launch the macOS application with:

```bash
open "./build/PluginDNA.app"
```

The command-line measurement target is also built:

```text
plugin_measure_grid_cli
```

## Downloads

Pre-built application downloads will be published through **GitHub Releases**.

Until a release build is available, PluginDNA should be considered development software and built from source.

## Project origins and attribution

PluginDNA began as an extension of the open-source **Plugin Analyser** project by Conceptual Machines / Luca Romagnoli:

https://github.com/Conceptual-Machines/plugin-analyser

Plugin Analyser provided the original JUCE-based plugin-hosting and automated measurement foundation, including grid-based parameter testing, transfer-curve analysis, RMS/peak measurement, harmonic/THD analysis, linear-response measurement and CSV export.

PluginDNA has subsequently expanded that foundation toward a broader DSP evidence system, including additional analyser families, operating-point fingerprints, residual analysis, harmonic fingerprinting, stereo and summing measurements, alias stress testing, richer evidence exports and planned multi-stage signal-path analysis.

Many thanks to the Plugin Analyser project for making the original work available as open source.

## Licence

PluginDNA contains work derived from **Plugin Analyser**, which is distributed under the MIT License.

The original Plugin Analyser copyright and licence notice must remain with the portions of the project derived from that work.

See `LICENSE` and any included third-party attribution notices for details.

## Status

PluginDNA is experimental research software.

Measurements and derived interpretations should be treated as evidence to investigate, not as absolute judgements about audio quality.

A processor producing more harmonics, aliasing, compression or non-additivity is not automatically better or worse. The purpose of PluginDNA is to show **what changed, where it changed, and how that behaviour develops**.
