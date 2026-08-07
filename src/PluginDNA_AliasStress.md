# PluginDNA Alias Stress Test

A dedicated Nyquist-aware alias validation pass.

## Stimulus

Ten sequential sine tones, each held for 0.5 seconds:

500, 1000, 2000, 5000, 8000, 10000, 12000, 15000, 18000 and 20000 Hz.

## Method

For every parameter and input-gain run, PluginDNA:

1. Measures the final FFT window of each tone segment.
2. Classifies legal H2-H10 components below Nyquist.
3. Predicts the mirrored frequency of every H2-H10 component above Nyquist.
4. Measures energy at each unambiguous predicted fold location.
5. Exports per-frequency alias power and a compact feature summary.

Fold locations that overlap a legitimate harmonic are excluded because their energy cannot be separated reliably from the valid component.

## Export

`grid_alias_stress.csv`

Key outputs include alias power, first detected alias fundamental, highest clean fundamental, strongest alias point, detected fold count and ladder integrity.

The default detection threshold is -90 dB relative to the fundamental.
