# PluginDNA Frozen Five — implemented

1. Open Plugin button for primary and optional serial processor.
2. Stereo Suite: L-only, R-only, Mid and Side tests; exports `grid_stereo_dna.csv`.
3. Summing Suite: deterministic A, B and A+B test; exports `grid_summing_dna.csv` and non-additivity metrics.
4. Serial Stage Analysis: optional second VST3 processed after the primary plugin on every pass. Both plugin states are restored into fresh pass instances.
5. Feature Extraction Audit v5: existing static sensors plus compact Stereo and Summing descriptors. Evidence schema 11.

The Summing Suite compares `P(A+B)` with `P(A)+P(B)` using phase-reset deterministic waveforms, so the residual represents nonlinear summing interaction rather than mismatched test material.
