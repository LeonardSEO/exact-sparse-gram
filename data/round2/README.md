# Frozen Round2 measurements

These three files are byte-identical copies from `gram_research/results/`.
The JSONL holds all 230 runs, each with 30 paired observations of six methods.
The pre-timing manifest fixes the jobs, method selection, sources and binary.
The completed environment receipt adds the actual host, power, thermal and
completion metadata. Their creation timestamps differ; all other shared fields
match. See `../../docs/PROVENANCE.md` for provenance and statistical limits.

Reproduce summaries from the release root with `python3 analysis/reproduce.py`.
This reads the frozen measurements and writes only `analysis/generated/`.
