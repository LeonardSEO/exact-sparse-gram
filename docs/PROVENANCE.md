# Source, data and validation provenance

Leonard van Hemert confirmed authorship and licensing for the project-local code
(MIT) and retained data (CC BY 4.0). No external repository, third-party source,
FHE work, deletion work, flip-track code, cache, or historical executable is
included in this focused artifact. The code retains original comments; source
attribution and a byte-level copy mapping are in
`evidence/original/source_mapping.json`.

## Frozen implementation and public interface

`src/lattice_sequential.c` is an unchanged copy of
`gram_research/src/lattice_sequential.c`, SHA-256
`7550990a7f6dd82fbd0e58580c9a6639eef7ffc07d841b56ae9cf189b394bde5`.
The focused public interface is `include/gram_exact.h`; the released core has no
vendor BLAS dependency. `scripts/build.sh` produces `build/libgram_exact.a` with
strict floating-point flags. IEEE binary64 and correctly rounded operations,
including FMA, are required. This release does not change the frozen algorithm.

## Retained experiment

`data/round2/lattice_final_all.jsonl` contains all 230 frozen Round2 comparisons,
including the losing cases. `lattice_environment_all.json` is the completed
experiment receipt; `lattice_frozen_manifest.json` is the pre-timing manifest.
They retain original source hashes, original binary hashes, original local
paths, compiler/environment metadata and power/thermal reports as historical
evidence. Embedded historical paths are not runtime dependencies.

The six methods are Accelerate DSYRK (0), native Accelerate DGEMM (4), earlier
certified floating (9), grouped lattice (11), merged-support substitution (12),
and the preselected sequential lattice (13). The original benchmark drivers and
all project-local sources needed to link them are byte-preserved under
`benchmarks/round2/src/`. Driver `benchmark.c` has unused historical entry points,
so their local dependencies are retained for faithful linkage. They are support
code, not additional claims of this release.

The merged-support method12 already supplies a specialized structured comparator
under the same dense-input/validation/materialization contract. A separate frozen
20-run H+E specialization was subsequently measured and is retained under
`data/specialized_2026-10-19/`; see [the focused comparator](SPECIALIZED_COMPARATOR.md).
It is not pooled with or substituted for the original 230 measurements. No q=0
STOMP timing or candidate selection was performed. q=0 is covered by correctness
tests; the recorded 230 jobs use q>0.

To build the original six-method benchmark on macOS, run
`scripts/build.sh --benchmarks`. The adapted runner is
`benchmarks/round2/run.py`; `--dry-run` prepares the manifest, and `--grid all`
starts the full timing campaign only when explicitly invoked. The runner also
writes a fresh pre-timing manifest. Analyze fresh runs with
`python3 analysis/reproduce.py --data-dir data/reproduced --output-dir analysis/reproduced`.
Its fresh outputs
are `data/reproduced/lattice_final_all.jsonl` and
`data/reproduced/lattice_environment_all.json`. It refuses to overwrite those
outputs. Run timing without competing work. No performance timing is part of
build, tests or analysis. Original sources and measurements stay frozen.

## Independent numerical tests

`scripts/test.sh` verifies the source/data/receipt hashes, builds native and forced
scalar libraries, and checks the sequential implementation against direct Python
integer dot products and exact Fraction scaling. No recurrence is used by that
oracle. Coverage includes vector tails, all tested single-defect lanes, extreme
exponents, coefficient caps, overlapping supports, q=0/q=64, empty boundaries,
padded strides, cancellation, nonfinite/off-grid inputs and rejection behavior.
The test also verifies the actual-FMA grouped-order counterexample, 6,561 finite
invariant states, an ordered polynomial identity and 371 CSR group prefixes.
`tests/sequential_verify.py` is adapted from the original verifier only to place
support sources, build products and fresh receipts in release-relative folders.
Its successful fresh receipt is `evidence/validation/sequential_results.json`.
The suite checks 2,044 lower entries per implementation. This finite coverage
supplements the proof; it does not establish universal correctness by itself.

`scripts/test_full_oracle.sh` is optional and requires macOS Accelerate. Its
independent integer-first C generator/oracle checks every requested output entry
for four large cases, each in contiguous and padded layouts, against sequential,
DSYRK and DGEMM calls. The adapted wrapper writes binaries to `build/` and fresh
receipts to `evidence/validation/algebra_lattice_full_entry_*`. On other platforms
it reports SKIP; it never claims vendor-oracle validation there. Historical
receipts under `evidence/original/verification/` are unmodified and distinct from
these fresh results.

## Analysis reproduction

From the release root run `python3 analysis/reproduce.py`. It uses only the Python
standard library. It checks that every non-timestamp field in the frozen manifest agrees with
the completed environment receipt, verifies frozen comparator source hashes,
and rejects missing, duplicate, mismatched or incomplete jobs, wrong method/sample
sets, and nonfinite/nonpositive measurements. Exactly 230 jobs with 30 samples
and six specified methods are required. Results go to:

- `analysis/generated/lattice_summary.csv`: all per-run times and paired ratios.
- `analysis/generated/lattice_selected_measurements.md`: six headline examples,
  their block intervals, every group including losses and aggregate win counts.
- `analysis/generated/lattice_statistics.json`: aggregate counts and group ranges.

The estimator is the median of paired per-trial ratios. Intervals use the
original deterministic seed and 5,000 resamples of five contiguous six-trial
blocks (sorted draws 124 and 4874). They are within-run conditional intervals,
not simultaneous or cross-hardware coverage. Native GEMM is included throughout;
no per-case lattice candidate is selected. The original summarizer and original
summary outputs remain in `evidence/original/` for comparison.

The preserved pre-timing manifest and execution manifest differ only in
`created_utc`: freezing and execution occurred at different times. All protocol,
jobs, source, binary and runner fields match; analysis checks those fields.

Fresh validation builds ran with Apple clang on the current macOS host. All
validation outputs describe correctness or data analysis, never new performance
measurements. Build products are disposable and excluded from the release.
