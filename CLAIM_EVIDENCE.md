# Claim–evidence map

Each claim below has a specific proof, implementation or measurement scope.

| Claim | Primary evidence | Boundary |
|---|---|---|
| Sequential left-then-right factor replacement returns the real Gram dot | `paper.tex`, `lean/GramExact/Core.lean` | Recurrence and reconstruction over the stated mathematical objects |
| Every correction prefix is a sum of at most n bounded factor products | `paper.tex`, `lean/GramExact/Bounds.lean` | Actual old/new input factors; includes overlapping supports |
| Bounded integer-grid intermediates are exactly reconstructed under the rounding model | `lean/GramExact/Arithmetic.lean`, `lean/GramExact/Bridge.lean`, `docs/LEAN_VERIFICATION.md` | Explicit primitive-rounding/binary64 bridge assumptions; C and hardware remain external |
| No recurrence-length drift for admitted bounded-grid inputs | Paper exactness proof; rounded reconstruction theorem; independent proof audit | Requires IEEE arithmetic assumptions and the admitted domain |
| Gram work is mn + 2(P−m) + mQ FMAs | Paper work count; `lean/GramExact/Bounds.lean` | Counts Gram arithmetic, not all validation/control operations or physical traffic |
| Delivered source validates actual dense input and returns exact lower output on success | `src/lattice_sequential.c`, `include/gram_exact.h`, `tests/sequential_verify.py` | Implementation theorem assumptions and finite tests; no verified compiler/C memory proof |
| Frozen kernel is unchanged | `tests/integrity.py`, `evidence/original/source_mapping.json` | SHA-256 `7550990a7f6dd82fbd0e58580c9a6639eef7ffc07d841b56ae9cf189b394bde5` |
| Portable native/scalar oracle checks pass | `evidence/validation/sequential_results.json` | 2,044 lower entries per implementation plus rejection/boundary and finite invariant checks |
| Representative full-entry oracle comparisons are all exact | `evidence/validation/algebra_lattice_full_entry_results.json`; historical receipt under `evidence/original/verification/` | Four cases in two layouts; sequential 1,263,744, SYRK 1,263,744, GEMM 2,523,136 entries; 5,050,624 total |
| Headline paired SYRK speedup is 2.169 [2.164,2.175] | `data/round2/lattice_final_all.jsonl`, `analysis/generated/lattice_selected_measurements.md` | One specified full-output 768×8192 q2 e−15 uniform innovation case, seed20261201, Apple M4 Accelerate |
| All 230 comparisons and all six method series are retained | `data/round2/lattice_frozen_manifest.json`, `lattice_environment_all.json`, `lattice_final_all.jsonl` | Frozen jobs, 30 paired samples each; reused/scaled/deterministic inputs are disclosed |
| 150 paired wins against both vendors; all 60 small-shape comparisons lose | `analysis/generated/lattice_statistics.json`, `lattice_summary.csv` | Fixed evaluation grid; candidate13 is preselected, not chosen per case |
| The analysis reproduces the historical per-run CSV exactly | `analysis/reproduce.py`, `evidence/validation/release_checks.json` | All 230 rows byte-identical; deterministic five-block bootstrap; within-run intervals |
| Same-contract specialized comparator is retained | Method12 in `benchmarks/round2/src/`, `docs/SPECIALIZED_COMPARATOR.md` | Merged-support sparse substitution; additional specialized comparators have their own documented scope |
| Lean formalization built and assumptions were inspected | `lean/evidence/`, `docs/LEAN_VERIFICATION.md`, `docs/LEAN_INDEPENDENT_AUDIT.md` | Pinned Lean/mathlib; distinguishes checked propositions from external machine correspondence |

The artifact makes no unrestricted dense speedup, universal hardware/runtime,
new matrix-multiplication exponent, or absolute historical-priority claim. Timing,
finite correctness tests and universal mathematical statements are separate forms
of evidence. The original records remain byte-preserved; fresh outputs never
rewrite historical receipts.
