# Verification report — version 1.0.0

The focused release preserves the frozen sequential implementation and all 230
Round2 comparisons. The following local validation receipts are retained. This
report distinguishes fresh checks, historical evidence and configured CI.

| Check | Result and retained evidence |
|---|---|
| Strict portable C archive build | Passed; `evidence/validation/release_checks.json` |
| Original-file and kernel integrity | Passed; 47 byte-preserved mapped files at assembly; `tests/integrity.py` |
| Native/scalar direct Python integer and Fraction oracle | Passed; 2,044 entries per implementation; `evidence/validation/sequential_results.json` |
| Rejection, dimensions, strides, grid caps and arithmetic witnesses | Passed in the portable suite; details in the same receipt |
| Original six-method Accelerate benchmark build | Passed without timing; fresh check receipt |
| Independent full-entry macOS vendor/candidate oracle | Passed; 5,050,624 comparisons, zero exact errors, missing writes, nonfinite outputs, input changes or padding violations; `evidence/validation/algebra_lattice_full_entry_results.json` |
| All-230 statistical reproduction | Passed; generated per-run CSV byte-identical to historical CSV; 150 paired wins and 150 dual interval wins; `evidence/validation/release_checks.json` |
| Pinned Lean build and theorem/axiom inspection | See original build and independent audit receipts in `lean/evidence/`, with scope in `docs/LEAN_VERIFICATION.md` and `docs/LEAN_INDEPENDENT_AUDIT.md` |
| Independent fresh Lean project replay | Passed from byte-identical sources with no copied project proof artifacts; all five modules rebuilt, all 53 printed types/closures byte-identical; dependency/runtime cache only; `lean/evidence/independent/` |
| Clean isolated source-release validation | Passed portable strict build, both 2,044-entry native/scalar oracle checks, 47-file integrity check, all-230 analysis, and byte-identical historical CSV reproduction; `evidence/validation/independent_clean/receipt.json` |
| New 20-case specialized recorded-data reproduction | Passed with no specialized executable present, all frozen source hashes and job receipts checked, CSV byte-identical; `evidence/validation/independent_clean/specialized-integrity-replay.json` |
| Strict Lean axiom CI guard | All 53 expected closures accepted; corrupted missing/duplicate/malformed/custom-axiom/`sorryAx` receipts rejected; `scripts/check_lean_axioms.py` and independent clean negative-case receipt |
| Metadata and maintained Markdown relative links | Version, date, author and title match; all 17 inspected local links resolve; PDF outputs, DOI/deposition and historical Markdown excluded; `evidence/validation/independent_clean/metadata-links.json` |
| Manuscript review and compilation | See `docs/MANUSCRIPT_AUDIT.md` and manuscript build log; compilation status belongs to the final manuscript build |
| Hosted continuous integration | Workflow supplied in `.github/workflows/verify.yml`; adding the workflow is not evidence of a hosted run |

The original performance campaign contains 230 historical Round2 comparisons.
During release preparation a separate, explicitly frozen 20-case specialized
comparator campaign was newly timed; its receipts are under
`data/specialized_2026-10-19/`. These are distinct experiments. The subsequent
packaging validation and isolated checks performed no performance timing: they
rebuilt correctness targets and recomputed statistics from the retained receipts.
The new specialized analysis succeeded without its benchmark executable present;
historical binary execution provenance relies on the matching frozen manifest,
environment and per-job receipts, while all frozen source hashes were checked
against available source bytes. Compilation, numerical correctness and successful
analysis do not extend measurements to a different machine.

To repeat local checks, run `scripts/build.sh`, `scripts/test.sh`,
`python3 analysis/reproduce.py`, and on macOS `scripts/test_full_oracle.sh`.
Lean and paper commands are in the root README. New local correctness receipts
replace files only under `evidence/validation/`; original evidence stays under
`evidence/original/`. Build products and fetched dependency caches are disposable
and excluded from the source release.

The formal proof verifies mathematical algebra, bounds and rounded reconstruction
under an explicit primitive rounding model and binary64 bridge. It does not prove
that a compiler or CPU implements that model, or verify the C memory behavior.
The finite independent oracles and source inspection address that separate
implementation boundary without claiming a fully verified executable.

Final manuscript: Tectonic 0.17.0 and BibTeX successfully produced the eight-page
`paper.pdf`. All eight rendered pages were visually inspected, with enlarged
inspection of the plot and formalization/comparator tables. No unresolved
references, citations, missing characters or overfull boxes remain. Two underfull
line-spacing warnings in the Lean table are benign. See
`evidence/validation/manuscript.json` and `evidence/latex-build.log`.

Publication completed: GitHub release `v1.0.0`, frozen commit `b52224b446e2b5f164791e2556774bc279b61927`, and published version DOI [10.5281/zenodo.23475920](https://doi.org/10.5281/zenodo.23475920). Hosted CI passed all three jobs. The downloaded Zenodo ZIP matches all 150 frozen files and its displayed checksum; [publication verification](docs/PUBLICATION.md) records the immutable version and subsequent citation metadata update.
