# Six-method Round2 benchmark

`src/` preserves the project-local original drivers and their link dependencies.
Method13 is the released sequential kernel; its copy is byte-identical to the
public source. Methods0/4 are vendor SYRK/GEMM, 9 the earlier certified member,
11 grouped lattice, and 12 merged-support substitution. Method12 is already a
specialized sparse-shift comparator with the same call/output contract.

Build on macOS: `scripts/build.sh --benchmarks` from the release root.
Prepare a fresh run manifest: `python3 benchmarks/round2/run.py --dry-run`.
Run all 230 jobs explicitly: `python3 benchmarks/round2/run.py --grid all`.
Run timing only without competing compute. The adapted runner never overwrites
frozen data; fresh outputs go to `data/reproduced/`. Analyze a complete fresh run with
`python3 analysis/reproduce.py --data-dir data/reproduced --output-dir analysis/reproduced`.
The analysis still requires exactly 230 jobs and never silently mixes runs.

The recorded run used Apple Accelerate, one thread, warm repeated operands and
full per-call validation/extraction/allocation/materialization. Portability of
the public kernel does not imply that these performance results transfer to
another architecture. See `docs/PROVENANCE.md` and the paper for limitations.
