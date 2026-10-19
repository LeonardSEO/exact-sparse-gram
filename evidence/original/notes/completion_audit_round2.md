# Independent completion audit, round 2

Date: 2026-10-19. Scope: source and existing evidence only. No benchmark,
compiler, numerical verification suite, or competing CPU workload was run.
Small standard-library reads recomputed raw sample counts, medians and SHA-256
digests. `proof-audit` and `contribution-review` supplied the review criteria.
This audit does not infer completion from the manifest or previous audit votes.

## Verdict

The delivered sparse-shift construction has a complete exact-real argument, a
credible source-matched universal numerical guarantee on its explicit domain,
and a concrete inclusive practical win against both measured vendor methods.
I found no new algebraic or numerical defect in the reconstructed sparse path.
The broader objective of a **new** provably correct, practically faster family
is not yet fully supported as a research originality claim. The outstanding
task is a precise, bounded contribution assessment against the closest existing
results, rather than proof of absolute first-in-history priority.

Correctness, local authorship, measured utility, and originality are separate
facts. The first three are substantially evidenced here. Newly written source
and the absence of an identical search result do not establish the fourth.
Conversely, established individual components do not refute novelty of a new
quantitative guarantee or a new integrated construction.

## Requirement-by-requirement evidence

| Requirement | Evidence actually inspected | Assessment |
| --- | --- | --- |
| A general family rather than one small identity | `notes/algebra.md:8-75`; `REPORT.md:63-123`; arbitrary row-sparse displacement, no dimension divisibility, optional direct anchors, zero/one-dimensional edges | Established exact-real construction for arbitrary dimensions. The executable certified family deliberately restricts binary64 range, n, and q. This is not an accelerated algorithm for every dense matrix. |
| Exact real proof, with valid ordered/block extension | Reconstructed `U_i U_j* + D_i A_j* + U_i D_j* = A_i A_j*` without commuting factors; remove the old head and add the new tail; induction on row; direct anchors compute the same target | Proof established for this recurrence and its compatible block version. The block version is a mathematical extension, not a delivered recursive block performance result. |
| Structure is validated, not provided free | `src/numerical_sparse_fast.c:122-194`; vector/tail range checks, equality tests, TwoDiff extraction, per-row count check before every residual store | All original dense input coordinates are scanned and actual residuals extracted. Each residual is represented by both hi and lo. The API rejects excess q or invalid domain. |
| Predetermined accuracy of the original represented data | `PROTOCOL.md:17-23`; `src/numerical_sparse_fast.c:195-249`; `notes/numerics.md:53-151`; `notes/anchor_totality.md:21-151` | The row-norm-scaled absolute target is fixed at 2e-12. Eight FMA streams plus three tree levels give depth <=8195; anchor and norm enclosures use outward/downward rounding. The recurrence error has coefficient one on incoming error, with each FMA contributing <=rho times its rounded magnitude. The exact TwoDiff limbs preserve the target matrix. |
| Accuracy is total on the supported family, not only conditional successful output | `notes/anchor_totality.md:46-151`; actual helper operations at `src/numerical_sparse_fast.c:23-42,195-247` | Reconstructed the loose constants: stored fresh-anchor certificate <1.01e-12 times the norm product, implemented target >1.98e-12 times it. Thus one fresh repair always suffices at fixed tolerance after validation/allocation. Zero rows pass exactly. The input/product lattices exclude underflow in recurrence arithmetic. No new missing totality lemma identified. |
| Independent executable evidence | `verification/numerical_sparse_fast_verify.py:20-90`; saved `verification/numerical_sparse_fast_results.json`; `results/algebra_exact.json`, `results/kernel_exact.json`, `results/final_verification.log` | The sparse oracle forms exact Fraction products and compares squared error against squared exact norms. Its saved 25 cases/1308 entries include short tails, q=64, nonzero low limbs, and 21 repairs. The ordered algebra receipt reports 2112 identities, 1200 rational cases and 1440 compiled integer calls. These are inspected receipts, not fresh executions and not replacements for the universal proof. |
| Strong same-operation comparator including SYRK | `src/benchmark.c:27-57,88-113,129-139`; `build.sh:4-7`; preserved v3 harness | SYRK is row-major lower/no-transpose/alpha=1/beta=0. All methods use identical A and preallocated G. Full mode includes the same blocked mirror for triangular methods; native GEMM is correctly exempt. Fixed-tolerance validation, packing, allocation/free, seed/norm work and repairs occur inside timed calls. Vendor single-thread mode is selected and checked. |
| Reproducible practical advantage | Read all three `results/final_*_grid.jsonl` files and recomputed headline medians; inspected `run_benchmarks.py`, `summarize.py`, `results/environment.json` | 234 comparison records, each with 24 samples, four preflight methods and recorded BLAS mode 1: 57 lower v1, 57 full v3, 120 evaluation records (30 corrected v2/90 v1). These contain 230 distinct mode/shape/seed/distribution/output tuples because main and evaluation jobs overlap; they are 234 runs, not 234 unique inputs. |
| Concrete full-output win over both measured vendor options | The 768x4096 sparse raw full record and `results/final_summary.json` | Recomputed medians: SYRK 7.2983229766 ms; native GEMM 7.6152085094 ms; candidate 5.2285259881 ms. Recomputed median paired SYRK/candidate ratio 1.3536322468. Stored four-trial block-bootstrap interval is [1.3089332811,1.4353316973]. This is the strongest directly supported practical result. |
| Honest performance regime and retained failures | Raw full records for sparse384x2048, dense1024x1 and dense512x32; lower evaluation grids; `REPORT.md:28-54,344-353` | Sparse384x2048 has paired full speedup 0.7913; dense1024x1 loses to GEMM despite beating SYRK; n=32 dense cases lose. Full headline data are bounded-uniform. The six-distribution/five-order-seed sweeps are lower output; impulses repeat inputs. No full-output all-distribution, all-input or cross-hardware guarantee is established or needed for the bounded reported win. |
| Numerical, memory and temporary analysis | `REPORT.md:125-249`; `src/numerical_sparse_fast.c:139-148,215-252`; `src/joint_kernel.c:29-42,47-82,136-153` | Scratch capacity is 24(m-1)min(max_q,n-1)+24m+8 bytes on this ABI, giving 55,256 bytes at m=768,q=2. Input/output allocations are separate and shared. Work includes O(mn) extraction and O(Rn) repairs; worst case returns to dense work. Traffic is symbolic, not a falsely measured DRAM/RSS improvement. Dense packing and its full allocation are included. |
| Source/evidence provenance | Rehashed `src/numerical_sparse_fast.c`, `src/communication_guard_fast.c`, `src/joint_kernel.c`, `src/benchmark_v2.c`, `src/benchmark_v3.c`, `src/benchmark.c` | All match their frozen kernel or relevant harness-revision hashes in `results/environment.json`. Corrected generator/mirror results are separated from archived earlier timings. The 57-entry full-output receipt records successful upper-and-lower checks. |
| Existing knowledge versus the original contribution | `notes/literature.md:184-229,250-354`; `REPORT.md:309-342`; `CHECKPOINT.md:21-34` | Known sliding recurrences, TwoDiff, restarts, arithmetic outward bounds and tiling are credited. The q-row-sparse/full-displacement-rank distinction is valid. But the closest broader residual/displacement subsumption comparison is explicitly unfinished (`literature.md:225-229`). The current record supports an independently derived integration and measured artifact; it does not yet establish which exact new theorem or guarantee survives the nearest-prior comparison. |

The source hashes most relevant to this verdict are:

```
7a706721832076901a1845031e6e5ea2e8620ad821b40ca24eddc02b0746b5d0  src/numerical_sparse_fast.c
40b33b64d12eac5fa5a3ca2e2acd00b317abd2c1c790b8f36c20430edcd912ae  src/communication_guard_fast.c
e3aa02beff68e108dc750058e54e7c2b2f495638489544c3a65c231f192e4306  src/joint_kernel.c
21764ac40f83c93158ed2407c44c9c1be22e6ccd94346b34967814c16893afab  src/benchmark.c
```

## Strongest gap and concrete corrective action

`notes/literature.md:225-229` explicitly leaves the relation to more general
residual/displacement frameworks open. Its search log at lines 325-347 reports
no identical complete integration, and names priority unresolved. That is
honest reporting, but it is not yet a positive argument identifying the research
delta requested by “new.” The non-low-rank example only excludes the inference
“row sparsity implies low rank”; it does not exclude a general framework that
can already exploit sparse residual application.

Do the following bounded comparison before either declaring the original
research objective achieved or deciding that further invention is necessary:

1. Freeze a single central statement: arbitrary-support q-row-sparse shift
   residuals, including full displacement rank, exact represented-input Gram
   output, validated dense-input extraction, O(mn+m^2(q+1)+Rn) work,
   O(mq+m) scratch, and the explicit norm-scaled certificate/totality guarantee.
2. For the nearest primary sparse-residual/displacement methods and corrected or
   masked sliding-window methods, record the exact old hypotheses, operations,
   work/storage bounds, and accuracy guarantees. Test whether their residual
   application cost specializes directly to this statement. A difference in
   vocabulary or absence of the literal source code is insufficient.
3. Identify at least one precise, proved added capability or guarantee, or
   explicitly classify the recurrence as a specialization/reformulation and
   assess the numerical/scheduling integration as the candidate contribution.
   A new integrated algorithm can qualify without each primitive being new;
   the evidence must explain the actual nontrivial delta.
4. If the closest known method already supplies the entire target claim, retain
   the verified kernels as baselines and develop a materially different
   construction or guarantee. If it does not, state a scoped evidence-based
   originality conclusion with citations and remaining search limits. Do not
   demand a logically impossible certificate that no earlier work anywhere exists.

This is a finite theorem-level review obligation. It is not a requirement for
peer acceptance, a different model vote, Lean formalization, or absolute priority
proof. It preserves the user's original request for new work while avoiding an
unfinishable standard of historical certainty.

## Small evidence-maintenance correction

`notes/audit_benchmark.md:173-192` describes the earlier individual-ratio
bootstrap and its old `summarize.py` digest `2c6e8c6a...`. Current
`summarize.py:8-10` correctly resamples contiguous four-trial balanced blocks,
and has SHA-256
`7c70e294744ec236559ab3e3e36278ee264655c8773beb7a48cee99c310bc80a`.
`REPORT.md` and `PROTOCOL.md` already describe the new analysis. Add a dated
audit addendum for that change and hash, preserving the old audit's chronology.
This does not require new timings and does not undermine the raw practical win.
