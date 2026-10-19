# Frozen experiment protocol (2026-10-19, before candidate timing)

Target: real AA^T, FP64 inputs and outputs, preallocated row-major lower
triangle. Every method reads the same A and writes every lower entry. beta=0.
No approximate rank truncation and no hidden preprocessing supplied to candidate.
All per-call validation, packing and temporary allocation are timed. Input/output
allocation, generation and reference calculation are outside the timed interval.

Hardware: Apple M4, macOS 27.2 (26B5091g), arm64, 10 logical CPUs, 16 GiB.
Compiler: Apple clang 21, -O3, no fast-math, explicit fma permitted.
Baseline: system Accelerate cblas_dsyrk, plus cblas_dgemm diagnostic. Fixed
VECLIB_MAXIMUM_THREADS=1 and BLAS single-thread API where available. Single benchmark
process, no competing research benchmarks. Thread control behavior to be audited.

Independent mechanism regimes selected before timing:
1. Dense skinny: m=128,256,512,1024 and n=1,2,4,8,16,32 (includes failures).
2. Exact sliding-window/Hankel: A_ik=x_(i+k), m=128,256,512 and
   n=128,256,512,1024. Dense validation and generator extraction included.

Accuracy acceptance fixed before timing: for every measured entry,
|Ghat_ij-Gexact_ij| <= 2e-12 ||A_i||_2 ||A_j||_2, with zero rows exact zero.
Nearly orthogonal dot products use the same norm-based absolute requirement;
relative error on arbitrarily small nonzero inner products is not promised.
Adversarial accuracy checks include cancellation, near dependent rows,
row scaling, signal scaling/dynamic ranges up to 2^200, zeros and odd dimensions.
Tests exclude arithmetic overflow/underflow from the rounding-error theorem and
must explicitly reject or fall back when its prerequisites cannot be certified.

Development seeds: 20261019, 20261020. Holdout seeds: 20261101..20261105,
not used to choose kernel or anchor parameters. Accuracy distributions: bounded
pseudorandom, Gaussian-like sum-of-uniforms, nearly constant/periodic, impulses,
alternating signs, and powers of two. Benchmark warmups >=5 calls, >=21 paired
samples in balanced randomized order, each sample batches >=10 ms where feasible.
Report median, IQR and paired speedup bootstrap 95% CI, raw samples and environment.
Final timings must be rerun after implementation freezes; development measurements
are explicitly labeled. Competing search methods get 30 CPU-second exploratory
budgets; exact verification is separate and unbounded by that search allowance.

Novelty: distinguish known algebraic/communication mechanisms from locally authored
algorithm/guard/kernel. Literature absence is not a priority proof.

## Development-to-evaluation freeze extension

Exploratory measurements exposed scalar validation/certificate overhead; they are
not final evidence. In particular, independent correctness work sometimes ran
while development timings were collecting. No such overlap is permitted in final
runs. Development results must never be substituted for final results.

Before final measurement, add exact sparse-shift inputs formed by independently
zeroing one coordinate in each dense window (hence at most two displacement
nonzeros per row). Include m=128,256,512 and n=1024,4096,16384. This tests window-
specific missing/corrected samples; editing only the shared generator would not.

Reserve unseen dimension cases m=192,384,768 with n=3,6,12 for dense evaluation,
and (m,n)=(384,2048),(384,8192),(768,4096) for sparse evaluation. Parameter/tile
choices must be frozen before these cases are run. Use seeds 20261101..20261105.
Final samples: 24 per dense/sparse case (4 methods), 25 per Hankel case
(5 methods), arranged in randomly permuted balanced cyclic method orders.

Primary output contract remains lower triangle, as SYRK. Also time full symmetric
materialization by applying the same in-place mirror to each method inside the
clocked call, with GRAM_FULL=1, and report those results separately.
Final requirements do not change: same FP64 accuracy threshold, same library,
thread model, input/output layout, inclusive preprocessing and allocations.

Final harness freeze correction: final Hankel comparison uses four methods:
SYRK, GEMM diagnostic, fixed classical sliding recurrence, and final SIMD guard.
Thus all final samples=24, balanced over4positions. Sparse uses SYRK, GEMM,
initial scalar certificate, and final fused sparse certificate. Frozen hashes and
actual hardware/power settings are in results/environment.json.

## Discovered harness issues and corrected final provenance

The v1dense distribution generator ignored dist for non-Hankel matrices. All30
finaldense evaluationcases were rerun with the correctedv2generator. Original
rows are retained only in results/pre_correction_dense_evaluation.jsonl; current
final_evaluation_lower_grid.jsonl contains correctedrows. Uniformmaincases and
allHankel/sparse cases were unaffected. The impulse distribution is deterministic:
five seeds vary runordering, not the impulsevalues. Do not call these fiveunique
inputdraws.

The v1/v2fullmirror made stridedstores and redundantly mirrored alreadyfullGEMM.
Finalfullcases were entirely rerun using a shared32x32blockedmirror for triangular
methods and nativefullGEMM without addedmirroring. Originalfulltimings are in
results/pre_blocked_mirror_full_grid.jsonl and are excluded from final summaries.

Rawfinalrows include harness_revision:1originaluniformlower/structuredevaluation,
2correcteddenseevaluation,3blockedfull. Frozen C sources/binaries benchmark_v1,
benchmark_v2,benchmark_v3 preserve those versions. Currentv4only broadens the
untimedfullpreflight to checkupperentries too; all57fullcases passed that check.
The implementationkernels and their hash-linked proofs never changed in these
harnessrepairs. Current default commands regeneratecorrectinputs/completechecks.

The runner now explicitly setsGRAM_FULL0/1 rather than inheriting a possibly
conflicting environment value. The summary groups evaluationbyoutputcontract.
Confidence intervals use5000bootstrapresamples of contiguousfour-trialbalanced
blocks, not falsely independent individualratios. They remain conditional
within-run descriptions and do not guarantee cross-run/hardwarecoverage.

## Round 2: exact-lattice specialization, frozen before any new timing

The existing general, wide-range certified family and all its receipts remain in
scope. Add an exact-arithmetic branch for validated A=2^e X with integer |X|≤32768,
not approximate quantization. Every value and sparse shift residual must be
validated from the same dense A on each timed call. A failed exact-lattice test
must reject or use the original guarded path; it must not change input values.
The existing accuracy threshold2e-12 stays fixed. New branch aims at exactly
represented output (errorzero), by a proved bound on every intermediate.

Keep Apple M4/Accelerate FP64/single-thread settings, row-major layout and the
inclusive per-call clock boundary. Include full matrix output with a shared
32x32mirror for triangular methods and exempt nativefullGEMM. Competitors:
AccelerateSYRK, nativeGEMM, existing certified floating path, exact-lattice path.
Minimum5warmups,24balancedpairedsamples,10msinitialbatchcalibration.

Development only: (m,n)=(256,1024),(512,4096), q=0and2, e=-15,
seed20261021. Uniform signed16bit values. Two independent structure mechanisms:
(a) Hankel with independent per-windowzero replacement (q<=2), and
(b) accumulated sparse innovations: shift previousrow, insert a new tail, then
replace q independently selected interior coordinates by signed16bitvalues.
Case(b) is not assumed Hankel plus onlyO(q) edits per entire row.

Reserved dimensions: (192,768),(384,2048),(768,8192), q=1,2,4,
e=0and-15; evaluationseeds20261201..20261205. Keep cases separate from development.
Final repeated distributionchecks at a selected developmentpromising shape also
use the reservedseeds. Signals: bounded signed16bitrandom, nearconstantinteger,
alternating/impulses, saturationextremes and sparseinnovations. Deterministic
signals repeatedwithdifferentorderseeds are explicitly marked as repeatedinputs.
No implementationparameter may depend on heldoutperformance inspection.

This is an additional exactness mechanism within the original goal, not a
replacement for the existing general-real correctness/numerical obligations.
Prior-art distinction must include known exactfloating arithmetic/BLAS use;
integerexactness alone is not a noveltyclaim. No new timing has been run when
this extension was written.

Before the first round2timing, add the coordinate-substitution exactness variant
alongside the grouped exact implementation. Both are compared in the same case,
so round2uses five methods and25samples (five balancedblocks of five), not24.
Neither candidate is selected percase afterseeingfinal results. The grouped
implementation is an explicit comparator; the substitution schedule aims to
match the classical nB²<=2^53 exact-dot sufficientcondition independentlyofq.
No originalkernel has been changed. Substitutiononlychanges sparsecorrection
order and admissionbound; the real polynomialtarget remains AA^T.

Analyticalrepair before timing the sequentialvariant: the same nB² bound can be
obtained without supportmerge by updating the left coefficient first with old
right, then the right coefficient with fullynewleft. Add this distinct sequential
schedule to the two existing exact comparators. Round2final timing will therefore
have six methods and30samples, balanced in blocks of six; bootstrap uses the same
six-trial blocks. The first four developmentcases (five-method25samples) remain
archived separately and cannot be relabeled as final. The parameter and reserved
input lists above stay fixed. The sequentialschedule and batchedvalidation are
repairs motivated by its invariant and the developmentcosts, not heldoutdata.

Round2 final job manifest is now frozen before evaluation: 180 uniform
accumulated-innovation runs = three reserved shapes × q{1,2,4} × e{0,-15}
× five reserved seeds × lower/full output; plus 50 distribution runs = two
families (corrected windows and accumulated innovations) at512×4096,q2,e-15
× five listed distributions × five reserved seeds, full output. Method13 is the
preselected candidate throughout; no per-case choice among methods11/12/13.
The randomized global order is fixed in run_lattice_benchmarks.py; the complete
230-job list and hashes are in results/lattice_frozen_manifest.json. The source
kernel remains7550990a7f6dd82fbd0e58580c9a6639eef7ffc07d841b56ae9cf189b394bde5.

Untimed preflight was independently repaired before final timing: reject
nonfinite outputs explicitly; poison the reused output before each method's
warmups; scan every required entry for finiteness; retain independent exact int64
sample checks and full-output symmetry checks. Record actual residual counts.
These repairs change no timed call or kernel. Final harness SHA256 is
5f89c6a6bf16a6f4dd53e8592de79f3e550a92219b7cceff72a02dfb97f4c1e5.
Bootstrap5000 resamples of the five contiguous six-trial balanced blocks. Its
intervals describe within-run variability; they do not imply cross-machine or
long-term scheduling coverage. Every source/binary hash is checked between runs.
