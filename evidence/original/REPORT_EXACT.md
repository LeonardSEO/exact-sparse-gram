# Exact Gram computation by sequential sparse factor updates

Research date: 19 October 2026. This report describes the second construction
round. The broader real-input certificate implementation and its earlier
measurements remain in [REPORT.md](REPORT.md). Original files in the parent
directory have been preserved.

## Construction and claim boundary

The family computes `G=AA^T` for matrices whose successive rows differ from a
one-coordinate shift at a small number of positions. The new member uses a
specific order of factor replacement. For validated bounded integer-grid inputs,
that order makes **every Gram anchor, residual and update intermediate exactly
representable in binary64**, independently of row count or recurrence length. No
rounding onto a grid, lower-precision multiplication, or approximate target is
used. The implementation reads and validates the actual dense input on every call.

The construction applies to arbitrary dimensions mathematically. The delivered
member supports `n<=65536`, at most 64 defects per shifted row, and inputs

```text
A[i,k] = 2^e X[i,k],     X[i,k] integer,
|X[i,k]| <= B_n = floor_sqrt(floor(2^53/n)),
-400 <= e <= 200.
```

These are real matrices, with a real, characteristic-zero correctness proof.
The integer grid supplies a numerical representability invariant, rather than
a finite-field identity. Bounded fixed-point/sample data converted exactly to
double is a relevant input class. The benchmarks use signed 16-bit coefficients;
no claim about a particular deployed sensor or audio workload is made.

The algebraic displacement product rule and exact floating-point integer
arithmetic are established knowledge. The contribution constructed here is the
particular two-loop schedule, its complete intermediate invariant, admission
guarantee, and validated implementation. This is a modest algorithmic and
numerical result, with no claim of a new matrix-multiplication exponent or an
absolute historical first. A finite theorem-level prior comparison appears below.

All 230 frozen comparisons completed. The candidate won against both vendor
routines in 150 runs, including every tested 512×4096 and 768×8192 case.
All 60 comparisons at 192×768 lost. The result is a measured structured regime,
not a universal speedup. [All selected/grouped measurements](results/lattice_selected_measurements.md),
[all per-run statistics](results/lattice_summary.csv), and
[raw paired data](results/lattice_final_all.jsonl) retain the failures.

## Measured result

These examples all use full output, uniform signed 16-bit coefficients and
seed 20261201. Times include each method's complete call and materialization.
Speedups are medians of paired ratios, rather than ratios of the time medians.

| Family | Shape | q cap | e | SYRK ms | Native GEMM ms | Sequential ms | SYRK/seq, 95% block interval |
|---|---:|---:|---:|---:|---:|---:|---:|
| Accumulated innovations | 768×8192 | 2 | -15 | 19.6336 | 20.5861 | 9.0460 | 2.169 [2.164,2.175] |
| Accumulated innovations | 768×8192 | 4 | 0 | 30.4593 | 31.1561 | 15.2613 | 2.045 [1.889,2.229] |
| Corrected windows | 512×4096 | 2 | -15 | 3.2515 | 3.6209 | 2.0573 | 1.596 [1.554,1.612] |
| Accumulated innovations | 512×4096 | 2 | -15 | 2.9476 | 3.3151 | 1.9707 | 1.541 [1.438,1.588] |
| Accumulated innovations | 384×2048 | 4 | -15 | 1.2021 | 1.3621 | 1.2632 | 0.952 [0.943,0.982] |
| Accumulated innovations | 192×768 | 2 | -15 | 0.0820 | 0.0858 | 0.1353 | 0.612 [0.606,0.614] |

For the first row the paired native-GEMM speedup is 2.277 [2.272, 2.282]. The earlier
general-real certified member took 13.2805 ms; the grouped and merged lattice
members took 11.8929 ms and 13.6443 ms. All six complete time series are retained.
The sequential member was selected before evaluation, rather than choosing a
different candidate after seeing each case.

Across all 60 reserved 768×8192 comparisons, both scales, q caps 1/2/4, five seeds
and lower/full output, paired SYRK speedups range 1.610–2.840. For full output,
the median run speedups at q=1/2/4 are 2.477/2.166/1.811, respectively. At
512×4096, the 50 full-output distribution comparisons all beat both vendors;
paired SYRK speedups range 1.455–2.002. At 384×2048, q=1/2 win in every tested run,
while q=4 loses to SYRK. At192×768 every tested case loses. These statements
describe the fixed evaluation distributions and parameters, not all possible
support placements or an unmeasured continuous crossover curve.

The reduction in dot-product work explains the opportunity; exact arithmetic
eliminates certificate/norm/repair work in this member, while batched validation
and fused input traversal keep the dense scan affordable. The smaller cases
show the remaining scan, gather and dependency costs. There is no timing
ablation attributing a separate speed percentage to the schedule alone.
Numerical exactness of that schedule and end-to-end speed are separately proved
and measured properties.

## Algorithm

For `i>=1` and `k<n-1`, define

```text
U[i,k] = A[i-1,k+1],
D[i,k] = A[i,k] - U[i,k],
S_i = supp(D[i,:]),  q_i = |S_i|,  Q = sum_(i>=1) q_i.
```

The first row has no residual list. Compute every first-column entry `G[i,0]`
as a direct dot product. The C implementation uses eight independent FMA
chains and a reduction tree combining disjoint coordinate subsets. Validate
the grid, extract the residual lists, and compute these anchors in one pass.
For increasing `i=1..m-1`, then `j=1..i`, execute exactly:

```text
v = G[i-1,j-1]
v = fma(-A[i-1,0], A[j-1,0], v)
v = fma(A[i,n-1], A[j,n-1], v)
for k in S_i: v = fma(D[i,k], U[j,k], v)
for k in S_j: v = fma(A[i,k], D[j,k], v)
G[i,j] = v
```

Complete the whole first correction loop before starting the second. Overlapping
supports require no merge or special branch. Reflect the lower triangle if a
full matrix is required. For `m=0` there is no output; for `n=0` write lower
zeros. For `n=1` both support lists are empty and the endpoint updates suffice.

## Complete real correctness proof

The predecessor dot is

```text
G[i-1,j-1] = A[i-1,0]A[j-1,0]
              + sum_(k<n-1) U[i,k]U[j,k].
```

Removing the head and inserting the tail therefore leaves

```text
A[i,n-1]A[j,n-1] + sum_(k<n-1) U[i,k]U[j,k].
```

For any prefix `P` of the first support loop, let `L_P[k]=A[i,k]` when
`k in P`, and `L_P[k]=U[i,k]` otherwise. The accumulator invariant is

```text
v(P) = A[i,n-1]A[j,n-1] + sum_(k<n-1) L_P[k]U[j,k].       (1)
```

Initially the prefix is empty. Adding `(A[i,k]-U[i,k])U[j,k]` changes just
one factor in (1) to `A[i,k]`. This proves the invariant after every first-loop
step. At its end `L_P[k]=A[i,k]` everywhere, because coordinates outside the
support already agree.

For a prefix `T` of the second loop, define `R_T[k]=A[j,k]` when `k in T`,
and `R_T[k]=U[j,k]` otherwise. Its invariant is

```text
v(T) = A[i,n-1]A[j,n-1] + sum_(k<n-1) A[i,k]R_T[k].       (2)
```

Adding `A[i,k](A[j,k]-U[j,k])` changes the next right factor, proving (2).
At the end, (2) is exactly the dot product of rows i and j. The predecessor
`(i-1,j-1)` belongs to an earlier computed row. Exact first-column anchors and
induction on i prove every lower entry. The empty and one-column cases above
close all dimension boundaries.

No factors were commuted. Equivalently,
`Ui Uj + Di Uj + Ai Dj = (Ui+Di)Uj + Ai Dj = Ai(Uj+Dj)`.
With the appropriate adjoints, this is also an ordered block identity. The
delivered algorithm is a scalar-entry recurrence, however; no unimplemented
recursive block speedup is inferred from that observation.

## Complete binary64 exactness proof

Assume `A=2^e X`, `|X|<=B`, `nB^2<=2^53`, and the exponent range above.
All anchor-chain states and reduction nodes combine disjoint subsets of at
most n input products. Their coefficient in units `2^(2e)` is an integer of
magnitude at most `nB^2`. The endpoint deletion leaves n-1 such products;
the insertion restores n. Every prefix in (1) or (2) contains n products
of actual old or new input factors, each bounded by B in coefficient units.
Thus the same bound holds after **every** arithmetic operation, including
overlapping supports and diagonal entries.

Every integer of magnitude at most `2^53` is representable in binary64.
Scaling these coefficients by `2^(2e)` preserves representability: nonzero
accumulator results are at least `2^-800`, and all magnitudes are at most
`2^453`. They are normal and finite. Correctly rounded FMA and addition
therefore return the exact result. This argument inducts through every
operation and then through every recurrence row. There is no accumulated
rounding drift, cancellation exception, or dependence on chain length.

Residual coefficients have magnitude at most `2B<2^28`, so their subtractions
are exact, normal when nonzero, and finite. Whenever a residual exists,
`n>=2`, and even its product coefficient is bounded by `2B^2<=2^53`.
At `n=1` no residual product occurs. The proof does not require that FMA rescue
an unrepresentable multiplication in this admitted class.

The C function obtains `B_n` using an integer quotient and corrected integer
square comparisons; a rounded square root is only an initial estimate. It
validates both `x=A*2^-e == trunc(x)` and `x*2^e == A`, along with the magnitude
cap. The reconstruction test rejects tiny invalid values that scale to zero.
It rejects NaNs, infinities, off-grid values, excess defects, unsupported sizes,
and unsafe allocation arithmetic. It never changes input values to make them
eligible. Strict floating-point compilation is required.

Consequently every successful exact-member call satisfies

```text
computed G[i,j] = sum_k A[i,k]A[j,k] exactly over the reals.
```

Its numerical error is zero, including nearly dependent rows and exact
cancellation zeros. This exceeds the unchanged external requirement
`|error_ij| <= 2e-12 ||A_i||_2 ||A_j||_2`. Unsupported inputs are rejected;
the separate general-real certified member can be called under its own domain.
There is no automatic fallback dispatcher in the measured implementation.

### Why ordering matters

Let `n=2`, `B=2^26`, `e=0`, rows `[0,-B]` and `[B-1,B]`.
The final diagonal `9007199120523265` is exactly representable and
`nB^2=2^53`. The earlier grouped order `Di Aj + Ui Dj` nevertheless creates
the odd intermediate `18014398308155393`, which rounds down by one. Its final
answer is `9007199120523264`. The sequential order passes through `67108864`
and returns the exact diagonal. Thus a final-dot magnitude bound alone does
not justify an arbitrary update order. The actual FMA witness is independently
checked in the supplied verifier.

## Work, storage, and traffic

Let `P=m(m+1)/2`, `D=P-m`, and `c=min(max_q,n-1)`. The exact number of Gram
FMAs for nonempty dimensions is

```text
F = mn + 2D + sum_(1<=j<=i<m) (q_i+q_j)
  = mn + 2(P-m) + mQ.
```

Residual row t appears t times as the first correction row and m-t times as
the second, giving exactly m uses. Validation separately performs `2mn` grid
multiplications, comparisons and integrality tests. The anchor tree adds `7m`
source-level additions; extraction performs Q subtractions. With small q the
work is `O(mn+m^2+mQ)`, compared with `nP` classical SYRK product terms.
This is a specialization of a known displacement-product mechanism, not a
new asymptotic result for unrestricted dense matrices.

On this 64-bit ABI requested heap workspace is

```text
16(m-1)c + 8(m+1) bytes, if (m-1)c > 0;
0 bytes otherwise.
```

There are no norm/error arrays, packed input, second Gram matrix, or repair
dots. Dense A/G allocations remain `8mn` and `8m^2` bytes. At768×8192,q2,
these are 48 MiB and 4.5 MiB, with 30,696 bytes of scratch; if every later row has
two defects, the Gram body uses 8,058,624 FMAs versus 2,419,064,832 classical SYRK
product terms. Admission still performs 12,582,912 multiplications.

The fused scan loads the current row, first-row anchor factors, and shifted
predecessor factors. Changed vector blocks can require scalar rereads. A
reuse-aware preparation operand model is
`8[2mn+(m-1)(n-1)]` bytes, plus at most 128 bytes per changed eight-value block.
The recurrence logically reads D predecessor G values and writes P outputs;
their combined bytes are `8m^2`. Endpoint factors, dense correction gathers,
and CSR records add logical operand demands `32D`, `8mQ`, and `16mQ` bytes.
Offsets and control operations add further costs. Caching and register reuse
can eliminate many physical loads. These are source-access counts, **not
measured DRAM traffic**. Full reflection adds `16D` logical bytes. Vendor-internal
packing/workspace and process peak memory were not instrumented; no numerical
claim comparing total vendor and candidate peak memory is made.

[The detailed traffic analysis](notes/sequential_traffic.md) includes padded
strides, all empty boundaries, changed-block rereads and comparison with the
TwoDiff member. For example, its sparse scratch at 768 rows, q=2 is 55,256 bytes,
versus 30,696 here. On lattice inputs its low residual limbs are zero, so its
core correction FMA count can equal the new member's count. Removing norm and
certificate work and batching admission predicates are also material. No claim
that a hundreds-fold arithmetic reduction implies a comparable runtime gain
is made, and no timing ablation isolates every individual optimization.

## Experimental protocol

The kernel was frozen before reserved evaluation; its SHA256 is
`7550990a7f6dd82fbd0e58580c9a6639eef7ffc07d841b56ae9cf189b394bde5`.
The development cases, including losing initial kernels, are retained as
`results/lattice_dev*.json`. No kernel parameter depends on held-out results.
The complete pre-timing manifest is `results/lattice_frozen_manifest.json`.

All six methods receive the same dense row-major double input, double output,
precision, one-thread setting, and output contract: Accelerate DSYRK, native
DGEMM, the earlier certified floating member, grouped lattice, merged-support
substitution, and the preselected sequential member13. No per-case selection
among the three lattice candidates is allowed. The old floating member uses
its Hankel-specialized path at q=0 and sparse path otherwise. All final jobs
here have q>0. Same-operation SYRK is the main comparator; GEMM is also reported
so a loss to a faster full-output vendor alternative cannot be hidden.

The environment is Apple M4, macOS 27.2 build 26B5091g, 16 GiB, Apple clang 21,
Accelerate FP64. `BLASSetThreading` sets single-thread mode, verified by
`BLASGetThreading`, with `VECLIB_MAXIMUM_THREADS=1`. Threads are not pinned.
The per-call clock includes input validation, residual extraction, anchors,
temporary allocations/frees, all recurrence work, and full reflection where
requested. A shared 32×32 blocked mirror is used for triangular methods; native
GEMM already computes full output and receives no redundant mirror.

Each case has at least five warmups per method, batch calibration to at least
10ms, and 30 paired samples in five balanced randomized cyclic blocks of six.
The operands are warm and repeated. Other research compute is paused during
timing. Source and binary hashes are checked between cases. Power, OS, linked
libraries, compiler and end-of-run thermal reports are retained in the environment
receipt; ordinary macOS scheduling variability remains. This Round2 run used
battery power (39% at start, 31% at end), unlike the earlier AC-powered Round1.
No thermal/performance warning was reported by the final system query. Absolute
times from the two rounds should not be treated as a controlled power comparison.

There are 180 reserved uniform-innovation comparisons: shapes 192×768, 384×2048,
768×8192; q=1/2/4; e=0/-15; five seeds 20261201..20261205; lower/full output.
Another 50 full-output comparisons use 512×4096, q=2, e=-15, both independent corrected
windows and accumulated innovations, with uniform signed 16-bit, nearly constant
integer, saturated alternating, saturated impulse, and ternary coefficients.
These are 230 comparison runs, not 230 independent input draws: lower/full reuse
inputs, scales reuse integer coefficients, and deterministic window patterns
repeat across seeds. Accumulated innovations are generated by shifting each
previous row and inserting/replacing coordinates; they need not be a Hankel
matrix with only a few edits per complete row.

Before timing, NaN poisoning and a full finiteness scan check output completeness;
an independent int64 oracle checks first 128 diagonals plus 256 seeded pairs, and
reflected triangular outputs are checked for symmetry. For full-output jobs, native GEMM is checked
against the oracle samples in both triangles. The exact-member comparisons demand zero
error. Vendor/floating comparisons use the original 2e-12 norm-scaled threshold.
With signed 16-bit coefficients, all oracle sums fit in signed64bits with
absolute value at most2^46, so conversion/scaling is exact. This sampled
preflight supplements the universal proof and independent exact verification;
it is not a full-entry oracle for every timed matrix.

Reported speedups are medians of paired per-trial ratios. Intervals use 5000
resamples of the five contiguous six-trial balanced blocks. They describe
within-run variability, with only five blocks; they are not guarantees across
machines, repeated days, or all 230 comparisons simultaneously. Both vendor
comparisons are reported. A best-vendor interval additionally conditions on
the vendor selected by that run's median, not a per-sample oracle.

## Established results and contribution

| Source | Established content | Contribution in this investigation |
|---|---|---|
| [Pan–Murphy–Rosholt,2008,§2.4Theorem 2.3](https://www.researchgate.net/publication/227235518_Unified_Nearly_Optimal_Algorithms_for_Structured_Integer_Matrices) | Arbitrary-rank product displacement; sparse Gram algebra and cost specialize from it | A particular exact binary64 execution schedule and its complete prefix guarantee |
| [Dumas–Giorgi–Pernet,FFLAS/FFPACK,§§2.1.2,3.2.3](https://arxiv.org/pdf/cs/0601133) | Bounded integer arithmetic in floating point; algorithm-specific intermediate bounds | `nB^2<=2^53` for this sparse schedule independently of q/m, with checked admission |
| [STOMP,2016,Eq4](https://mcyeh.github.io/paper/2016_icdm_exploiting_novel_algorithm.pdf) | q=0 sliding-dot recurrence | Arbitrary row-specific sparse factor changes under the stated exactness contract |
| [SCAMP,2019,§§3.2.4,3.3](https://www.cs.ucr.edu/~eamonn/SCAMP-camera-ready-final1.pdf) | Tile seeds, stability-oriented scheduling and centered updates | This no-restart exact lattice schedule; not numerical scheduling as a general idea |
| [flip-cpd](https://github.com/khoruzhii/flip-cpd), [Khoruzhii–Gelss–Pokutta v2](https://arxiv.org/html/2511.10786v2) | Published 34-product ordered identity and mixed upper-product recursion | Independent certificates/reference; no novel rank or measured speed claim for this route |

The displacement theorem has **no low-rank prerequisite**; a full-rank residual
does not make the algebra new. [Its explicit subsumption derivation](notes/priority_displacement.md)
corrects the earlier open comparison. The schedule is an elementary constructive
corollary of known bilinear and representability principles, whose useful extra
content is the simultaneous sparse work bound, ordinary dot capacity, no-merge
execution and fully checked software. [The independent bounded contribution
review](notes/contribution_sequential.md) reconstructs that argument and the
strict grouped-order counterexample. It does not infer historical priority
from an unsuccessful search.

The investigation began with distinct algebraic, flip, shared-expression,
communication, numerical and storage routes. [The registry](REGISTRY.md)
groups their proofs, failures and reopening conditions. The flip experiments
used equal two-second single-thread CPU budgets in F2/F3; their modular
rank-17 outputs were independently checked and never treated as real algorithms.
The published characteristic-zero scheme has a separately verified ordered
certificate and executable recursion. No search-method superiority is inferred
from those short equal-budget walks.

The earlier five-product scalar candidate was independently checked. Its saved
timing ratio is 1.2088068, confirming the reported 21% slowdown in that historical
benchmark, not a new Mac measurement. It returns zero for a diagonal whose exact
value is 2 in the large-scale cancellation example, and noncommuting block
counterexamples reject an automatic block recursion. Those failures concern
that construction and are not an impossibility result for this research.

## Verification and reproduction

Source: [src/lattice_sequential.c](src/lattice_sequential.c).
Public declaration: [gram.h](gram.h).
Full arbitrary-prefix proof, pointer/allocation audit and operation details:
[notes/lattice_sequential.md](notes/lattice_sequential.md).

The independently written verifier imports no search procedure. Native and
forced-scalar builds each pass 2,044 exact lower entries, all short vector tails,
coefficient caps, exponent endpoints, q=64, padded strides, untouched upper
entries/padding, invalid grid/NaN/subnormal cases, and the strict witness.
It additionally checks 6,561 finite invariant states, an ordered polynomial
certificate and 371 CSR-group prefixes. These finite checks supplement the proof.
The whole original suite plus all three new members passed after integration;
the transcript is `results/round2_verification.log`.

A second independently written integer-first generator/oracle checks every
entry of four representative large matrices, with contiguous and padded
layouts, after all timing is complete. The new member passes 1,263,744 lower-entry
comparisons; SYRK passes the same count; native GEMM passes 2,523,136 full entries.
All 5,050,624 comparisons are exactly equal, with zero nonfinite outputs, missing
writes, input changes or padding violations. Thus the vendor baselines also
achieve zero error on these cases; the comparison does not give the candidate
an easier accuracy target. [The full-entry receipt](verification/algebra_lattice_full_entry_results.json)
records all calls, lossless error checking, compiler commands and stable hashes.
The direct oracle is integer multiplication/summation, not another displacement
recurrence. These four matrices do not replace the universal proof or claim
full-entry oracle coverage for every timed input.

```sh
cd /Users/leonard/Downloads/onderzoek_drie_problemen/gram_research
./build.sh
./verify.sh
clang -O3 -std=c11 -fno-fast-math -ffp-contract=off \
  example_exact.c src/lattice_sequential.c -lm -o example_exact
./example_exact
python3 verification/algebra_lattice_full_entry_verify.py --run

# Preserve the supplied timing receipts before reproducing the full grid.
# Run timing separately from verification and other compute.
mkdir -p results/recorded_round2
mv results/lattice_final_all.jsonl results/lattice_environment_all.json results/recorded_round2/
VECLIB_MAXIMUM_THREADS=1 python3 run_lattice_benchmarks.py --grid all
python3 summarize_lattice.py
```

No Python packages are needed. The numerical verifier builds scalar/native
libraries with warning-as-error strict flags. Accelerate is required for the
recorded benchmark; scalar implementation portability does not imply a speedup
on other platforms. The runner freezes environment/hashes and 230 jobs before
execution. The summarizer rejects incomplete manifests, missing/duplicate jobs,
incorrect method/sample sets and overlapping run files.

The API writes the lower triangle only. A/G must be disjoint valid allocated
objects, with strides in doubles. Status 0 certifies exact output; negative
statuses report invalid objects/dimensions, invalid values or excess defects,
allocation failure, or unsupported lattice/parameter class. Discard partially
written output on error. The separate general-real certificate API retains its
own 2e-12 contract and numerical domain; rejection by either member never establishes
eligibility for the other.

## Limits

The speed claim is limited to completed measurements in the specified structured,
warm, single-thread FP64 regime. Small shapes and larger defect costs may lose.
There is no unrestricted dense asymptotic improvement, new rank bound, universal
runtime win, or tested claim for multithreading, GPUs, other hardware, cold
operands, energy consumption, or application-level throughput. Actual physical
memory traffic was not instrumented. Unsupported off-grid or out-of-range inputs
are not projected. Absolute historical first occurrence remains unclaimed; the
concrete theorem, implementation and finite prior comparison define the result.
