# Reconstructed Hankel plus sparse correction comparator

This is an independently implemented comparator derived from established sliding
inner products (STOMP) and the elementary Hankel-plus-sparse decomposition described below. It is **not** a published
state-of-the-art implementation, a historical-priority finding, or a claim that
this particular Hankel reconstruction is optimal. No external implementation was
copied. Original kernels and all original 230 measurements remain unchanged.

## Input and fairness

`benchmarks/specialized/src/hankel_sparse.c` accepts the same dense, row-major,
strided matrix, lattice exponent, local-shift defect cap and lower-output contract
as `src/lattice_sequential.c`. It validates actual dense input against
`A[i,k]=2^e X[i,k]`, integer `X`, `|X|<=floor_sqrt(floor(2^53/n))`,
`-400<=e<=200`, `n<=65536`, and at most `max_q<=64` actual interior shift
defects per later row. It preserves the original status classes and empty-size
semantics. It does not quantize. As for the original C API, callers supply
disjoint valid array objects and strict IEEE binary64 operations.

Every call reconstructs a generator of length `m+n-1`:

```
h[k] = A[0,k]                     (0 <= k < n)
h[n+i-1] = A[i,n-1]              (1 <= i < m)
H[i,k] = h[i+k]
E[i,k] = A[i,k] - H[i,k].
```

Thus every H factor is an actual validated A value with the same coefficient
bound B. The code counts actual mismatches, allocates exactly that many sparse
records, and fills E from A. No generating signal, edit list, cached analysis,
or latent metadata is supplied. Validation, H reconstruction, two dense
mismatch traversals, local-defect admission, temporary allocations/frees, core
work and full reflection all belong to the timed call. Separate untimed E
metadata reporting does not feed either algorithm.

The simple reconstruction intentionally matches the explicit protocol request.
It may choose a different base from the latent benchmark generator: a first-row
edit propagates along its anti-diagonal. It does not use an anti-diagonal mode,
minimum-edit base, or a supplied signal. Therefore measured conclusions must
refer to this concrete comparator. The anti-diagonal mode alternative in the
original note remains unimplemented and unmeasured.

## Independent correctness and exactness argument

First compute `B=HH^T`. Direct first-column seeds use eight disjoint FMA chains
and a disjoint reduction tree. For `i>=1, 1<=j<=i`, the established sliding rule is

```
B[i,j] = B[i-1,j-1] - h[i-1]h[j-1] + h[i+n-1]h[j+n-1].
```

Removing the head leaves n-1 products and adding the tail restores n. The entire
H-only lower triangle is finished before any B entry is corrected, so an
in-place correction cannot corrupt a future predecessor. For each lower pair,
start `v=B[i,j]`, then run **in this order**:

```
for k in supp(E_i): v = fma(E[i,k], H[j,k], v)
for k in supp(E_j): v = fma(A[i,k], E[j,k], v).
```

For a prefix P of the first loop, v is `sum_k L_P[k] H[j,k]`, where L uses A
at replaced positions and H elsewhere. At its end v is `A_i dot H_j`. For a
prefix T of the second loop, v is `sum_k A[i,k] R_T[k]`, with R using A at
replaced positions and H elsewhere. Its end is exactly `A_i dot A_j`.
This proves the real identity for **all E**, including dense E and overlapping
supports. No sparsity assumption enters correctness.

Every intermediate seed, reduction, deletion, insertion and correction prefix
is a sum of at most n products of actual A/H factors. Its integer coefficient
in units `2^(2e)` has magnitude at most `nB^2<=2^53`. Every such value is exactly
representable. Nonzero states lie between `2^-800` and `2^453`, hence are normal
and finite. Each exact residual coefficient has magnitude at most `2B<2^28`,
so residual subtraction is exact too. At n=1 the reconstruction equals A and
E is empty; otherwise residual product coefficients are bounded by `2B^2<=2^53`.
Correctly rounded FMA/add/subtract consequently return the exact real values
at every operation. This is an execution-order proof, not an inference from a
small final answer. Cancellation and nearly constant rows do not weaken it.

The dynamic bound uses an integer quotient and integer-square corrections after
a sqrt estimate. Admission checks both scaled integrality and reconstruction;
the latter rejects nonzero invalid values that scale to zero. Tests independently
check all entries using scalar integer products, compare the original sequential
kernel, check untouched padding/upper output, and cover arbitrary dense E,
exponent endpoints, q=0, n=0/1, the grouped-order boundary witness, NaN/infinity,
off-grid/subnormal values, excessive q, coefficient overflow and invalid strides.

## Cost and limitations

Write `R=nnz(E)`, `P=m(m+1)/2`. For nonempty dimensions the core uses

```
mn + 2(P-m) + (m+1)R FMAs
```

and `7m` reduction additions. Each E row appears `(i+1)` times in the first
correction loop and `(m-i)` times in the second, totaling m+1 uses. Preparation
is O(mn), so total work is O(mn+m^2+mR). Requested heap workspace is
`16R + 8(m+n-1) + 8(m+1)` bytes on the current 64-bit ABI, besides A and G.
There is no second Gram matrix. Full output uses exactly the same blocked-32
mirror as the original benchmark; native GEMM receives no redundant mirror.

Accumulated innovations can make R large even when every local transition has
q=1 or q=2. This is a real parameter difference: the present measurements record many reconstructed corrections despite local q=2. For such
inputs this comparator can require O(m^2 n) correction work and substantial
scratch. That obstacle must be disclosed alongside results; it does not imply
that all specialized structured algorithms have this cost. Even window inputs
may have more reconstructed edits than their latent one-edit representation.

## Separate, predeclared experiment

The driver compares exactly Accelerate SYRK (0), native GEMM (4), the original
preselected sequential member (13), and this H+E comparator (14). The input
generation is copied from the original project-owned Round2 helper, preserving
its RNG, draw order and coefficient distributions. The small grid has 20 jobs:
512x4096, q cap 2, e=-15, full output, corrected windows/accumulated innovations,
uniform signed16/nearly constant integer coefficients, seeds 20261201..20261205.
No results from it may be merged into the original 230 without a separate label.

The four-method schedule uses 32 paired samples in eight randomized balanced
cyclic blocks of four, at least five warmups and >=10ms batch calibration per
method. Four methods require a different block size from the original six.
Use block bootstrap with these eight four-sample blocks for uncertainty.
All methods receive identical dense A and double output, single-thread BLAS,
NaN poisoning, full finite-output checks and the original independent int64
sample oracle. Exact equality is required for methods 13/14; vendors retain
the common 2e-12 norm-scaled error threshold. Original small full-entry tests
supply additional exact checks outside timing. Implementation checks were untimed. After source freeze, the separate 20-case
measurement campaign completed; its files are under `data/specialized_2026-10-19/`.

From `benchmarks/specialized`:

```
make all check
python3 run.py                         # preview only, never times
python3 run.py --freeze frozen.json     # exclusive creation; no timing
# Only after the operator approves the frozen grid and quiet serial window:
python3 run.py --run frozen.json --output-dir results
python3 summarize.py results/specialized_raw.jsonl --manifest frozen.json --output results/summary.csv
```

The runner fingerprints its sources, unchanged candidate, build flags (via the
Makefile), itself and the binary; checks fingerprints between cases; refuses
receipt overwrites; and records environment, power, libraries and final thermal
status. It never freezes or times by default. The original benchmark directory
and artifacts are untouched. Compilation/correctness testing can happen while
other work runs; timing must wait until competing compute is paused.

Primary antecedent: Zhu et al., *Matrix Profile II* (ICDM 2016), Eq. (4),
https://mcyeh.github.io/paper/2016_icdm_exploiting_novel_algorithm.pdf . The
source-specific comparison and limitations are stated above; the sequential
mixed-factor exactness argument above
is independently reconstructed for this comparator rather than delegated to
STOMP's floating-point behavior.

### Harness regression check before freeze

A pre-freeze review caught an argument-index regression while reducing six method
arrays to four: the seed was inadvertently read from q's argument. This was fixed
before any freeze or timing. Every argument now has checked numeric parsing.
`make check` runs a no-clock input diagnostic that exits before kernel calls,
warmup or calibration. An independent Python implementation reconstructs the RNG,
coefficient draws, sparse-coordinate draw order, complete binary64 input
fingerprint and balanced method order on a small grid. All 20 original seeds/
family/distribution combinations print the requested seed, match that reference,
and produce five distinct fingerprints per family/distribution. The saved
`tests/harness_diagnostic.json` contains these untimed diagnostic records.

Freeze also captures actual full-grid input fingerprints and method orders. Each
measurement must match its frozen job, seed, input fingerprint, full-output flag,
method order and hashes. The summarizer requires the original frozen manifest and
completed environment receipt, verifies historical source hashes and recorded
binary-hash correspondence, and rejects
incomplete grids, duplicate job indices, changed seeds, inconsistent receipt
hashes, invalid samples or unbalanced schedules before emitting statistics.

### Portable analysis of the completed 20-case campaign

The 20 additional measurements completed under
`data/specialized_2026-10-19/`; they remain separate from the original 230.
Reproducing their CSV is standard-library Python analysis and needs no macOS
executable, Accelerate library, compiler, or new timing:

```
python3 benchmarks/specialized/summarize.py \
  data/specialized_2026-10-19/specialized_raw.jsonl \
  --manifest data/specialized_2026-10-19/frozen.json \
  --output /tmp/specialized-reproduced.csv \
  --integrity-report /tmp/specialized-integrity.json
```

This command checks every historical source SHA256 against available file bytes.
The pre-portability frozen summarizer is preserved under
`benchmarks/specialized/frozen_sources/<its-sha256>/summarize.py`, so updating the
analysis entry point does not discard the original frozen timing/analysis source.
Every other frozen source must likewise match available current or archived
bytes; a missing or altered historical source is an error.

The binary SHA256 must agree between the frozen manifest, completed environment
receipt and all 20 individual measurement receipts. A local executable is
optional: if present, its hash is reported independently and may differ after
rebuilding or changing platform. That does not change the archived data. The
recorded receipts attest to the historical binary and execution; hash agreement
alone does not independently prove that execution occurred, or reproduce the
historical compiler/Accelerate environment. The optional integrity JSON makes
that trust boundary explicit. The runner continues to demand matching local
source **and binary** hashes for any new timing campaign.
