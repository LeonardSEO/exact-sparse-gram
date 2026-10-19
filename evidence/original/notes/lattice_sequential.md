# Exact sequential factor updates without a support merge

Research date: 19 October 2026. This is a separate implementation of the
sparse row-shift Gram construction. Updating the two factors sequentially
gives the same sufficient bound `nB^2<=2^53` as coordinate substitution,
using the two original CSR support loops and no merge branch. There is no
historical-priority or performance claim in this note.

## Theorem and exact operation order

Let `m,n` be nonnegative integers, `e` an integer in `[-400,200]`, and `B` a
nonnegative integer. For `n>=1`, assume

    A[i,k] = 2^e X[i,k],       X[i,k] integral, |X[i,k]| <= B,
    n B^2 <= 2^53.

For `i>=1` and `0<=k<n-1`, let

    U[i,k] = A[i-1,k+1],
    D[i,k] = A[i,k] - U[i,k],
    S_i = supp(D[i]).

Compute the first-column dot products with FMA chains followed by a
reduction tree combining disjoint product subsets. For every `1<=j<=i<m`,
in row order, compute

    v = G[i-1,j-1]
    v = fma(-A[i-1,0], A[j-1,0], v)
    v = fma(A[i,n-1], A[j,n-1], v)
    for k in S_i: v = fma(D[i,k], U[j,k], v)
    for k in S_j: v = fma(A[i,k], D[j,k], v)
    G[i,j] = v.

Under correctly rounded IEEE binary64 operations, without reassociation
that changes this order or weakens validation, all residual subtractions,
anchor operations, and recurrence FMAs are exact over the reals, with finite
results and normal nonzero results. The output is exactly the lower triangle
of `AA^T`, including true zero inner products. No mathematical support-size
or row-count bound is needed for the numerical theorem. These parameters
only control implementation work and storage.

The coordinate identity is

    Ui Uj + Di Uj + Ai Dj
      = (Ui+Di) Uj + Ai Dj
      = Ai (Uj+Dj)
      = Ai Aj.

This differs from the old grouped order `Di Aj + Ui Dj`. Both have the same
final polynomial. Their intermediate coordinate factors differ, which is
decisive for binary64 exactness.

## Complete CSR-group prefix invariants

Assume `G[i-1,j-1]` is exact. Removing the head product leaves
`sum_(k<n-1) U[i,k]U[j,k]`, a sum of `n-1` bounded products. Adding the new
tail gives `n` bounded products. In coefficient units each product has
magnitude at most `B^2`.

For a prefix `P` of the first CSR row `S_i`, define

    L_P[k] = A[i,k] if k is in P, else U[i,k].

The first-group invariant is exactly

    v(P) = A[i,n-1] A[j,n-1]
           + sum_(k<n-1) L_P[k] U[j,k].                 (1)

Initially `P` is empty. At the next support coordinate `k`, adding
`D[i,k]U[j,k]` changes its factor from `U[i,k]` to `A[i,k]`, proving (1)
after that FMA. Every factor in (1) is an actual old or new input element
with coefficient magnitude at most `B`; hence the total coefficient is
bounded by `nB^2`. Once this group finishes, `L_P[k]=A[i,k]` at every
coordinate: outside `S_i`, the factors were equal already.

For a prefix `Q` of the second CSR row `S_j`, define

    R_Q[k] = A[j,k] if k is in Q, else U[j,k].

The second-group invariant is exactly

    v(Q) = A[i,n-1] A[j,n-1]
           + sum_(k<n-1) A[i,k] R_Q[k].                 (2)

Adding `A[i,k]D[j,k]` changes the next right factor to `A[j,k]`, proving
(2) after every second-group FMA. The same `nB^2` bound holds. At the end,
all right factors equal `A[j,k]`, so the entry equals the true inner product.

Overlapping supports require no branch or special formula: the first group
changes the left factor while the right one is old; the second changes the
right one while the left one is new. Each coordinate occurs once in each
applicable CSR row. The order within a group can be arbitrary; the
implementation uses the increasing scan order. The completion of the first
group before starting the second is load-bearing.

Anchor chain states and reduction nodes are sums of disjoint subsets of at
most `n` original products. Every anchor and recurrence prefix thus has a
coefficient of magnitude at most `nB^2<=2^53`, an exactly representable
integer. After scaling by `2^(2e)`, every nonzero result is at least `2^-800`
and every accumulator has magnitude at most `2^453`, so all are normal and
finite. Correct rounding returns their exact values. Induction on operations
proves exactness within an entry, and row induction supplies the exact
predecessor for every next entry. There is no drift with recurrence length.

The bound implies `B<2^27`, so every residual coefficient has magnitude
`<=2B<2^28` and is exactly representable. Its nonzero scaled magnitude is
at least `2^-400` and its magnitude below `2^228`, proving exact residual
subtraction. If a residual coordinate exists then `n>=2`, and its product
coefficient satisfies `2B^2<=nB^2<=2^53`. At `n=1` supports are empty.
Thus residual products also remain representable in this admitted class.

## Strict witness against the old grouped order

Set `n=2`, `e=0`, `B=2^26`, and use

    A[0] = [0, -B],
    A[1] = [B-1, B].

There is one residual `D[1,0]=2B-1`, and `nB^2=2^53`. Both schedules reach
`v=2B^2` after the boundary operations. The old first correction computes

    2B^2 + (2B-1)(B-1) = 4B^2-3B+1 = 18014398308155393,

which round-to-nearest ties-to-even changes to `18014398308155392`. The
second correction therefore gives `9007199120523264`, one below the exactly
representable true value

    (B-1)^2+B^2 = 2B^2-2B+1 = 9007199120523265.

The sequential order first uses the old right factor `-B`:

    2B^2 + (2B-1)(-B) = B = 67108864.

Its second correction adds `(B-1)(2B-1)`, giving the exact true diagonal.
Both FMAs are exact. The verifier checks this witness with actual FMA calls
and through the compiled kernel. The new `nB^2` condition is strictly more
permissive than the old grouped sufficient bound `(n+2q)B^2<=2^53`, even
at `q=1`; it matches the standard bounded-subset sufficient condition for
classical dot accumulation. Neither condition is necessary for every
individual input, since larger integers with trailing binary zeros may also
be exactly representable.

## Batched validation, extraction and pointer bounds

`src/lattice_sequential.c` exports

```c
int lattice_sequential_gram(const double *A, size_t m, size_t n, size_t lda,
                           double *G, size_t ldg, int scale_exp, size_t max_q);
```

The supported nonempty class has `n<=65536`, `max_q<=64`, `e` in
`[-400,200]`, and the n-dependent coefficient cap

    B_n = floor_sqrt(floor(2^53/n)).

As in `lattice_substitute.c`, the quotient and final square comparisons use
`uint64_t`; `sqrt` provides only an initial estimate. Every candidate square
is below `2^54`. The adjusted integer inequality, not a rounded root alone,
establishes `nB_n^2<=2^53`. The scalar/vector scan checks every actual input
for the cap and exact lattice roundtrip, without projection or trusted
metadata.

On Arm64 the main loop condition is `k+8<n`, leaving a scalar tail of one
through eight elements. Four two-lane loads cover `a[k]` through `a[k+7]`.
Their `abs(a)<=bound` comparison masks are combined before **one** range
branch. Because bound is finite, the ordered comparisons also reject NaN
and infinities. Four pairs of integrality and reconstruction masks are
combined before **one** grid branch. The grid test checks

    x = a*2^-e,  x == trunc(x),  x*2^e == a.

Reconstruction rejects tiny invalid inputs that scale to zero. Only after
a block mask fails does a scalar scan of those eight elements determine an
error class. Its defensive final return rejects if scalar/vector predicates
were ever inconsistent; it cannot accept a partially scanned row.

For a nonfirst row, four shifted loads cover `previous[k+1]` through
`previous[k+8]`. The strict loop predicate gives `k+8<=n-1`, so all lanes
are valid. Eight shifted equality comparisons combine into one sparse-change
branch. If any differ, scalar extraction checks all eight positions in
increasing order, storing only actual nonzero residuals and enforcing the
row cap before each store. An all-equal block has no scalar extraction.
The tail checks `k+1<n` before reading a shifted predecessor. It never
reads `previous[n]`. Every coordinate is visited exactly once, and the last
coordinate has no shifted residual.

The eight independent exact anchor chains are preserved. Four NEON vectors
store to eight scalar chains before the tail; a tail of eight uses indices
`0..7` only. The forced-scalar path uses the same block/tail division. The
final tree combines disjoint feature subsets. Previous rows and row zero
are fully validated before use by a later row.

CSR storage has `(m-1)c` entries and `m+1` offsets when its capacity is
positive, with `c=min(max_q,n-1)`. `rows[i]..rows[i+1]` is exactly row i's
support, including an empty row-zero interval. Validation failure bypasses
the recurrence, so uninitialized later offsets are never used. Every append
checks its per-row count; the total capacity follows from at most `c`
appends in each of `m-1` rows. Residual coordinates obey `k<n-1`, so the
recurrence's `uj[k+1]` and `ai[k]` accesses are valid. A/G extents and
allocation byte arithmetic are checked before use. Valid allocated objects
and nonaliasing are caller preconditions.

Each recurrence entry uses exactly `2+q_i+q_j` FMAs, including overlapping
supports. The total product count remains

    mn + sum_(1<=j<=i<m) (2+q_i+q_j).

This holds for q=0, q=1 and every admitted support pattern. No merge, dense
transform, norm buffer, error buffer or additional per-entry storage is used.
Workspace on a 64-bit ABI is `16(m-1)c+8(m+1)` bytes for positive capacity,
and zero otherwise. All validation, extraction, bound calculation and
allocation occur inside the API call and belong in any later timing.

Returns are `0` success, `-1` invalid pointer/size/stride/extent, `-2`
nonfinite input or too many defects, `-3` allocation failure, and `-4`
outside the admitted lattice/parameter class. Discard partial output on
error. `m=0` is an unconditional no-op; `n=0` checks output extents and
parameter caps and writes lower zeros. At `n=1` the predecessor cancels,
then the sole new product is inserted. Only the lower triangle is written;
strides are measured in doubles.

## Source provenance and untimed verification

The new source derives from the project-local `src/lattice_substitute.c` at
SHA-256 `d1400c003886f0e57102e75bb7e442b7329cbad26d6d796b844b9bf877f1ba95`,
which derived from `src/lattice_exact.c` at SHA-256
`3402d7c6531fdb17178493df4205b54aaa580fcc6f6d99aa3a66cffa3e3da540`.
Both previous files were preserved. No external implementation or additional
license grant is introduced; the project owner's licensing decisions remain
applicable. Changes are the sequential correction factors, combined block
validation/extraction predicates, reserved scalar tail, and distinct API/
forced-scalar macro names.

`python3 verification/sequential_verify.py` builds native and forced-scalar
implementations with `-O3 -std=c11 -fno-fast-math -ffp-contract=off` and
`-Wall -Wextra -Wconversion -Werror`. The receipt
`verification/sequential_results.json` records source/verifier hashes,
commands, cases, rejection classes and exact witness values. Both paths
passed **2,044 lower-entry comparisons each** against direct Python integer
Gram products scaled by exact Fractions. The cases include n-dependent
maximum caps, all lengths `n=1..25`, exponent endpoints, q=64, a 37-row
chain, maximum n, signed extremes, zero rows, padded strides, and unchanged
upper/output padding. Another 86 cases put exactly one defect at every
shifted lane for `n=9,16,17,24,25`. Invalid grid/range values exercise every
lane across scalar-only, complete-block and reserved-tail layouts.

The script exhausts **6,561** states at `n=2,B=1`, checking every group
prefix and product count. An ordered-word polynomial certificate checks the
four-monomial identity without commutation. A separate factor-list model
checks **371 CSR-group prefixes in 96 cases**, covering unique-left,
unique-right and overlapping coordinates. These finite artifacts supplement
the arbitrary-n analytical proof and are not proof-assistant certificates.
The verifier checks the frozen hashes of both predecessor sources and uses
no timer, benchmark or vendor call. Historical novelty and measured speed
are separate questions not established by this proof or implementation.

A separate fresh-context, read-only audit reconstructed both group invariants
and checked arithmetic ranges, vector loads/tails, admission, CSR capacity,
extents and empty cases at source SHA-256
`7550990a7f6dd82fbd0e58580c9a6639eef7ffc07d841b56ae9cf189b394bde5`.
It found no substantive defect under the stated IEEE/storage contract.
One non-defect difference is rejection precedence for inputs with multiple
violations: native batching can detect a later invalid coefficient before
an earlier excessive defect, whereas scalar scanning interleaves the two.
Both reject; the API prescribes no precedence among simultaneous violations.
This audit is additional reasoning evidence, not a proof-assistant artifact.
