# Arithmetic, traffic, and storage of the frozen sequential lattice algorithm

Read-only source analysis, 2026-10-19. The analyzed file is
`src/lattice_sequential.c`, SHA-256
`7550990a7f6dd82fbd0e58580c9a6639eef7ffc07d841b56ae9cf189b394bde5`.
The hash was read and confirmed. No source changes, compilation, tests, or
timing were performed for this analysis. Numbers below are mathematical or
source-access counts; none are measured hardware traffic counters.

## Notation and exact arithmetic count

For a successful call with m,n>0, let

- `P=m(m+1)/2`, the number of lower-triangle outputs;
- `D=P-m=m(m-1)/2`, the number of recurrence entries after first-column anchors;
- `c=min(max_q,n-1)` and `C=(m-1)c`, the reserved residual capacity;
- `q_i` be the actual number of shifted defects in row i, with `q_0=0`;
- `Q=sum_{i=1}^{m-1}q_i`, so `0<=Q<=C`.

Each first-column anchor performs exactly n scalar FMAs across its eight chains,
for mn in total. Each recurrence entry `(i,j)`, `1<=j<=i<m`, performs exactly
`2+q_i+q_j` FMAs: two endpoint updates, q_i left-factor substitutions, and q_j
right-factor substitutions. There is no support intersection or merge term.
Consequently the exact number of Gram scalar-product/FMA evaluations is

`F = mn + 2D + sum_{i=1}^{m-1} sum_{j=1}^{i}(q_i+q_j)`

`  = mn + 2(P-m) + mQ`.

Indeed, residual row t contributes t times in the left group (i=t), and m-t
times in the right group (j=t), hence exactly m appearances. This also counts
both substitutions on a diagonal entry; they cannot be discarded as duplicates.

For a uniform positive per-row defect count q, `Q=q(m-1)` and
`F=mn+(q+1)m(m-1)`. For a capacity cap q with smaller actual support, this is
an upper bound rather than the actual count. Conventional lower-triangle SYRK
uses nP scalar product terms; actual optimized instruction schedules, redundant
diagonal-tile work, and hardware throughput are separate matters.

F is **not the complete floating-point instruction count**. On a successful
call, input lattice admission computes `a*inverse` and then `x*scale` for every
input: another 2mn scalar multiplications, implemented in vector lanes where
available. There are also range/equality/integrality predicates, conversions to
comparison masks, and branches. The eight-chain anchor tree performs seven
scalar additions per row, 7m total at the source-algorithm level, including
zero chains for short n. Each extracted defect performs one exact subtraction,
Q in total. The coefficient-bound calculation uses a square-root estimate and
verified integer square comparisons once per call. The integer indexing, loop
control, allocation/free costs, and memory accesses are not represented by F.

Counting one FMA as two conventional FLOPs gives 2F for the FMA body, but that
notation must not hide the 2mn admission multiplications or other preparation.

## Actual allocations and output layout

On the stated 64-bit ABI, `struct lattice_residual {size_t k; double d;}` is
16 bytes: an eight-byte coordinate and eight-byte exact residual. When C>0 the
function allocates exactly:

`residual array: 16C bytes`, `CSR offsets: 8(m+1) bytes`.

Thus requested heap workspace is

`W = 16(m-1)c + 8(m+1)` bytes if C>0, and `W=0` if C=0.

There are two allocations and two frees for positive capacity. Actual support
Q does not shrink the reserved capacity. Unused residual records are not
initialized. There is no generator copy, packed dense input, norm array, error
array, anchor array, per-entry certificate, or second Gram buffer. Anchors are
written directly into G's first column. The q8 chains and four SIMD accumulators
use O(1) stack/register state; actual compiler frame size is not claimed here.

For contiguous input/output, caller-owned dense allocations are `8mn` bytes for
A and `8m^2` bytes for G. With strides, the last addressed footprints are
`8((m-1)lda+n)` and `8((m-1)ldg+m)` bytes; allocated padding can be larger.
Only the lower triangle is written: `8P` result bytes. The full dense G
allocation must still be included when comparing this API with SYRK. A packed
triangle interface would be a different contract.

## Preparation: validation, extraction, and q8 anchors

The NEON preparation loop handles eight features at a time and deliberately
uses `k+8<n`. Let `b=floor((n-1)/8)` and `t=n-8b`, so `1<=t<=8`. Every row has b
vector blocks followed by a scalar tail of t values. For example, n=4096 has
511 blocks and an eight-value tail; n=8192 has 1023 blocks and an eight-value
tail. The strict inequality keeps every shifted predecessor load in bounds.

Within a successful vector block:

- Eight current-row A values are loaded once into four vectors. The same
  register values feed range/grid validation, shifted equality comparisons,
  and the left factors of the anchor FMAs.
- Eight first-row anchor factors are loaded.
- For i>0, eight shifted predecessor factors are loaded and compared. Across
  the complete row, these predecessor comparisons number n-1, since the final
  feature has no shifted predecessor.
- If any shifted comparison differs, scalar extraction checks all eight
  positions. Its source expressions may reread eight current and eight shifted
  values even when only one lane changed. Thus one changed block exposes up
  to 128 bytes of additional source operand reads. Compilers can reuse values;
  this is not an asserted extra cache miss or DRAM transfer.

Let H be the number of changed vector blocks over all nonfirst rows. Every such
block contains an actual defect, so `H<=Q`, as well as `H<=(m-1)b`. Extra
changed-block source reads are at most `128H` bytes under this reread model.
For max_q=2, at most two vector blocks per nonfirst row require extraction.
Scalar tail expressions likewise allow local repeated reads of the same A
element; ordinary register reuse can eliminate them. Invalid block replay is
not part of the successful-call counts given here.

A useful **reuse-aware operand-role model**, counting current input, first-row
anchor factors, and shifted predecessor factors separately, is

`L_pre = mn + mn + (m-1)(n-1)` doubles,

or `8L_pre` bytes, before changed-block rereads. Roles can alias, especially for
row zero, and scalar-tail loads can be commoned. Therefore this is a transparent
logical-access model, not an exact retired-load count. It also is not DRAM
traffic: repeatedly used first and adjacent rows may stay in cache. Conversely,
large rows, cache conflicts, or cold input can increase physical transfers.

The preparation writes Q residual records, `16Q` bytes, and, for positive
capacity, m+1 CSR offsets, `8(m+1)` bytes. It writes m first-column anchors,
`8m` bytes. No dense input conversion or external preparation is excluded from
the API timing; validation and residual extraction occur on every call.

## Recurrence operand traffic and reuse

For each of the D recurrence entries, the source uses:

- One previous-row G value and one final G store.
- Four endpoint A factors for the two endpoint FMAs. The left endpoints are
  constant across a fixed output row and may be hoisted/retained.
- For each of q_i+q_j corrections, one 16-byte CSR residual record and one
  eight-byte dense A factor at its stored coordinate. Correction factors are
  gathers from dense A, rather than a packed row transform.
- CSR interval endpoints for row i and row j. Four offsets per output is a
  direct loop-endpoint operand model; compilers can hoist the two row-i bounds,
  and retain loop bounds in registers. It is not an exact machine load count.

The following logical components follow from the source, prior to reuse:

| Component | Logical bytes |
|---|---:|
| G stores, including m anchors | `8P` |
| Previous-row G reads | `8D` |
| Four endpoint A operands | `32D` |
| Dense A correction-factor operands | `8mQ` |
| CSR residual-record operands | `16mQ` |
| Four CSR loop endpoints, when rows exists | `32D` |

The first two sum to `8(P+D)=8m^2` bytes, despite only half the matrix being
written: the recurrence also reads predecessors. These reads can often hit
the immediately preceding G row in cache. No initial G zeroing/read is needed.
For c=0, both correction and CSR-offset components vanish.

The CSR lists occupy O(mc) bytes and are reused many times; multiplying record
size by mQ is an operand-demand count, not a claim that each record is fetched
from DRAM m times. Similarly, dense correction factors may hit cached rows or
may cause noncontiguous cache-line fetches. The exact pattern depends on support
coordinates, strides, cache geometry, and replacement. Scalar dependency chains
within each output also constrain FMA throughput, irrespective of cache traffic.

If the caller requires a full symmetric G, mirroring the strict lower triangle
adds D value reads and D stores: `16D` logical bytes and no floating arithmetic.
That same materialization must be included for lower-only SYRK. GEMM directly
produces a full result and has a different arithmetic/output schedule.

## Representative sizes and counts

The table takes max_q=2 and q_i=2 for every later row, hence Q=2(m-1).
For smaller actual support, the allocation stays the same while F and correction
operand counts decrease. These are analytic counts, not benchmark results.

| Quantity | m=512, n=4096 | m=768, n=8192 |
|---|---:|---:|
| P lower outputs | 131,328 | 295,296 |
| D recurrence outputs | 130,816 | 294,528 |
| Q actual defects in this example | 1,022 | 1,534 |
| Gram FMAs F | 2,882,048 | 8,058,624 |
| Additional admission multiplies | 4,194,304 | 12,582,912 |
| q8 reduction additions | 3,584 | 5,376 |
| Conventional SYRK product terms nP | 537,919,488 | 2,419,064,832 |
| Dense A allocation | 16,777,216 B = 16 MiB | 50,331,648 B = 48 MiB |
| Dense G allocation | 2,097,152 B = 2 MiB | 4,718,592 B = 4.5 MiB |
| Required lower result stores | 1,050,624 B | 2,362,368 B |
| Sequential heap workspace | 20,456 B | 30,696 B |
| Preparation operand-role model | 50,294,792 B | 150,923,272 B |
| Changed-block reread bound 128Q | 130,816 B | 196,352 B |
| G logical reads plus stores | 2,097,152 B | 4,718,592 B |
| Endpoint A operand bytes | 4,186,112 B | 9,424,896 B |
| Correction A operand bytes | 4,186,112 B | 9,424,896 B |
| CSR record operand bytes | 8,372,224 B | 18,849,792 B |
| Four-offset operand model | 4,186,112 B | 9,424,896 B |
| Additional full-mirror bytes | 2,093,056 B | 4,712,448 B |

## Empty, zero-defect, and narrow boundaries

- m=0 returns immediately before pointer/domain checks: no accesses, arithmetic,
  or allocation.
- For m>0,n=0, output extents and scale/q parameter caps are checked, then the
  P lower entries are set to zero. There is no A access, coefficient-bound
  division, residual allocation, or FMA. The n>0 formulas must not be extended
  blindly to this case.
- For c=0, residual capacity is zero and both r and rows are null. Accepted
  n>1 input has no shifted defects. The arithmetic is `F=mn+2D`; the full
  validation/equality scan still occurs.
- For n=1, c=0 and Q=0 for any admitted max_q. Any bounded one-column A is
  admissible. The endpoint recurrence cancels the previous one-term product,
  then inserts the new one. The actual implementation uses `F=m+2D=m^2` FMAs,
  rather than P direct products. This boundary is correct but is not a promised
  fast regime.
- For m=1,n>0, capacity is zero irrespective of c. Only the n-term anchor and
  its admission checks execute; there are no recurrence entries.

## Fair comparison with the earlier TwoDiff certificate implementation

Both `numerical_sparse.c` and `numerical_sparse_fast.c` reserve residual records
`{size_t k; double hi,lo;}` of 24 bytes on the same ABI. They allocate CSR offsets
even at zero capacity, plus m norm doubles and m error doubles. For nonempty
input their requested scratch is

`W_old = 24C + 8(m+1) + 16m` bytes.

For C>0, the new source saves `8C+16m` bytes. At the two representative sizes,
the older scratch is 36,824 B and 55,256 B, versus 20,456 B and 30,696 B. At
C=0 the new source allocates nothing, whereas that sparse TwoDiff source still
allocates `24m+8` bytes. The benchmark's separate zero-defect Hankel guard may
have a different generator/norm/error allocation and must not be conflated with
this sparse-source formula.

This is not the removal of storage from an equally broad arbitrary-double
algorithm without compensation. TwoDiff supports a wider floating input class
and supplies an approximate norm-scaled certificate with possible repairs.
The new source enforces an explicit integral lattice/range admission rule and
uses a proof of exact arithmetic, making norm/error arrays unnecessary in that
admitted class. The comparison must use inputs accepted by both algorithms and
the same required external accuracy; the new algorithm's exactness is stronger
than the earlier tolerance criterion.

On those lattice inputs, TwoDiff's low residual component is zero. Its core
recurrence can then have the **same** `2+q_i+q_j` FMA count as the new source.
The differences include the earlier three-sum norm/anchor/absolute prepass
(3mn FMAs), TwoDiff extraction work and 24-byte records, error-magnitude
accumulation, norm targets, bound updates, and possible repairs. The new source
has only mn anchor FMAs, but also pays its lattice-admission multiplications.
Thus any speedup over the older certificate code should not be explained solely
by the nP-versus-F comparison with conventional SYRK.

## Limits of performance predictions

The large reduction in Gram product terms is a valid arithmetic property. It
does not establish the practical speedup. The new source still scans all dense
A entries on every call, performs lattice/grid predicates, branches on sparse
blocks, allocates/free CSR buffers, traverses sparse lists, gathers dense factors,
and runs short dependent scalar FMA chains. Optimized SYRK instead has regular
packed data, long vector/matrix kernels, and much higher achievable arithmetic
throughput. Cache retention and complete-output mirroring can change the ratio.

Only the controlled final timings can establish the practical regime. This note
provides no bandwidth, cache-miss, energy, or DRAM counter measurements, and the
logical-access tables must not be presented as such.
