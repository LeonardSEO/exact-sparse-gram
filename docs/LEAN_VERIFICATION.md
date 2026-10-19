# Lean verification

Author: Leonard van Hemert. License: MIT.

The delivered project builds with Lean **4.24.0**, compiler commit
`797c613eb9b6d4ec95db23e3e00af9ac6657f24b`, and mathlib revision
`f897ebcf72cd16f89ab4577d0c826cd14afaafc7`. The complete transitive dependency
revisions are pinned in `lean/lake-manifest.json`. The checked compiler reports
`arm64-apple-darwin23.6.0`; this identifies the Lean executable, rather than the
C benchmark environment.

**Verdict:** the arbitrary-dimension real recurrence, both complete sparse-loop
prefix invariants, coefficient bounds, work identity, and reconstruction are
kernel checked. The complete rounded arithmetic reconstruction is kernel
checked under an explicit primitive rounding model. Its interpretation as IEEE
binary64 remains conditional on the separately stated `Binary64Bridge`; the
project does not formalize machine bit patterns, C execution, compiler semantics,
validation code, allocations, SIMD loads, or hardware instructions.

## Principal result and exact statement

The final declaration is `GramExact.binary64_gram_exact` in
`lean/GramExact/Bridge.lean`. Its fully elaborated type is printed in
`lean/evidence/types-and-axioms.txt` along with its transitive axioms. In ordinary
notation, it says:

- Given a rounding function satisfying the explicit `Binary64Bridge` contract,
  an integer exponent `-400 ≤ e ≤ 200`, arbitrary `N : ℕ`, arbitrary integer
  coefficient rows `A : ℕ → ℕ → ℤ`, `B ≥ 0`, `|A i k| ≤ B`, and
  `N B² ≤ 2^53`, every `j ≤ i` satisfies
  `roundedGram (2^e) round N A i j = Σ(k<N) (A i k·2^e)(A j k·2^e)`.
- `roundedGram` constructs its anchors with actual abstract rounded FMAs and
  reconstructs all later lower entries with the two endpoint FMAs followed by
  the complete first and second correction groups. No exactness conclusion
  about anchors, prefixes, updates, or Gram values is assumed.
- `N=0` returns zero explicitly. `N=n+1` uses `n` interior coordinates and the
  final coordinate `n`; `n=0` therefore includes the one-column boundary case.
  The row index is unbounded. A finite `m×N` matrix is represented by zero
  extension outside its dimensions; restricting the theorem to `i<m` gives
  its entire lower triangle. No support-size or recurrence-length bound is
  imposed.

The lower-level `rounded_gram_exact` has precisely the same reconstruction
conclusion in coefficient units, assuming only `RoundingModel (s*s)` for an
arbitrary real scale `s`. `scaled_dot` proves the equality between that target
and the real dot product; it is not an identification by naming convention.

## Statement mapping

| Paper/source obligation | Checked declaration(s) | Scope and conventions |
|---|---|---|
| First CSR prefix invariant (1) | `first_prefix`, `first_step` | Any finite prefix set inside `range n`; mixed left factors, old right factors; generic commutative ring, hence real scalars |
| Second CSR prefix invariant (2) | `second_prefix`, `second_step` | Any finite prefix set; new left factors, mixed right factors |
| Full support completion | `mixed_support`, `groups_join`, `groups_complete` | Support is defined by actual coefficient inequality; no trusted support metadata |
| Delete head, then insert tail | `boundary_identity`, `deleted_head`, `rounded_boundary_exact` | Source operation order retained; deletion separately bounded before insertion |
| Principal real recurrence | `recurrence` | Arbitrary real rows and arbitrary feature count `n+1` |
| Exact-anchor lower-triangle reconstruction | `reconstruction` | Induction follows `(i-1,j-1)` until a first-column anchor |
| Anchor-chain/reduction bound | `product_subset_bound`, `rounded_anchor_exact`, `rounded_anchor_merge` | Every distinct-coordinate subset; disjoint union at tree nodes |
| Both complete intermediate bounds | `first_prefix_bound`, `second_prefix_bound`, `first_prefix_cap`, `second_prefix_cap` | Integer coefficient magnitude `(n+1)B²`; includes overlaps and arbitrary within-group order |
| Residual subtraction | `residual_bound`, `residual_cap`, `rounded_residual_exact` | Admission derives the cap for `a-u`; it is not assumed separately |
| Exact rounded loops | `rounded_first_exact`, `rounded_second_exact`, `rounded_groups_exact` | Genuine folds of rounded FMAs; distinct support coordinates; first group fully finishes before the second |
| Full rounded recurrence | `rounded_update_exact` | Includes endpoint FMAs and both loops; exact predecessor premise is discharged by reconstruction |
| Constructed exact anchors | `computed_anchor_exact` | Direct finite FMA chain; `rounded_anchor_merge` also supports the C implementation's disjoint eight-chain reduction tree |
| All rounded lower entries | `rounded_reconstruction`, `rounded_gram_exact` | No exact intermediate premise; the latter also constructs its anchors |
| Real Gram target | `scaled_dot`, `binary64_gram_exact` | Integer-grid coefficient theorem transported to the real matrix target |
| Gram FMA work | `support_work`, `total_work` | `mn + 2D + mQ`, with `D = Σ(i<m)i`, `Q = Σ(i<m)q_i`, and `q_0=0`; no division convention hidden in `D` |

The work model counts one endpoint or correction FMA per listed operation and
`mn` anchor FMAs. It does not formalize C instruction counts, validation work,
allocation, wall-clock time, or workspace-byte arithmetic. `entry_work` merely
unfolds the per-entry count; the substantive triangular multiplicity proof is
`support_work`, used by `total_work`.

## Rounding model and the machine boundary

`GridRepresentable scale x` is the concrete proposition
`∃ z : ℤ, |z| ≤ 2^53 ∧ x = (z:ℝ)*scale`.
`RoundingModel scale` provides a function `round : ℝ → ℝ` and the primitive
property that it fixes those values. The coefficient bounds are independently
proved in `Bounds.lean`. `roundedFMA_coefficients` then derives exactness from
the polynomial identity for `a*b+c`; `roundedSub_coefficients` similarly derives
subtraction exactness. Thus the model does not postulate any algorithmic
conclusion. `identityRoundingModel` supplies a consistency witness, without
claiming to be a machine rounding implementation.

`Binary64SafeDyadic` explicitly describes the subset consisting of a bounded
integer coefficient times `2^e` with `-1022 ≤ e ≤ 970` (including zero).
`Binary64Bridge` assumes only that its interpreted machine rounding fixes this
subset. The input and doubled accumulator exponent containment, and
`(2^e)^2=2^(2e)`, are proved. Connecting that subset to a particular IEEE bit
representation and correctly rounded implementation is the disclosed external
mathematical/implementation correspondence. No imported floating-point library
or bit-level IEEE theorem is claimed.

Residual operands in `roundedUpdate` are coefficient differences scaled exactly;
`rounded_residual_exact` separately proves that the actual abstract rounded
subtraction computes those operands under the input-scale model. The same
`Binary64Bridge` supplies input-scale and accumulator-scale models. Negating an
input is exact at the mathematical representation boundary.

A fold uses `List.foldr`, so its operation order is the reverse of the supplied
list. Theorems accept **every** distinct-coordinate support enumeration; the C
increasing scan is represented by supplying its reversed coordinate list. The
canonical `roundedGram` definition chooses one such enumeration. The proof does
not identify that canonical enumeration with the C memory traversal; order
independence inside each group establishes the common exact result. The
completion of group one before group two is retained throughout.

## Execution and axiom audit

Successful commands, run from `lean/`:

```sh
lake build
lake env lean evidence/Inspect.lean
```

`lean/evidence/build.log` records the successful **944-job** project build.
`lean/evidence/types-and-axioms.txt` is the actual `#print` and `#print axioms`
output for the principal definitions and theorems. `lean/evidence/environment.json`
records exact versions, dependency revisions, source SHA-256 hashes, commands,
and results. All audited declaration closures use only the standard Lean axioms
`propext`, `Classical.choice`, and `Quot.sound` (some use a subset). There is no
`sorryAx`, custom mathematical axiom, native-evaluation proof, or kernel-bypass
setting in the delivered proof. `Binary64Bridge` is an ordinary explicit theorem
parameter and remains visible in the final theorem type; it is not a hidden
axiom. The exact-anchor premise of the intermediate reconstruction theorem is
discharged in `rounded_gram_exact` using constructed anchors.

The build reused a local cached **runtime/dependency environment only** through
an ignored `.lake/packages` symlink. No mathematical source from the unrelated
sibling project is imported. `.lake/` is ignored and excluded from the release.
For a fresh checkout with Lean's pinned toolchain available, Lake obtains the
locked dependency sources; `lake exe cache get` can obtain mathlib's matching
precompiled cache before `lake build`. No global toolchain upgrade is required.

The receipt verifies this source and these declarations in this environment.
The analytical C/source audit and empirical exact-entry tests remain separate
evidence for the delivered implementation; the Lean result is not a verified-C
or verified-compiler claim.
