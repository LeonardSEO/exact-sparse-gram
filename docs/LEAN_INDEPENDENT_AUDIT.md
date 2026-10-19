# Independent Lean assurance audit

Audit date: 19 October 2026. Author of the mathematical artifact: Leonard van
Hemert. License: MIT.

## Verdict and scope

No blocking mathematical statement mismatch, circular rounding premise, or
unfinished proof was found in the delivered source and supplied elaborated
declaration/axiom receipt. The artifact establishes the arbitrary-dimension
abstract rounded lower-triangle theorem under its explicit primitive rounding
model. Its binary64 interpretation remains conditional on `Binary64Bridge`.
That condition is a representation/primitive-rounding correspondence, not an
assumption that any Gram entry, anchor, loop, or recurrence is exact.

This independent pass inspected every delivered Lean source file, the principal
elaborated types and transitive axiom reports, the pinned build configuration,
`docs/LEAN_VERIFICATION.md`, the sibling `gram_research/REPORT_EXACT.md`, and the
released `src/lattice_sequential.c`. It also compared the manuscript's theorems,
formal mapping table and trust-boundary discussion in `paper.tex`. No mismatch
was found. During benchmark quiet time it ran no Lean build, numerical test, or
CPU benchmark. Once that lock was released, it rebuilt every project module in
a fresh isolated project and reran the complete declaration/axiom inspection.

## Intended statement and actual target

The research claim is exact Gram reconstruction for integer-grid rows
`A[i,k] = X[i,k] 2^e`, `|X[i,k]| <= B`, `N B^2 <= 2^53`, and
`-400 <= e <= 200`, using head deletion, tail insertion, the entire first sparse
factor-replacement group, then the entire second group. Row count and recurrence
length must not introduce a further error or admission condition.

The inspected target is `GramExact.binary64_gram_exact`. Its actual type has
arbitrary `N : Nat`, `B : Int`, and `A : Nat -> Nat -> Int`; explicit coefficient
bounds and admission inequality; arbitrary row indices `i j` with `j <= i`;
and a real dot-product conclusion. The only extra interpretation parameter is
`I : Binary64Bridge`. There is no support-size, row-count, restart, chain-length,
or exact-intermediate hypothesis. The stronger global bound on the infinite row
function does not exclude any finite admitted matrix: zero extension satisfies
it, and only the requested finite coordinates are used.

`scaled_dot` supplies a proved equality to the real target; `dot` is an ordinary
finite sum of real products. This is not an exactness conclusion obtained by
renaming the algorithm as its specification.

## Detailed assurance checks

| Obligation | Independent finding |
|---|---|
| Primitive rounding assumption | `GridRepresentable s x` means an integer coefficient of magnitude at most `2^53` times `s`. `RoundingModel.fixes` fixes exactly that explicit subset. It contains no reference to the Gram algorithm. |
| External binary64 binding | `Binary64SafeDyadic` has exponent interval `[-1022,970]`; bounded coefficients describe finite representable dyadics, including zero. `Binary64Bridge.fixes_safe` is the disclosed external binding. Input and doubled exponents lie inside this interval. No bit-level IEEE proof is claimed. |
| Nonvacuity | The identity rounding function provides `identityRoundingModel`. It is correctly identified as a consistency witness, not a machine implementation. |
| First and second prefixes | `first_prefix` and `second_prefix` independently rewrite the accumulator into a sum with one actual old/new bounded factor at every coordinate. The bounds therefore have `N` products, including overlapping supports. |
| Endpoint order | `rounded_boundary_exact` proves the deletion result bounded before proving insertion; it does not use the final-dot bound to skip the deletion intermediate. |
| Support completion | `support` is the actual inequality filter over the interior coordinate range. `mixed_support`, `groups_join`, and `groups_complete` justify omission of equal coordinates. |
| Enumeration | Rounded-loop lemmas quantify over every duplicate-free list with the specified support. `foldr` executes the reverse order. The increasing C scan is covered by supplying its reversal; within-group order is unrestricted, but group order is preserved. |
| Anchors | `rounded_anchor_exact` proves every distinct-coordinate chain. `rounded_anchor_merge` proves a disjoint-union reduction node. Thus the eight C chains and their disjoint reduction tree have the needed local mathematical obligations. The final definition itself uses a direct chain. |
| Reconstruction | `rounded_reconstruction` has an intermediate exact-anchor premise, then `rounded_gram_exact` discharges it using `computed_anchor_exact`. The exact-predecessor hypothesis in the single-update lemma is discharged by induction. |
| Upper fallback | `roundedReconstruct` contains an exact-dot branch for `(0,j+1)`. The lower-triangle condition makes that branch unreachable in the theorem's predecessor chain; it cannot secretly provide exact interior or anchor values. |
| Empty dimensions | `roundedGram` explicitly returns zero for `N=0`. At `N=1`, the interior range is empty, so both support lists are empty and only endpoint operations remain. A finite matrix with no rows has no requested lower entries. |
| Residuals | `residual_cap` derives the subtraction bound from admission. `rounded_residual_exact` proves the actual abstract rounded subtraction at input scale. The final loop definition uses the resulting exact coefficient differences, as documented. |
| Work | `support_work` proves the triangular multiplicity identity, rather than merely unfolding a chosen result. `total_work` derives `mN + 2D + mQ` for the listed nonempty-dimension schedule, with `q 0 = 0`. |

The C/source connection remains analytical rather than a verified execution
semantics: the formal definition substitutes exact coefficient differences for
stored residuals, chooses a canonical support enumeration, and uses direct-chain
anchors. Separate checked lemmas justify residual subtraction, arbitrary support
order, and disjoint anchor-tree merges. No theorem identifies C memory, SIMD
instructions, allocations, validators, or compiler output with this definition.
These limits are already disclosed in `LEAN_VERIFICATION.md`; they do not weaken
the theorem into a circular exactness assumption.

The work theorem is an identity for an explicitly defined operation-count model.
It is not a theorem about measured instructions or about the zero-column early
return: for `N=0`, the C algorithm writes zeros and performs no Gram FMAs, while
the symbolic nonempty-schedule count would still include endpoint terms. The
research report explicitly restricts that FMA formula to nonempty dimensions.

## Build and axiom evidence

The project pins Lean `4.24.0` (compiler commit
`797c613eb9b6d4ec95db23e3e00af9ac6657f24b`) and mathlib revision
`f897ebcf72cd16f89ab4577d0c826cd14afaafc7`; the manifest pins the transitive
dependencies. The project imports only its four delivered modules and mathlib.
No sibling mathematical project is imported.

The supplied `lean/evidence/build.log` ends with
`Build completed successfully (944 jobs).` Its environment receipt records exit
status 0 for `lake build` and `lake env lean evidence/Inspect.lean`.
The inspection source prints the actual principal definitions, statements, and
transitive axioms. The saved 53 axiom reports contain only the accepted standard
axioms `propext`, `Classical.choice`, and `Quot.sound` (or a subset).

In particular the saved reports for `rounded_gram_exact` and
`binary64_gram_exact` list exactly those three axioms. Source inspection found
no `sorry`, custom mathematical `axiom`, `native_decide`, unsafe declaration,
or kernel-bypass option. `pp.proofs false` affects printing only.
`Binary64Bridge` is visible as an ordinary theorem parameter, rather than being
hidden from `#print axioms` by an implicit global postulate.

## Independent fresh-project replay

The new receipts are in `lean/evidence/independent/`. `setup.json` records the
temporary project path, dependency-cache path, and SHA-256 of each copied input.
Only the exact project `.lean` sources, `lean-toolchain`, `lakefile.toml`,
`lake-manifest.json`, and `evidence/Inspect.lean` were copied. No project build
directory, `.olean`, or other project proof artifact was copied. The new
`.lake/packages` symlink reused only the existing pinned dependency/runtime
cache, without rebuilding all mathlib.

Using `/Users/leonard/.elan/bin` on `PATH`, commands from the isolated project
were:

```sh
lean --version
lake --version
lake build
lake env lean evidence/Inspect.lean
```

All four commands exited 0. Lean reported `4.24.0`, compiler commit
`797c613eb9b6d4ec95db23e3e00af9ac6657f24b`, and Lake reported
`5.0.0-src+797c613`. The fresh build completed successfully, compiling all five
project modules, with no warning or error diagnostics. The complete independent
type/axiom output is byte-for-byte identical to the original supplied output:
all 53 reports, including both final theorem types and their transitive closures,
match exactly. The exact final principal closure is:

```text
'GramExact.rounded_gram_exact' depends on axioms: [propext, Classical.choice, Quot.sound]
'GramExact.binary64_gram_exact' depends on axioms: [propext, Classical.choice, Quot.sound]
```

`verification.json` records equality of every source byte/hash before and after
replay, all dependency checkout revisions matching the locked manifest, the
newly produced five project `.olean` paths, and every inspected axiom closure.
The logs and `commands.json` retain the commands, statuses, timings and complete
diagnostics. This is a fresh project-source replay using cached dependencies;
it is not an independent rebuild of Lean or all dependency sources.

Analytical status: the inspected real recurrence, intermediate bounds, and
rounding-to-grid induction are sound within the stated abstraction. Formal
status: independently verified in the pinned environment for the named abstract
rounded targets, with binary64/C correspondence conditional as stated. No proof
repair is required by this audit.
