# Independent verification index

Run `./build.sh && ./verify.sh`. No Python dependency installation is needed.
The Round1 transcript is `results/final_verification.log`; the successful
integrated suite, including all three exact-lattice members, is
`results/round2_verification.log`.

| Check | Artifact | Scope |
|---|---|---|
| Dense kernel exact dyadic accuracy | verification/kernel_verify.py; results/kernel_exact.json | 560 cases,149600 entries; odd dimensions,paddedstrides,zeros,dynamicranges,near dependence,upper untouched |
| Independent skinny prototype | verification/communication_verify.py | 132 shapes,17435 exactinteger entries,memory guards |
| Earlier scalar formula | verification/previous_identity_verify.py | 9 zero polynomial residuals,720 Fraction cases,exact block obstruction,cancellation; saved timing ratio recalculated |
| Ordered sparse construction | verification/algebra_verify.py; results/algebra_exact.json | 2112 formal identities,374400 anchored formal outputs,1200 rational cases,1440 compiledcalls/94224integer entries; overflow rejection |
| Scalar certified Hankel | verification/numerical_verify.py | 1349 exactFraction entries;21repairs exercised |
| Scalar certified sparse | verification/numerical_sparse_verify.py | 1038 exactFraction entries;exactlow residuals/21repairs |
| Paired SIMD Hankel | verification/communication_guard_verify.py | 1469 exactFraction entries;14signals/21repairs |
| Final fused Hankel | verification/communication_guard_fast_verify.py | 1469 exactFraction entries;14signals/21repairs |
| Final fused sparse | verification/numerical_sparse_fast_verify.py | 1308 exactFraction entries;25cases including all short SIMD tails andq64/21repairs |
| Published ordered coefficients | verification/flip_verify_certificates.py | All8 suppliedGT schemes, rational coefficient arithmetic; no upstream search imports |
| Published v2Figure4 | verification/flip_extract_paper_scheme.py | 2560 ordered coefficients, characteristiczero,34products/12partial-support/8trueGram |
| Executed mixed recursion | verification/flip_test_reference.py | 18 independentX/Y and18Gram heldoutcases, arbitrary rectangular padding/depth |
| Full output completeness | results/full_output_verification.json | All57finalfull cases checked inbothtriangles, against mirroredroundedSYRK at4e-12; separate from exact numerical proof |
| Independent mathematical challenge | notes/audit_guard.md; notes/audit_fast.md | Sourcehashes, roundinghelper proofs, SIMDindexing,cachedseed invariant, full fixedtol totality |
| Benchmark source audit | notes/audit_benchmark.md | Contracts,clockboundaries,balancedorders,generators,statistical/reporting limits |
| Grouped exact-lattice member | verification/lattice_verify.py; verification/lattice_results.json | Native/scalar each34394exact entries, coefficient/grid/range/stride checks |
| Merged-support substitution member | verification/substitute_verify.py; verification/substitute_results.json | Native/scalar each1562exact entries;6561finite prefix states and strict grouped-order witness |
| Sequential exact-lattice member | verification/sequential_verify.py; verification/sequential_results.json | Native/scalar each2044exact entries;6561finite invariant states, ordered polynomial,371CSR-prefix checks; q64, cap/exponent endpoints, every vector lane/tail |
| Exact arithmetic and access proof | notes/lattice_sequential.md; notes/sequential_traffic.md | Arbitrary dimension/prefix proof, no-merge ordinary-dot bound, CSR/allocation/NEON safety, honest logical traffic |
| Round2 benchmark audit | notes/audit_lattice_benchmark.md | Six-method balance, timer boundary, fixed230jobs; poisoned/all-finite output plus exact sampled oracle; limitations |
| Theorem-level contribution challenge | notes/priority_displacement.md; notes/contribution_sequential.md | Algebra subsumption by Pan; known integer-FP principle; narrow surviving schedule/guarantee, no historical-first inference |
| Independent full-entry integer oracle | verification/algebra_lattice_full_entry_verify.py; verification/algebra_lattice_full_entry_results.json | Four large matrices, each contiguous/padded:1,263,744candidate entries, same SYRK count,2,523,136GEMM entries; all5,050,624exactly equal; no poison/input/stride failures |
| Fresh final source/data audit | notes/final_exact_audit.md | Independent prefix/CSR/NEON/extent reconstruction; all230manifest jobs and41,400timings checked; source/binary hashes stable; bounded prior comparison |

The Fraction oracle converts each actual binary64 input exactly, forms exact
products/sums, and compares squared error with τ² times exact squared row norms.
This avoids using rounded reference norms. Those tests are independent of the
implementation's certificate calculation; source-level proofs establish the
universal statements beyond these finite checks. No Lean formalization is claimed.

The timing harness's4e-12 check against roundedSYRK is a preflight consistency
check, **not** the exact2e-12 accuracy gate. The2e-12 claim rests on the universal
proof and independent exactFraction checks. Freshanchor totality is separately
proved in notes/anchor_totality.md and independently audited.
