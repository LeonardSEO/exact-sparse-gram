# Independent manuscript audit

Date: 19 October 2026. Author of the artifact: Leonard van Hemert.
License of this audit: CC BY 4.0.

## Scope and verdict

This bounded audit reads the manuscript, candidate C source, formal-verification
mapping and retained original/new numerical evidence. It runs no performance
measurement, changes no kernel, and makes no historical-priority finding.
PDF compilation and rendering inspection belong to the separate release check.

The real recurrence, complete first/second mixed-factor invariants and binary64
representability argument are coherent. The manuscript adequately distinguishes
the mathematical theorem, the restricted C admission API, the abstract Lean
rounding proof and the external IEEE/C/compiler correspondence. The empirical
claims inspected below agree with retained data and preserve losing cases.

One concrete arithmetic-count defect was reported to the release coordinator:
the reconstructed-Hankel comparator costs `mn+2D+(m+1)R` FMAs, not
`mn+2D+mR`. It corrects column-zero entries too: row i appears i+1 times on
the left and m-i times on the right. The standalone comparator document already
has the correct count. The manuscript now contains the corrected `(m+1)R` expression and explicitly
notes first-column correction. This finding is resolved. The assembled new
table was subsequently checked as recorded below; no remaining mathematical
or data-transcription blocker was found in this bounded manuscript audit.

## Mathematical and implementation correspondence

- Deletion leaves n-1 old-factor products; insertion restores n. Every first
  prefix uses new/old left factors against old right factors, and every second
  prefix uses new left factors against new/old right factors. Overlap and
  diagonal entries cause no exception. A complete first group is essential.
- `n B^2 <= 2^53` implies `B < 2^27` for n>=1, hence a residual coefficient
  has magnitude at most 2B<2^28. At n>=2 its product coefficient has magnitude
  at most 2B^2<=2^53; at n=1 no residual product exists. There is no hidden
  defect-count or recurrence-length capacity condition.
- At scale `2^(2e)`, nonzero integer-grid accumulator states are at least
  `2^-800` and no greater than `2^453`. These are normal and finite binary64
  values with representable significands. Correct rounding fixes them.
  Residuals at scale `2^e` are likewise normal when nonzero. Cancellation zero
  is a real equality claim, with signed-zero encoding explicitly excluded.
- The eight anchor chains contain disjoint coordinate sets; every reduction
  node is a sum over a disjoint union. The same nB^2 bound therefore applies
  to intermediate additions, not only to final anchors.
- The grouped witness independently evaluates to the printed integer
  18014398308155393 before rounding and the exact diagonal
  9007199120523265. The alternative grouped prefix exceeds the direct-dot
  capacity while the specified sequential prefix equals 67108864.
- The C correction loops match the pseudocode, including shifted old-right
  lookup at k+1 and unshifted new-left lookup at k. The NEON strict loop
  inequality reserves a scalar tail of 1..8 values, keeping shifted loads
  in the row. The theorem's unbounded dimensions are not attributed to C:
  n<=65536, max_q<=64 and supported exponent restrictions are disclosed.
- The integer quotient/square corrections establish the admission cap; sqrt
  is only a starting estimate. Integrality and reconstruction are both tested.
  Empty m/n and n=1 branches, partial output on failure, disjoint objects,
  and absence of automatic fallback are separately disclosed.
- The candidate FMA count mn+2D+mQ and requested scratch count agree with
  source. Scaling, reduction additions, validation, reflection, and source
  traffic are not misrepresented as part of that FMA count or measured DRAM.

The Lean statement mapping accurately exposes the ordinary theorem parameter
`Binary64Bridge`, canonical chain versus eight-chain correspondence, reversed
foldr enumeration, separate residual-subtraction lemma, and unverified machine
execution. This audit does not replace the independent Lean build receipt.

## Retained empirical evidence and generated assets

The original raw JSONL has 230 records, and the summary CSV has 230 records.
Every method median (0/4/9/11/12/13) and both vendor/candidate paired-ratio
medians were independently recomputed from all raw samples and matched exactly.
There are 150 paired-median wins against both vendors and 80 losses to SYRK.
The original JSONL SHA-256 is
`a3e585c5f496309574ebdc1630e278d35a20e2b759b5e355e25edcdead86b948`,
matching `evidence/original/MANIFEST.json`'s original file record. No original
raw measurement was replaced by the specialized campaign.

The original `selected_results.tex` and `speedup.tex` were reproduced in memory
from `analysis/make_paper_assets.py`; both outputs matched their delivered
contents byte for byte. Their selections include losing cases. The figure
shows minimum/maximum run-level ranges, not bootstrap confidence intervals.
The bootstrap code resamples contiguous balanced blocks; the manuscript
correctly calls its intervals conditional within-run and disclaims independent
cross-machine or simultaneous coverage. Repeated seeds/scales/output layouts
are not advertised as independent input draws.

The full-entry receipt totals independently sum to 1,263,744 candidate entries,
1,263,744 SYRK entries and 2,523,136 GEMM entries: 5,050,624 total. The manuscript
correctly allocates the checks across methods and distinguishes these untimed
all-entry checks from sampled preflight checks on timed matrices. Vendor zero
error on those particular cases is preserved.

The separate specialized CSV contains 20 records, five seeds in each of four
family/distribution groups. Reconstructed-Hankel mismatch counts are 518 for
windows and 213,176..234,577 for innovations; Q is 955..1022. These directly
support the distinction between accumulated sparse transitions and this
specific reconstructed-base comparator, without establishing a universal
structured-algorithm lower bound. The pre-freeze harness seed fix is disclosed;
the seed/fingerprint-bound frozen campaign is the relevant new evidence.

## Primary bibliography checks

The following primary records were read on the audit date:

- [CUNY institutional report record](https://academicworks.cuny.edu/gc_cs_tr/308/):
  Pan, Murphy and Rosholt; TR-2008003; title and 2008 year match. The institutional
  PDF was not retrievable in this pass. Section/theorem matching is supported
  by the separately retained source-comparison note rather than claimed as a
  fresh primary-PDF reinspection here.
- [FFLAS/FFPACK author manuscript record](https://arxiv.org/abs/cs/0601133v3):
  Dumas, Giorgi and Pernet; version 3 dated 14 January 2009; ACM TOMS 35(3),
  article 19 (2009). Manuscript bibliography metadata match.
- [STOMP author-hosted paper](https://mcyeh.github.io/paper/2016_icdm_exploiting_novel_algorithm.pdf):
  title and eight authors match; Section IV-A Eq. (4) is the sliding-dot
  remove/add recurrence, and the output is a matrix profile rather than a
  materialized Gram matrix. The manuscript's attribution is appropriately scoped.
- [Gilda version-specific Zenodo record](https://zenodo.org/records/20599316):
  title, author, DOI, v1 publication 8 June 2026 and modification 29 June 2026
  match. The manuscript makes no independently verified formal-proof claim
  for that source.

SCAMP's DOI and author-hosted PDF were not retrievable by the browser in this
pass. Its bibliographic and section-level support therefore remains the earlier
retained comparison's check, not a new live verification. This is an access
limitation, not evidence against the cited work. No exhaustive priority hunt or
claim of first publication was attempted. Dense SYRK/GEMM comparisons and the
newly implemented focused comparator are explicitly separated from published
structured state of the art.

## Release follow-up

The assembled `figures/specialized_results.tex` was reproduced in memory
from its delivered generator and matches byte for byte. All four table rows
were independently compared with the five corresponding retained CSV rows:
R ranges, paired medians, minima and maxima match at printed precision. The
20/20 wins against SYRK, GEMM and this specific H+E comparator are supported;
all three conditional interval lower endpoints also exceed one in every run,
though the table correctly labels its displayed brackets as observed ranges.
The prose rounds H+E ratios to about 1.7--1.8 (windows) and 42--46
(innovations), and explicitly leaves better Hankel reconstruction and other
published structured implementations open. This supports the narrow
comparison and does not establish structured state of the art.

Verdict: bounded manuscript content and numerical-asset review passes after
the resolved comparator-count correction. Artifact-availability URLs, immutable
publication receipts, checksums and rendered PDF quality remain separate
release-coordinator checks. The SCAMP primary-source re-fetch limitation above
remains disclosed; it is not relabeled as freshly verified.
