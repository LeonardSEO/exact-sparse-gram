import GramExact.Bounds
import GramExact.Arithmetic

namespace GramExact
noncomputable section
open scoped BigOperators

/-- Coefficient scaling used throughout the accumulator proof. -/
def scaled (s : ℝ) (z : ℤ) : ℝ := (z : ℝ) * (s*s)

/-- Actual rounded first-loop schedule. `foldr` executes the listed coordinates
from right to left; every CSR enumeration is represented by its reversed list. -/
def roundedFirst (s : ℝ) (round : ℝ → ℝ) (a u v : ℕ → ℤ)
    (initial : ℝ) (ks : List ℕ) : ℝ :=
  ks.foldr (fun k c => roundedFMA round (((a k-u k : ℤ) : ℝ)*s)
    ((v k : ℝ)*s) c) initial

def roundedSecond (s : ℝ) (round : ℝ → ℝ) (a b v : ℕ → ℤ)
    (initial : ℝ) (ks : List ℕ) : ℝ :=
  ks.foldr (fun k c => roundedFMA round ((a k : ℝ)*s)
    (((b k-v k : ℤ) : ℝ)*s) c) initial

/-- Every first-loop rounded operation is exact, for arbitrary distinct lists. -/
theorem rounded_first_exact (s : ℝ) (M : RoundingModel (s*s))
    (n : ℕ) (x y B : ℤ) (a u v : ℕ → ℤ)
    (hB : 0 ≤ B) (hx : |x| ≤ B) (hy : |y| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hu : ∀ k, |u k| ≤ B) (hv : ∀ k, |v k| ≤ B)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53)
    (ks : List ℕ) (hnd : ks.Nodup) (hks : ks.toFinset ⊆ Finset.range n) :
    roundedFirst s M.round a u v (scaled s (firstState n (x*y) a u v ∅)) ks =
      scaled s (firstState n (x*y) a u v ks.toFinset) := by
  induction ks with
  | nil => rfl
  | cons k ks ih =>
    obtain ⟨hk, hnd⟩ := List.nodup_cons.mp hnd
    have hsub : ks.toFinset ⊆ Finset.range n := by
      intro j hj
      exact hks (by simp only [List.toFinset_cons]; exact Finset.mem_insert_of_mem hj)
    change roundedFMA M.round _ _ (roundedFirst s M.round a u v _ ks) = _
    rw [ih hnd hsub]
    have hstep := first_step n (x*y) a u v ks.toFinset k (by simpa using hk)
    have hcap := first_prefix_cap n x y B a u v (insert k ks.toFinset)
      (by simpa using hks) hB hx hy ha hu hv hadmit
    rw [hstep] at hcap
    have hc : |(a k-u k)*v k + firstState n (x*y) a u v ks.toFinset| ≤ coefficientCap := by
      unfold coefficientCap
      simpa [add_comm] using hcap
    unfold scaled
    rw [roundedFMA_coefficients M (a k-u k) (v k) _ hc]
    simp only [List.toFinset_cons]
    rw [hstep]
    congr 2
    ring

/-- Every second-loop rounded operation is exact, independently of overlap. -/
theorem rounded_second_exact (s : ℝ) (M : RoundingModel (s*s))
    (n : ℕ) (x y B : ℤ) (a b v : ℕ → ℤ)
    (hB : 0 ≤ B) (hx : |x| ≤ B) (hy : |y| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hb : ∀ k, |b k| ≤ B) (hv : ∀ k, |v k| ≤ B)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53)
    (ks : List ℕ) (hnd : ks.Nodup) (hks : ks.toFinset ⊆ Finset.range n) :
    roundedSecond s M.round a b v (scaled s (secondState n (x*y) a b v ∅)) ks =
      scaled s (secondState n (x*y) a b v ks.toFinset) := by
  induction ks with
  | nil => rfl
  | cons k ks ih =>
    obtain ⟨hk, hnd⟩ := List.nodup_cons.mp hnd
    have hsub : ks.toFinset ⊆ Finset.range n := by
      intro j hj
      exact hks (by simp only [List.toFinset_cons]; exact Finset.mem_insert_of_mem hj)
    change roundedFMA M.round _ _ (roundedSecond s M.round a b v _ ks) = _
    rw [ih hnd hsub]
    have hstep := second_step n (x*y) a b v ks.toFinset k (by simpa using hk)
    have hcap := second_prefix_cap n x y B a b v (insert k ks.toFinset)
      (by simpa using hks) hB hx hy ha hb hv hadmit
    rw [hstep] at hcap
    have hc : |a k*(b k-v k) + secondState n (x*y) a b v ks.toFinset| ≤ coefficientCap := by
      unfold coefficientCap
      simpa [add_comm] using hcap
    unfold scaled
    rw [roundedFMA_coefficients M (a k) (b k-v k) _ hc]
    simp only [List.toFinset_cons]
    rw [hstep]
    congr 2
    ring

/-- The two genuinely rounded loops compose to the exact target dot coefficient.
No exactness of any intermediate appears among the hypotheses. -/
theorem rounded_groups_exact (s : ℝ) (M : RoundingModel (s*s))
    (n : ℕ) (x y B : ℤ) (a u b v : ℕ → ℤ)
    (hB : 0 ≤ B) (hx : |x| ≤ B) (hy : |y| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hu : ∀ k, |u k| ≤ B)
    (hb : ∀ k, |b k| ≤ B) (hv : ∀ k, |v k| ≤ B)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53)
    (left right : List ℕ) (hleft : left.Nodup) (hright : right.Nodup)
    (hls : left.toFinset = support n a u) (hrs : right.toFinset = support n b v) :
    roundedSecond s M.round a b v
      (roundedFirst s M.round a u v
        (scaled s (firstState n (x*y) a u v ∅)) left) right =
      scaled s (x*y + ∑ k ∈ Finset.range n, a k*b k) := by
  rw [rounded_first_exact s M n x y B a u v hB hx hy ha hu hv hadmit left hleft
    (by rw [hls]; exact Finset.filter_subset _ _)]
  rw [hls, groups_join n (x*y) a u b v]
  rw [rounded_second_exact s M n x y B a b v hB hx hy ha hb hv hadmit right hright
    (by rw [hrs]; exact Finset.filter_subset _ _)]
  rw [hrs, groups_complete]

/-- Integer direct dot, used to expose anchor and boundary coefficients. -/
def intDot (n : ℕ) (a b : ℕ → ℤ) : ℤ := ∑ k ∈ Finset.range n, a k*b k

/-- Head deletion is bounded by n products, before the new tail is inserted. -/
theorem deleted_head (n : ℕ) (p q : ℕ → ℤ) :
    -p 0*q 0 + intDot (n+1) p q =
      ∑ k ∈ Finset.range n, p (k+1)*q (k+1) := by
  unfold intDot
  rw [Finset.sum_range_succ']
  ring

/-- Both endpoint FMAs are exact; the initial predecessor is a direct dot value. -/
theorem rounded_boundary_exact (s : ℝ) (M : RoundingModel (s*s))
    (n : ℕ) (B : ℤ) (p q a b : ℕ → ℤ) (hB : 0 ≤ B)
    (hp : ∀ k, |p k| ≤ B) (hq : ∀ k, |q k| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hb : ∀ k, |b k| ≤ B)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53) :
    roundedFMA M.round ((a n : ℝ)*s) ((b n : ℝ)*s)
      (roundedFMA M.round (((-p 0 : ℤ) : ℝ)*s) ((q 0 : ℝ)*s)
        (scaled s (intDot (n+1) p q))) =
    scaled s (firstState n (a n*b n) a (fun k => p (k+1)) (fun k => q (k+1)) ∅) := by
  have hs := product_subset_bound (Finset.range n) (fun k => p (k+1))
    (fun k => q (k+1)) B hB (fun k => hp (k+1)) (fun k => hq (k+1))
  simp only [Finset.card_range] at hs
  have hdel : |-p 0*q 0 + intDot (n+1) p q| ≤ coefficientCap := by
    rw [deleted_head]
    unfold coefficientCap
    have hnon : 0 ≤ B^2 := sq_nonneg B
    have hlen : (n : ℤ)*B^2 ≤ ((n+1 : ℕ) : ℤ)*B^2 := by
      push_cast
      nlinarith
    exact hs.trans (hlen.trans hadmit)
  have hins := tail_sum_bound n (a n) (b n) B (fun k => p (k+1))
    (fun k => q (k+1)) hB (ha n) (hb n) (fun k => hp (k+1)) (fun k => hq (k+1))
  have hcap : |a n*b n + (-p 0*q 0 + intDot (n+1) p q)| ≤ coefficientCap := by
    rw [deleted_head]
    exact hins.trans hadmit
  unfold scaled
  rw [roundedFMA_coefficients M (-p 0) (q 0) _ hdel]
  rw [roundedFMA_coefficients M (a n) (b n) _ hcap]
  rw [deleted_head]
  simp only [firstState, Finset.sum_empty, add_zero]

/-- A concrete rounded FMA chain for an arbitrary anchor coordinate enumeration. -/
def roundedAnchor (s : ℝ) (round : ℝ → ℝ) (a b : ℕ → ℤ) (ks : List ℕ) : ℝ :=
  ks.foldr (fun k c => roundedFMA round ((a k : ℝ)*s) ((b k : ℝ)*s) c) 0

/-- All direct-anchor FMA chain states are exact from the same subset bound. -/
theorem rounded_anchor_exact (s : ℝ) (M : RoundingModel (s*s))
    (n : ℕ) (B : ℤ) (a b : ℕ → ℤ) (hB : 0 ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hb : ∀ k, |b k| ≤ B)
    (hadmit : (n : ℤ)*B^2 ≤ 2^53) (ks : List ℕ)
    (hnd : ks.Nodup) (hks : ks.toFinset ⊆ Finset.range n) :
    roundedAnchor s M.round a b ks = scaled s (∑ k ∈ ks.toFinset, a k*b k) := by
  induction ks with
  | nil => simp [roundedAnchor, scaled]
  | cons k ks ih =>
    obtain ⟨hk, hnd⟩ := List.nodup_cons.mp hnd
    have hsub : ks.toFinset ⊆ Finset.range n := by
      intro j hj
      exact hks (by simp only [List.toFinset_cons]; exact Finset.mem_insert_of_mem hj)
    change roundedFMA M.round _ _ (roundedAnchor s M.round a b ks) = _
    rw [ih hnd hsub]
    have hs := product_subset_bound (insert k ks.toFinset) a b B hB ha hb
    have hc := Finset.card_le_card hks
    simp only [List.toFinset_cons, Finset.card_range] at hc
    have hle : ((insert k ks.toFinset).card : ℤ)*B^2 ≤ (n : ℤ)*B^2 :=
      mul_le_mul_of_nonneg_right (by exact_mod_cast hc) (sq_nonneg B)
    rw [Finset.sum_insert (by simpa using hk)] at hs
    unfold scaled
    rw [roundedFMA_coefficients M (a k) (b k) _ (hs.trans (hle.trans hadmit))]
    simp only [List.toFinset_cons, Finset.sum_insert (by simpa using hk : k ∉ ks.toFinset)]

/-- A reduction node merging disjoint chains is exact, whatever the tree shape. -/
theorem rounded_anchor_merge (s : ℝ) (M : RoundingModel (s*s))
    (n : ℕ) (B : ℤ) (a b : ℕ → ℤ) (hB : 0 ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hb : ∀ k, |b k| ≤ B)
    (hadmit : (n : ℤ)*B^2 ≤ 2^53) (P Q : Finset ℕ)
    (hPQ : Disjoint P Q) (hsub : P ∪ Q ⊆ Finset.range n) :
    roundedAdd M.round (scaled s (∑ k ∈ P, a k*b k))
      (scaled s (∑ k ∈ Q, a k*b k)) =
      scaled s (∑ k ∈ P ∪ Q, a k*b k) := by
  have hs := product_subset_bound (P ∪ Q) a b B hB ha hb
  have hc := Finset.card_le_card hsub
  simp only [Finset.card_range] at hc
  have hle : ((P ∪ Q).card : ℤ)*B^2 ≤ (n : ℤ)*B^2 :=
    mul_le_mul_of_nonneg_right (by exact_mod_cast hc) (sq_nonneg B)
  rw [Finset.sum_union hPQ] at hs ⊢
  unfold scaled
  exact roundedAdd_coefficients M _ _ (hs.trans (hle.trans hadmit))

/-- Admission implies the residual coefficient fits well inside the exact cap. -/
theorem residual_cap (n : ℕ) (B a u : ℤ) (hB : 0 ≤ B)
    (ha : |a| ≤ B) (hu : |u| ≤ B)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53) : |a-u| ≤ coefficientCap := by
  have hlen : (1 : ℤ) ≤ ((n+1 : ℕ) : ℤ) := by omega
  have hbb : B^2 ≤ (2 : ℤ)^53 := by
    have h := mul_le_mul_of_nonneg_right hlen (sq_nonneg B)
    simp only [one_mul] at h
    exact h.trans hadmit
  have hsmall : 2*B ≤ coefficientCap := by
    unfold coefficientCap
    norm_num at hbb ⊢
    nlinarith
  exact (residual_bound a u B ha hu).trans hsmall

/-- Actual rounded residual subtraction is exact under admission. -/
theorem rounded_residual_exact (s : ℝ) (M : RoundingModel s)
    (n : ℕ) (B a u : ℤ) (hB : 0 ≤ B) (ha : |a| ≤ B) (hu : |u| ≤ B)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53) :
    roundedSub M.round ((a : ℝ)*s) ((u : ℝ)*s) = ((a-u : ℤ) : ℝ)*s :=
  roundedSub_coefficients M a u (residual_cap n B a u hB ha hu hadmit)

/-- The arithmetic recurrence, including both rounded endpoint FMAs. Residual
operands use exact coefficient differences, justified by rounded_residual_exact. -/
def roundedUpdate (s : ℝ) (round : ℝ → ℝ) (n : ℕ)
    (p q a b : ℕ → ℤ) (g : ℝ) (left right : List ℕ) : ℝ :=
  roundedSecond s round a b (fun k => q (k+1))
    (roundedFirst s round a (fun k => p (k+1)) (fun k => q (k+1))
      (roundedFMA round ((a n : ℝ)*s) ((b n : ℝ)*s)
        (roundedFMA round (((-p 0 : ℤ) : ℝ)*s) ((q 0 : ℝ)*s) g)) left) right

/-- Full rounded recurrence from an exact predecessor, all operations included. -/
theorem rounded_update_exact (s : ℝ) (M : RoundingModel (s*s))
    (n : ℕ) (B : ℤ) (p q a b : ℕ → ℤ) (hB : 0 ≤ B)
    (hp : ∀ k, |p k| ≤ B) (hq : ∀ k, |q k| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hb : ∀ k, |b k| ≤ B)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53)
    (left right : List ℕ) (hleft : left.Nodup) (hright : right.Nodup)
    (hls : left.toFinset = support n a (fun k => p (k+1)))
    (hrs : right.toFinset = support n b (fun k => q (k+1))) :
    roundedUpdate s M.round n p q a b (scaled s (intDot (n+1) p q)) left right =
      scaled s (intDot (n+1) a b) := by
  unfold roundedUpdate
  rw [rounded_boundary_exact s M n B p q a b hB hp hq ha hb hadmit]
  rw [rounded_groups_exact s M n (a n) (b n) B a (fun k => p (k+1)) b
    (fun k => q (k+1)) hB (ha n) (hb n) ha (fun k => hp (k+1)) hb
    (fun k => hq (k+1)) hadmit left right hleft hright hls hrs]
  unfold intDot
  rw [Finset.sum_range_succ]
  congr 1
  ring

/-- Canonical support enumeration (execution follows its reverse). -/
def supportList (n : ℕ) (a u : ℕ → ℤ) : List ℕ := (support n a u).toList

/-- A fully rounded diagonal reconstruction, parametrized only by actual anchors. -/
def roundedReconstruct (s : ℝ) (round : ℝ → ℝ) (n : ℕ)
    (A : ℕ → ℕ → ℤ) (anchor : ℕ → ℝ) : ℕ → ℕ → ℝ
  | i, 0 => anchor i
  | 0, j+1 => scaled s (intDot (n+1) (A 0) (A (j+1)))
  | i+1, j+1 => roundedUpdate s round n (A i) (A j) (A (i+1)) (A (j+1))
      (roundedReconstruct s round n A anchor i j)
      (supportList n (A (i+1)) (fun k => A i (k+1)))
      (supportList n (A (j+1)) (fun k => A j (k+1)))

/-- Universal rounded lower-triangle reconstruction. The sole arithmetic-model
hypothesis fixes concrete bounded integer-grid numbers, not Gram conclusions. -/
theorem rounded_reconstruction (s : ℝ) (M : RoundingModel (s*s))
    (n : ℕ) (B : ℤ) (A : ℕ → ℕ → ℤ) (hB : 0 ≤ B)
    (hA : ∀ i k, |A i k| ≤ B) (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53)
    (anchor : ℕ → ℝ)
    (hanchor : ∀ i, anchor i = scaled s (intDot (n+1) (A i) (A 0))) :
    ∀ i j, j ≤ i → roundedReconstruct s M.round n A anchor i j =
      scaled s (intDot (n+1) (A i) (A j)) := by
  intro i
  induction i with
  | zero =>
    intro j hj
    have hz : j = 0 := by omega
    subst j
    exact hanchor 0
  | succ i ih =>
    intro j hj
    cases j with
    | zero => exact hanchor (i+1)
    | succ j =>
      simp only [roundedReconstruct]
      rw [ih j (by omega)]
      apply rounded_update_exact s M n B (A i) (A j) (A (i+1)) (A (j+1))
        hB (hA i) (hA j) (hA (i+1)) (hA (j+1)) hadmit
      all_goals first | exact Finset.nodup_toList _ | simp [supportList]

/-- Constructed direct anchors satisfy the reconstruction premise. -/
theorem computed_anchor_exact (s : ℝ) (M : RoundingModel (s*s))
    (n : ℕ) (B : ℤ) (A : ℕ → ℕ → ℤ) (hB : 0 ≤ B)
    (hA : ∀ i k, |A i k| ≤ B) (hadmit : (n : ℤ)*B^2 ≤ 2^53) (i : ℕ) :
    roundedAnchor s M.round (A i) (A 0) (Finset.range n).toList =
      scaled s (intDot n (A i) (A 0)) := by
  have h := rounded_anchor_exact s M n B (A i) (A 0) hB (hA i) (hA 0) hadmit
    (Finset.range n).toList (Finset.nodup_toList _) (by simp)
  simpa [intDot] using h

/-- Arbitrary feature dimension, with the zero-column case explicit. -/
def roundedGram (s : ℝ) (round : ℝ → ℝ) (N : ℕ) (A : ℕ → ℕ → ℤ)
    (i j : ℕ) : ℝ :=
  match N with
  | 0 => 0
  | n+1 => roundedReconstruct s round n A
      (fun r => roundedAnchor s round (A r) (A 0) (Finset.range (n+1)).toList) i j

/-- Complete arbitrary-dimension abstract arithmetic theorem, including computed
anchors and every recurrence row. Model fixes-grid is the only rounding premise. -/
theorem rounded_gram_exact (s : ℝ) (N : ℕ) (B : ℤ) (A : ℕ → ℕ → ℤ)
    (M : RoundingModel (s*s)) (hB : 0 ≤ B) (hA : ∀ i k, |A i k| ≤ B)
    (hadmit : (N : ℤ)*B^2 ≤ 2^53) (i j : ℕ) (hji : j ≤ i) :
    roundedGram s M.round N A i j = scaled s (intDot N (A i) (A j)) := by
  cases N with
  | zero => simp [roundedGram, intDot, scaled]
  | succ n =>
    exact rounded_reconstruction s M n B A hB hA hadmit _
      (computed_anchor_exact s M (n+1) B A hB hA hadmit) i j hji

/-- The integer target is the real Gram product of the scaled input rows. -/
theorem scaled_dot (s : ℝ) (N : ℕ) (a b : ℕ → ℤ) :
    scaled s (intDot N a b) = dot N (fun k => (a k : ℝ)*s) (fun k => (b k : ℝ)*s) := by
  unfold scaled intDot dot
  push_cast
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro k hk
  ring

/-- Conditional binary64 interpretation of the complete abstract Gram theorem.
The `Binary64Bridge` records exactly the remaining machine correspondence. -/
theorem binary64_gram_exact (I : Binary64Bridge) (e : ℤ)
    (hlo : -400 ≤ e) (hhi : e ≤ 200) (N : ℕ) (B : ℤ)
    (A : ℕ → ℕ → ℤ) (hB : 0 ≤ B) (hA : ∀ i k, |A i k| ≤ B)
    (hadmit : (N : ℤ)*B^2 ≤ 2^53) (i j : ℕ) (hji : j ≤ i) :
    roundedGram (dyadicScale e) I.round N A i j =
      dot N (fun k => (A i k : ℝ)*dyadicScale e)
        (fun k => (A j k : ℝ)*dyadicScale e) := by
  exact (rounded_gram_exact (dyadicScale e) N B A (I.atAccumulator e hlo hhi)
    hB hA hadmit i j hji).trans (scaled_dot _ _ _ _)

end
end GramExact
