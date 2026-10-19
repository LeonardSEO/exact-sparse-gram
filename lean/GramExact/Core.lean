import Mathlib.Data.Real.Basic
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Abel
import Lean.Elab.Tactic.Omega

/-! Exact real algebra for the two CSR loops. A prefix is a finite set of
already visited distinct coordinates; this permits every enumeration order. -/
noncomputable section
namespace GramExact
open scoped BigOperators

def mixed {R : Type*} (P : Finset ℕ) (a u : ℕ → R) (k : ℕ) : R :=
  if k ∈ P then a k else u k

def firstState {R : Type*} [CommRing R] (n : ℕ) (t : R)
    (a u v : ℕ → R) (P : Finset ℕ) : R :=
  t + ∑ k ∈ Finset.range n, u k * v k + ∑ k ∈ P, (a k - u k) * v k

def secondState {R : Type*} [CommRing R] (n : ℕ) (t : R)
    (a b v : ℕ → R) (Q : Finset ℕ) : R :=
  t + ∑ k ∈ Finset.range n, a k * v k + ∑ k ∈ Q, a k * (b k - v k)

theorem prefix_sum {R : Type*} [CommRing R] (n : ℕ) (P : Finset ℕ)
    (hP : P ⊆ Finset.range n) (f : ℕ → R) :
    (∑ k ∈ Finset.range n, if k ∈ P then f k else 0) = ∑ k ∈ P, f k := by
  rw [Finset.sum_ite_mem, Finset.inter_eq_right.mpr hP]

/-- The complete first CSR-group prefix invariant, over any commutative ring. -/
theorem first_prefix {R : Type*} [CommRing R] (n : ℕ) (t : R)
    (a u v : ℕ → R) (P : Finset ℕ) (hP : P ⊆ Finset.range n) :
    firstState n t a u v P = t + ∑ k ∈ Finset.range n, mixed P a u k * v k := by
  unfold firstState
  rw [← prefix_sum n P hP (fun k => (a k - u k) * v k)]
  rw [add_assoc, ← Finset.sum_add_distrib]
  congr 1
  apply Finset.sum_congr rfl
  intro k hk
  unfold mixed
  split_ifs <;> ring

/-- The complete second CSR-group prefix invariant. -/
theorem second_prefix {R : Type*} [CommRing R] (n : ℕ) (t : R)
    (a b v : ℕ → R) (Q : Finset ℕ) (hQ : Q ⊆ Finset.range n) :
    secondState n t a b v Q = t + ∑ k ∈ Finset.range n, a k * mixed Q b v k := by
  unfold secondState
  rw [← prefix_sum n Q hQ (fun k => a k * (b k - v k))]
  rw [add_assoc, ← Finset.sum_add_distrib]
  congr 1
  apply Finset.sum_congr rfl
  intro k hk
  unfold mixed
  split_ifs <;> ring

theorem first_step {R : Type*} [CommRing R] (n : ℕ) (t : R)
    (a u v : ℕ → R) (P : Finset ℕ) (k : ℕ) (hk : k ∉ P) :
    firstState n t a u v (insert k P) =
      firstState n t a u v P + (a k - u k) * v k := by
  simp only [firstState, Finset.sum_insert hk]
  ring

theorem second_step {R : Type*} [CommRing R] (n : ℕ) (t : R)
    (a b v : ℕ → R) (Q : Finset ℕ) (k : ℕ) (hk : k ∉ Q) :
    secondState n t a b v (insert k Q) =
      secondState n t a b v Q + a k * (b k - v k) := by
  simp only [secondState, Finset.sum_insert hk]
  ring

def support {R : Type*} [DecidableEq R] (n : ℕ) (a u : ℕ → R) : Finset ℕ :=
  (Finset.range n).filter (fun k => a k ≠ u k)

theorem mixed_support {R : Type*} [DecidableEq R] (n : ℕ) (a u : ℕ → R) (k : ℕ)
    (hk : k ∈ Finset.range n) : mixed (support n a u) a u k = a k := by
  unfold mixed support
  by_cases h : a k = u k <;> simp [hk, h]

theorem groups_join {R : Type*} [CommRing R] [DecidableEq R] (n : ℕ) (t : R) (a u b v : ℕ → R) :
    firstState n t a u v (support n a u) = secondState n t a b v ∅ := by
  rw [first_prefix n t a u v _ (by exact Finset.filter_subset _ _)]
  simp only [secondState, Finset.sum_empty, add_zero]
  congr 1
  apply Finset.sum_congr rfl
  intro k hk
  rw [mixed_support n a u k hk]

theorem groups_complete {R : Type*} [CommRing R] [DecidableEq R] (n : ℕ) (t : R) (a b v : ℕ → R) :
    secondState n t a b v (support n b v) =
      t + ∑ k ∈ Finset.range n, a k * b k := by
  rw [second_prefix n t a b v _ (by exact Finset.filter_subset _ _)]
  congr 1
  apply Finset.sum_congr rfl
  intro k hk
  rw [mixed_support n b v k hk]

/-- Direct dot product, valid also for zero columns. -/
def dot (n : ℕ) (a b : ℕ → ℝ) : ℝ := ∑ k ∈ Finset.range n, a k * b k

/-- Boundary deletion and insertion, exactly in the implemented order. -/
def boundary (n : ℕ) (p q a b : ℕ → ℝ) (g : ℝ) : ℝ :=
  (g - p 0 * q 0) + a n * b n

theorem boundary_identity (n : ℕ) (p q a b : ℕ → ℝ) :
    boundary n p q a b (dot (n+1) p q) =
      a n * b n + ∑ k ∈ Finset.range n, p (k+1) * q (k+1) := by
  unfold boundary dot
  rw [Finset.sum_range_succ']
  ring

/-- One entry, with full first group completed before the second. -/
def update (n : ℕ) (p q a b : ℕ → ℝ) (g : ℝ) : ℝ :=
  boundary n p q a b g +
    ∑ k ∈ support n a (fun k => p (k+1)), (a k - p (k+1)) * q (k+1) +
    ∑ k ∈ support n b (fun k => q (k+1)), a k * (b k - q (k+1))

theorem recurrence (n : ℕ) (p q a b : ℕ → ℝ) :
    update n p q a b (dot (n+1) p q) = dot (n+1) a b := by
  unfold update
  rw [boundary_identity]
  have h := groups_join n (a n * b n) a (fun k => p (k+1)) b (fun k => q (k+1))
  unfold firstState at h
  simp only [secondState, Finset.sum_empty, add_zero] at h
  have h' := groups_complete n (a n * b n) a b (fun k => q (k+1))
  unfold secondState at h'
  rw [h]
  rw [h']
  unfold dot
  rw [Finset.sum_range_succ]
  ring

/-- Reconstruction follows the predecessor diagonal until a first-column anchor. -/
def reconstruct (n : ℕ) (A : ℕ → ℕ → ℝ) (anchor : ℕ → ℝ) : ℕ → ℕ → ℝ
  | i, 0 => anchor i
  | 0, j+1 => dot (n+1) (A 0) (A (j+1))
  | i+1, j+1 => update n (A i) (A j) (A (i+1)) (A (j+1))
      (reconstruct n A anchor i j)

/-- Arbitrary row count, all lower entries from exact first-column anchors. -/
theorem reconstruction (n : ℕ) (A : ℕ → ℕ → ℝ) (anchor : ℕ → ℝ)
    (hanchor : ∀ i, anchor i = dot (n+1) (A i) (A 0)) :
    ∀ i j, j ≤ i → reconstruct n A anchor i j = dot (n+1) (A i) (A j) := by
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
      simp only [reconstruct]
      rw [ih j (by omega)]
      exact recurrence n (A i) (A j) (A (i+1)) (A (j+1))

end GramExact
