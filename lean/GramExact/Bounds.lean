import GramExact.Core

namespace GramExact
open scoped BigOperators

/-- Every mixed factor is one of the two bounded input coefficients. -/
theorem mixed_bound (P : Finset ℕ) (a u : ℕ → ℤ) (B : ℤ)
    (ha : ∀ k, |a k| ≤ B) (hu : ∀ k, |u k| ≤ B) (k : ℕ) :
    |mixed P a u k| ≤ B := by
  unfold mixed
  split_ifs <;> [exact ha k; exact hu k]

/-- Any subset of input products; covers anchor chains and disjoint reductions. -/
theorem product_subset_bound (s : Finset ℕ) (a b : ℕ → ℤ) (B : ℤ)
    (hB : 0 ≤ B) (ha : ∀ k, |a k| ≤ B) (hb : ∀ k, |b k| ≤ B) :
    |∑ k ∈ s, a k * b k| ≤ (s.card : ℤ) * B^2 := by
  calc
    |∑ k ∈ s, a k * b k| ≤ ∑ k ∈ s, |a k * b k| :=
      Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _k ∈ s, B^2 := by
      apply Finset.sum_le_sum
      intro k hk
      rw [abs_mul, pow_two]
      exact mul_le_mul (ha k) (hb k) (abs_nonneg _) hB
    _ = (s.card : ℤ) * B^2 := by simp

/-- A tail product plus n mixed products has coefficient at most (n+1)B². -/
theorem tail_sum_bound (n : ℕ) (x y B : ℤ) (a b : ℕ → ℤ)
    (hB : 0 ≤ B) (hx : |x| ≤ B) (hy : |y| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hb : ∀ k, |b k| ≤ B) :
    |x*y + ∑ k ∈ Finset.range n, a k*b k| ≤ ((n+1 : ℕ) : ℤ)*B^2 := by
  have ht : |x*y| ≤ B^2 := by
    rw [abs_mul, pow_two]
    exact mul_le_mul hx hy (abs_nonneg _) hB
  have hs := product_subset_bound (Finset.range n) a b B hB ha hb
  simp only [Finset.card_range] at hs
  calc
    _ ≤ |x*y| + |∑ k ∈ Finset.range n, a k*b k| := abs_add_le _ _
    _ ≤ B^2 + (n : ℤ)*B^2 := add_le_add ht hs
    _ = _ := by push_cast; ring

/-- Complete first-loop intermediate bound, with no support-size hypothesis. -/
theorem first_prefix_bound (n : ℕ) (x y B : ℤ) (a u v : ℕ → ℤ)
    (P : Finset ℕ) (hP : P ⊆ Finset.range n) (hB : 0 ≤ B)
    (hx : |x| ≤ B) (hy : |y| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hu : ∀ k, |u k| ≤ B) (hv : ∀ k, |v k| ≤ B) :
    |firstState n (x*y) a u v P| ≤ ((n+1 : ℕ) : ℤ)*B^2 := by
  rw [first_prefix n (x*y) a u v P hP]
  exact tail_sum_bound n x y B (mixed P a u) v hB hx hy
    (mixed_bound P a u B ha hu) hv

/-- Complete second-loop intermediate bound, including overlapping supports. -/
theorem second_prefix_bound (n : ℕ) (x y B : ℤ) (a b v : ℕ → ℤ)
    (Q : Finset ℕ) (hQ : Q ⊆ Finset.range n) (hB : 0 ≤ B)
    (hx : |x| ≤ B) (hy : |y| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hb : ∀ k, |b k| ≤ B) (hv : ∀ k, |v k| ≤ B) :
    |secondState n (x*y) a b v Q| ≤ ((n+1 : ℕ) : ℤ)*B^2 := by
  rw [second_prefix n (x*y) a b v Q hQ]
  exact tail_sum_bound n x y B a (mixed Q b v) hB hx hy ha
    (mixed_bound Q b v B hb hv)

/-- The admission inequality implies the binary64 coefficient cap at every prefix. -/
theorem prefix_cap (n : ℕ) (B z : ℤ)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53)
    (hbound : |z| ≤ ((n+1 : ℕ) : ℤ)*B^2) : |z| ≤ 2^53 :=
  le_trans hbound hadmit

/-- Exact residual subtraction is governed by a much smaller coefficient bound. -/
theorem residual_bound (a u B : ℤ) (ha : |a| ≤ B) (hu : |u| ≤ B) :
    |a-u| ≤ 2*B := by
  calc
    |a-u| ≤ |a| + |u| := abs_sub _ _
    _ ≤ B+B := add_le_add ha hu
    _ = 2*B := by ring

/-- Every Gram entry uses two endpoint FMAs and one per listed support coordinate. -/
def entryWork (left right : Finset ℕ) : ℕ := 2 + left.card + right.card

theorem entry_work (left right : Finset ℕ) :
    entryWork left right = 2 + left.card + right.card := rfl

/-- The full operation count before algebraic simplification of triangular sums. -/
def totalWork (m n : ℕ) (q : ℕ → ℕ) : ℕ :=
  m*n + ∑ i ∈ Finset.range m, ∑ j ∈ Finset.range i, (2 + q i + q (j+1))

/-- Each residual row occurs m times across the triangular recurrence domain. -/
theorem support_work (m : ℕ) (q : ℕ → ℕ) (hzero : q 0 = 0) :
    (∑ i ∈ Finset.range m, ∑ j ∈ Finset.range i, (q i + q (j+1))) =
      m * ∑ i ∈ Finset.range m, q i := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [Finset.sum_range_succ, ih]
    rw [Finset.sum_add_distrib]
    have hshift : (∑ j ∈ Finset.range m, q (j+1)) =
        (∑ j ∈ Finset.range m, q j) + q m := by
      have h := Finset.sum_range_succ' q m
      rw [hzero, add_zero, Finset.sum_range_succ] at h
      exact h.symm
    rw [hshift, Finset.sum_range_succ]
    simp only [Finset.sum_const, Finset.card_range, smul_eq_mul]
    ring

/-- Full Gram FMA count, with the triangular cardinality left as an explicit sum.
    This avoids natural-division conventions and is equivalent to mn+2D+mQ. -/
theorem total_work (m n : ℕ) (q : ℕ → ℕ) (hzero : q 0 = 0) :
    totalWork m n q = m*n + 2*(∑ i ∈ Finset.range m, i) +
      m*(∑ i ∈ Finset.range m, q i) := by
  unfold totalWork
  have hsplit :
      (∑ i ∈ Finset.range m, ∑ j ∈ Finset.range i, (2 + q i + q (j+1))) =
      2*(∑ i ∈ Finset.range m, i) +
        ∑ i ∈ Finset.range m, ∑ j ∈ Finset.range i, (q i + q (j+1)) := by
    simp_rw [Nat.add_assoc, Finset.sum_add_distrib]
    simp only [Finset.sum_const, Finset.card_range, smul_eq_mul]
    rw [← Finset.sum_mul]
    ring
  rw [hsplit, support_work m q hzero]
  omega

/-- Each finite first-group prefix contains exactly its old and new input factors;
    the cap is independent of the enumeration order. -/
theorem first_prefix_cap (n : ℕ) (x y B : ℤ) (a u v : ℕ → ℤ)
    (P : Finset ℕ) (hP : P ⊆ Finset.range n) (hB : 0 ≤ B)
    (hx : |x| ≤ B) (hy : |y| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hu : ∀ k, |u k| ≤ B) (hv : ∀ k, |v k| ≤ B)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53) :
    |firstState n (x*y) a u v P| ≤ 2^53 :=
  prefix_cap n B _ hadmit (first_prefix_bound n x y B a u v P hP hB hx hy ha hu hv)

theorem second_prefix_cap (n : ℕ) (x y B : ℤ) (a b v : ℕ → ℤ)
    (Q : Finset ℕ) (hQ : Q ⊆ Finset.range n) (hB : 0 ≤ B)
    (hx : |x| ≤ B) (hy : |y| ≤ B)
    (ha : ∀ k, |a k| ≤ B) (hb : ∀ k, |b k| ≤ B) (hv : ∀ k, |v k| ≤ B)
    (hadmit : ((n+1 : ℕ) : ℤ)*B^2 ≤ 2^53) :
    |secondState n (x*y) a b v Q| ≤ 2^53 :=
  prefix_cap n B _ hadmit (second_prefix_bound n x y B a b v Q hQ hB hx hy ha hb hv)
end GramExact
