import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Lean.Elab.Tactic.Omega

/-!
# Exact operations on bounded integer grids

This is a real-arithmetic rounding abstraction, not an implementation or
verification of IEEE 754 bit encodings. `RoundingModel.fixes` is the explicit
implementation assumption. The kernel derives primitive exactness from integer
coefficient bounds and polynomial identities; no Gram recurrence is assumed.
`Binary64Bridge` separately packages the external IEEE interpretation, restricted
to a conservative finite normal dyadic subset (zero included). This file does
not certify a compiler, FMA instruction, rounding mode, or reassociation policy.
-/
namespace GramExact
noncomputable section

/-- Inclusive exact-integer coefficient cap. -/
def coefficientCap : ℤ := 2 ^ 53

/-- Concrete bounded integer coefficients at a common real scale. -/
def GridRepresentable (scale x : ℝ) : Prop :=
  ∃ z : ℤ, |z| ≤ coefficientCap ∧ x = (z : ℝ) * scale

/-- Primitive specification only, with no algorithmic conclusion. -/
structure RoundingModel (scale : ℝ) where
  round : ℝ → ℝ
  fixes : ∀ x, GridRepresentable scale x → round x = x

def roundedFMA (round : ℝ → ℝ) (a b c : ℝ) : ℝ := round (a * b + c)
def roundedSub (round : ℝ → ℝ) (a b : ℝ) : ℝ := round (a - b)
def roundedAdd (round : ℝ → ℝ) (a b : ℝ) : ℝ := round (a + b)
def roundedMul (round : ℝ → ℝ) (a b : ℝ) : ℝ := round (a * b)

theorem grid_of_coefficient {s : ℝ} {z : ℤ} (hz : |z| ≤ coefficientCap) :
    GridRepresentable s ((z : ℝ) * s) := ⟨z, hz, rfl⟩

theorem round_coefficient {s : ℝ} (M : RoundingModel s) {z : ℤ}
    (hz : |z| ≤ coefficientCap) :
    M.round ((z : ℝ) * s) = (z : ℝ) * s := M.fixes _ (grid_of_coefficient hz)

/-- Only the exact FMA result needs to be representable, not its product. -/
theorem roundedFMA_exact {s a b c : ℝ} (M : RoundingModel s)
    {z : ℤ} (hz : |z| ≤ coefficientCap)
    (hresult : a * b + c = (z : ℝ) * s) :
    roundedFMA M.round a b c = a * b + c := by
  exact M.fixes _ ⟨z, hz, hresult⟩

theorem roundedSub_exact {s a b : ℝ} (M : RoundingModel s)
    {z : ℤ} (hz : |z| ≤ coefficientCap)
    (hresult : a - b = (z : ℝ) * s) :
    roundedSub M.round a b = a - b := by
  exact M.fixes _ ⟨z, hz, hresult⟩

theorem roundedAdd_exact {s a b : ℝ} (M : RoundingModel s)
    {z : ℤ} (hz : |z| ≤ coefficientCap)
    (hresult : a + b = (z : ℝ) * s) :
    roundedAdd M.round a b = a + b := by
  exact M.fixes _ ⟨z, hz, hresult⟩

/-- Subtraction closure is derived from its integer coefficient bound. -/
theorem roundedSub_coefficients {s : ℝ} (M : RoundingModel s)
    (a b : ℤ) (hbound : |a - b| ≤ coefficientCap) :
    roundedSub M.round ((a : ℝ) * s) ((b : ℝ) * s) =
      ((a - b : ℤ) : ℝ) * s := by
  have heq : (a : ℝ) * s - (b : ℝ) * s = ((a - b : ℤ) : ℝ) * s := by
    push_cast
    ring
  unfold roundedSub
  rw [heq]
  exact round_coefficient M hbound

/-- Addition closure applies to disjoint-subset anchor reduction nodes. -/
theorem roundedAdd_coefficients {s : ℝ} (M : RoundingModel s)
    (a b : ℤ) (hbound : |a + b| ≤ coefficientCap) :
    roundedAdd M.round ((a : ℝ) * s) ((b : ℝ) * s) =
      ((a + b : ℤ) : ℝ) * s := by
  have heq : (a : ℝ) * s + (b : ℝ) * s = ((a + b : ℤ) : ℝ) * s := by
    push_cast
    ring
  unfold roundedAdd
  rw [heq]
  exact round_coefficient M hbound

/-- Input products lie on the squared-scale grid. -/
theorem roundedMul_coefficients {s : ℝ} (M : RoundingModel (s * s))
    (a b : ℤ) (hbound : |a * b| ≤ coefficientCap) :
    roundedMul M.round ((a : ℝ) * s) ((b : ℝ) * s) =
      ((a * b : ℤ) : ℝ) * (s * s) := by
  have heq : ((a : ℝ) * s) * ((b : ℝ) * s) =
      ((a * b : ℤ) : ℝ) * (s * s) := by
    push_cast
    ring
  unfold roundedMul
  rw [heq]
  exact round_coefficient M hbound

/-- FMA closure is derived from an integer coefficient bound and ring identity. -/
theorem roundedFMA_coefficients {s : ℝ} (M : RoundingModel (s * s))
    (a b c : ℤ) (hbound : |a * b + c| ≤ coefficientCap) :
    roundedFMA M.round ((a : ℝ) * s) ((b : ℝ) * s)
      ((c : ℝ) * (s * s)) = ((a * b + c : ℤ) : ℝ) * (s * s) := by
  have heq : ((a : ℝ) * s) * ((b : ℝ) * s) + (c : ℝ) * (s * s) =
      ((a * b + c : ℤ) : ℝ) * (s * s) := by
    push_cast
    ring
  unfold roundedFMA
  rw [heq]
  exact round_coefficient M hbound

/-- Nonvacuity witness only; identity is not an IEEE rounding implementation. -/
def identityRoundingModel (s : ℝ) : RoundingModel s where
  round := id
  fixes := by intro x _; rfl

def dyadicScale (e : ℤ) : ℝ := (2 : ℝ) ^ e

/-- Mathematical subset description, not a machine-bit representation predicate.
The exponent range leaves room for every coefficient through `2^53`. -/
def Binary64SafeDyadic (x : ℝ) : Prop :=
  ∃ e : ℤ, -1022 ≤ e ∧ e ≤ 970 ∧ GridRepresentable (dyadicScale e) x

/-- EXTERNAL IEEE BINDING ASSUMPTION. A client must justify that its interpreted
machine rounding fixes this subset; that correspondence is not proved here. -/
structure Binary64Bridge where
  round : ℝ → ℝ
  fixes_safe : ∀ x, Binary64SafeDyadic x → round x = x

def Binary64Bridge.atExponent (I : Binary64Bridge) (e : ℤ)
    (hlo : -1022 ≤ e) (hhi : e ≤ 970) : RoundingModel (dyadicScale e) where
  round := I.round
  fixes := by
    intro x hx
    exact I.fixes_safe x ⟨e, hlo, hhi, hx⟩

theorem input_exponent_safe {e : ℤ} (hlo : -400 ≤ e) (hhi : e ≤ 200) :
    -1022 ≤ e ∧ e ≤ 970 := by omega

theorem accumulator_exponent_safe {e : ℤ} (hlo : -400 ≤ e) (hhi : e ≤ 200) :
    -1022 ≤ 2 * e ∧ 2 * e ≤ 970 := by omega

/-- The squared input scale is the accumulator scale used in the report. -/
theorem dyadicScale_square (e : ℤ) :
    dyadicScale e * dyadicScale e = dyadicScale (2 * e) := by
  unfold dyadicScale
  rw [two_mul]
  exact (zpow_add₀ (by norm_num : (2 : ℝ) ≠ 0) e e).symm

/-- A single external rounding function handles the input and squared grids. -/
def Binary64Bridge.atAccumulator (I : Binary64Bridge) (e : ℤ)
    (hlo : -400 ≤ e) (hhi : e ≤ 200) :
    RoundingModel (dyadicScale e * dyadicScale e) where
  round := I.round
  fixes := by
    intro x hx
    rw [dyadicScale_square] at hx
    obtain ⟨hl, hu⟩ := accumulator_exponent_safe hlo hhi
    exact I.fixes_safe x ⟨2 * e, hl, hu, hx⟩

/-- Conditional IEEE bridge theorem: primitive exactness is kernel-derived from
an integer coefficient bound; the sole machine-binding assumption is `I`. -/
theorem binary64_fma_coefficients (I : Binary64Bridge) (e : ℤ)
    (hlo : -400 ≤ e) (hhi : e ≤ 200) (a b c : ℤ)
    (hbound : |a * b + c| ≤ coefficientCap) :
    roundedFMA I.round ((a : ℝ) * dyadicScale e) ((b : ℝ) * dyadicScale e)
      ((c : ℝ) * (dyadicScale e * dyadicScale e)) =
      ((a * b + c : ℤ) : ℝ) * (dyadicScale e * dyadicScale e) :=
  roundedFMA_coefficients (I.atAccumulator e hlo hhi) a b c hbound

theorem binary64_sub_coefficients (I : Binary64Bridge) (e : ℤ)
    (hlo : -400 ≤ e) (hhi : e ≤ 200) (a b : ℤ)
    (hbound : |a - b| ≤ coefficientCap) :
    roundedSub I.round ((a : ℝ) * dyadicScale e) ((b : ℝ) * dyadicScale e) =
      ((a - b : ℤ) : ℝ) * dyadicScale e := by
  obtain ⟨hl, hu⟩ := input_exponent_safe hlo hhi
  exact roundedSub_coefficients (I.atExponent e hl hu) a b hbound

set_option pp.proofs false
#print roundedFMA_coefficients
#print roundedSub_coefficients
#print Binary64Bridge.atExponent
#print axioms roundedFMA_coefficients
#print axioms roundedSub_coefficients
#print axioms Binary64Bridge.atExponent
#print binary64_fma_coefficients
#print axioms binary64_fma_coefficients
#print axioms binary64_sub_coefficients
end
end GramExact
