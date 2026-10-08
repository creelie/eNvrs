import Mathlib

/-!
# Elementary lemmas behind the blowup argument for the Katz–Pavlović dyadic model

These are the analytic and arithmetic lemmas that the computer-assisted step of the paper
"Blowup above the exponent one third in viscous Katz–Pavlović dyadic models with small
shell ratios" combines.  The interval-arithmetic verification itself is not formalized here.

* `ahead_tail`      : the induction that closes the tail of shells ahead of the window;
* `deep_tail`       : the Riccati bound for the shells behind the window;
* `supersolution`   : the algebra of the exponential supersolution in the comparison lemma;
* `section_time`    : the time shift of a transversal section crossing;
* `eps_le`          : the viscosity parameter never increases along the induction;
* `threshold_iff`   : `Λ^{2α} < r G` is equivalent to `α < 1/3 + log G / (2 log Λ)`;
* `time_sum`        : the step durations are summable, so the front arrives in finite time;
* `weighted_tendsto`: the weighted amplitudes `Λ^{(σ-1/3) n} Y_n` tend to infinity.
-/

namespace DyadicNS

open Real Filter Topology Set

/-! ## The ahead tail -/

/-- The ahead-tail induction (Lemma 4.2 of the paper).  For `j ≥ 0`, `q j` bounds the shell
`K_a + 1 + j` during one step, `θ / 2^j` is the hypothesis on that shell at the start of the step
and `R0 * r^j` its transfer rate, `qKA` bounds the last window shell. -/
theorem ahead_tail (r s θ G R0 qKA : ℝ) (q : ℕ → ℝ)
    (hr : 0 ≤ r) (hr2 : r ≤ 2) (hs : 0 ≤ s) (hθ : 0 < θ) (hR0 : 0 ≤ R0)
    (hqnn : ∀ j, 0 ≤ q j)
    (hq0 : q 0 ≤ θ + R0 * s * qKA ^ 2)
    (hq : ∀ j, q (j + 1) ≤ θ / 2 ^ (j + 1) + R0 * r ^ (j + 1) * s * q j ^ 2)
    (base : R0 * s * qKA ^ 2 ≤ θ)
    (cond1 : 4 * (R0 * r) * s * θ ≤ 1 / 2)
    (cond2 : 1 / 2 + 4 * (R0 * r) * s * θ ≤ G) :
    (∀ j, q j ≤ 2 * (θ / 2 ^ j)) ∧ (∀ j, q (j + 1) ≤ G * (θ / 2 ^ j)) := by
  -- a_j = 4 R0 r^{j+1} s θ / 2^j decreases in j, so a_j ≤ a_0
  have ha : ∀ j : ℕ, 4 * (R0 * r ^ (j + 1)) * s * θ / 2 ^ j ≤ 4 * (R0 * r) * s * θ := by
    intro j
    have h2 : (0 : ℝ) < 2 ^ j := by positivity
    have hpow : r ^ j ≤ 2 ^ j := pow_le_pow_left₀ hr hr2 j
    rw [div_le_iff₀ h2, pow_succ]
    have hc : 0 ≤ 4 * R0 * r * s * θ := by positivity
    nlinarith [mul_le_mul_of_nonneg_left hpow hc]
  -- the key one-step estimate
  have key : ∀ j, q j ≤ 2 * (θ / 2 ^ j) →
      q (j + 1) ≤ θ / 2 ^ j * (1 / 2 + 4 * (R0 * r) * s * θ) := by
    intro j hj
    have h2 : (0 : ℝ) < 2 ^ j := by positivity
    have hsq : q j ^ 2 ≤ (2 * (θ / 2 ^ j)) ^ 2 := pow_le_pow_left₀ (hqnn j) hj 2
    have hcoef : 0 ≤ R0 * r ^ (j + 1) * s := by positivity
    have h1 := hq j
    have h3 : R0 * r ^ (j + 1) * s * q j ^ 2 ≤ R0 * r ^ (j + 1) * s * (2 * (θ / 2 ^ j)) ^ 2 :=
      mul_le_mul_of_nonneg_left hsq hcoef
    have e1 : θ / 2 ^ (j + 1) = θ / 2 ^ j * (1 / 2) := by
      rw [pow_succ]; field_simp
    have e2 : R0 * r ^ (j + 1) * s * (2 * (θ / 2 ^ j)) ^ 2
        = θ / 2 ^ j * (4 * (R0 * r ^ (j + 1)) * s * θ / 2 ^ j) := by
      field_simp; ring
    have h4 : θ / 2 ^ j * (4 * (R0 * r ^ (j + 1)) * s * θ / 2 ^ j)
        ≤ θ / 2 ^ j * (4 * (R0 * r) * s * θ) :=
      mul_le_mul_of_nonneg_left (ha j) (by positivity)
    calc q (j + 1) ≤ θ / 2 ^ (j + 1) + R0 * r ^ (j + 1) * s * q j ^ 2 := h1
      _ ≤ θ / 2 ^ j * (1 / 2) + θ / 2 ^ j * (4 * (R0 * r) * s * θ) := by
          rw [e1] at *; linarith [h3, h4, e2.le, e2.ge]
      _ = θ / 2 ^ j * (1 / 2 + 4 * (R0 * r) * s * θ) := by ring
  have first : ∀ j, q j ≤ 2 * (θ / 2 ^ j) := by
    intro j
    induction j with
    | zero => simp; linarith
    | succ j ih =>
      have h := key j ih
      have hpos : 0 ≤ θ / 2 ^ j := by positivity
      have : θ / 2 ^ j * (1 / 2 + 4 * (R0 * r) * s * θ) ≤ θ / 2 ^ j * 1 :=
        mul_le_mul_of_nonneg_left (by linarith) hpos
      have e : 2 * (θ / 2 ^ (j + 1)) = θ / 2 ^ j := by rw [pow_succ]; field_simp
      rw [e]; linarith
  refine ⟨first, fun j => ?_⟩
  have h := key j (first j)
  have hpos : 0 ≤ θ / 2 ^ j := by positivity
  have : θ / 2 ^ j * (1 / 2 + 4 * (R0 * r) * s * θ) ≤ θ / 2 ^ j * G :=
    mul_le_mul_of_nonneg_left cond2 hpos
  linarith [mul_comm G (θ / 2 ^ j)]

/-! ## The deep tail -/

/-- Strict comparison with a Riccati barrier, for a perturbed barrier. -/
theorem deep_tail_aux (f f' : ℝ → ℝ) (R W b δ : ℝ) (hR : 0 ≤ R) (hW : 0 < W) (hδ : 0 < δ)
    (hb : (W + δ) * (R + δ) * b < 1)
    (hf : ContinuousOn f (Icc 0 b)) (hf' : ∀ x ∈ Ico 0 b, HasDerivWithinAt f (f' x) (Ici x) x)
    (hbound : ∀ x ∈ Ico 0 b, f' x ≤ R * f x ^ 2) (h0 : f 0 ≤ W) :
    ∀ x ∈ Icc 0 b, f x ≤ (W + δ) / (1 - (W + δ) * (R + δ) * x) := by
  set c := (W + δ) * (R + δ) with hc
  have hcpos : 0 < c := by positivity
  have hden : ∀ x ∈ Icc (0 : ℝ) b, 0 < 1 - c * x := by
    intro x hx
    have : c * x ≤ c * b := mul_le_mul_of_nonneg_left hx.2 hcpos.le
    linarith
  set B : ℝ → ℝ := fun x => (W + δ) / (1 - c * x)
  set B' : ℝ → ℝ := fun x => (W + δ) * c / (1 - c * x) ^ 2
  have hBd : ∀ x ∈ Icc (0 : ℝ) b, HasDerivAt B (B' x) x := by
    intro x hx
    have hne : 1 - c * x ≠ 0 := (hden x hx).ne'
    have h1 : HasDerivAt (fun x => 1 - c * x) (-c) x := by
      simpa using (hasDerivAt_id x).const_mul c |>.const_sub 1
    have h2 := (h1.inv hne).const_mul (W + δ)
    convert h2 using 1
    · funext y; simp [B, div_eq_mul_inv]
    · simp only [B']; field_simp
  have hB0 : f 0 ≤ B 0 := by simp only [B]; simp; linarith
  have hBc : ContinuousOn B (Icc 0 b) := fun x hx => (hBd x hx).continuousAt.continuousWithinAt
  have hBd' : ∀ x ∈ Ico (0 : ℝ) b, HasDerivWithinAt B (B' x) (Ici x) x :=
    fun x hx => (hBd x (Ico_subset_Icc_self hx)).hasDerivWithinAt
  have bound : ∀ x ∈ Ico (0 : ℝ) b, f x = B x → f' x < B' x := by
    intro x hx hfx
    have hxI : x ∈ Icc (0 : ℝ) b := Ico_subset_Icc_self hx
    have hd := hden x hxI
    have hBpos : 0 < B x := by simp only [B]; positivity
    have hB' : B' x = (R + δ) * B x ^ 2 := by
      simp only [B', B, hc]; field_simp
    calc f' x ≤ R * f x ^ 2 := hbound x hx
      _ = R * B x ^ 2 := by rw [hfx]
      _ < (R + δ) * B x ^ 2 := by nlinarith [pow_pos hBpos 2]
      _ = B' x := hB'.symm
  intro x hx
  exact image_le_of_deriv_right_lt_deriv_boundary' hf hf' hB0 hBc hBd' bound hx

/-- The deep-tail bound (Lemma 4.1 of the paper): if `f' ≤ R f^2` and `f 0 ≤ W` then
`f x ≤ W / (1 - W R x)` as long as `W R x < 1`. -/
theorem deep_tail (f f' : ℝ → ℝ) (R W b : ℝ) (hR : 0 ≤ R) (hW : 0 < W)
    (hb : W * R * b < 1)
    (hf : ContinuousOn f (Icc 0 b)) (hf' : ∀ x ∈ Ico 0 b, HasDerivWithinAt f (f' x) (Ici x) x)
    (hbound : ∀ x ∈ Ico 0 b, f' x ≤ R * f x ^ 2) (h0 : f 0 ≤ W) :
    ∀ x ∈ Icc 0 b, f x ≤ W / (1 - W * R * x) := by
  intro x hx
  -- the perturbed barriers converge to the unperturbed one as δ → 0
  have hcont : ContinuousWithinAt
      (fun δ : ℝ => (W + δ) / (1 - (W + δ) * (R + δ) * x)) (Ioi 0) 0 := by
    have hx0 : W * R * x < 1 := by
      have : W * R * x ≤ W * R * b := mul_le_mul_of_nonneg_left hx.2 (by positivity)
      linarith
    apply ContinuousAt.continuousWithinAt
    apply ContinuousAt.div (by fun_prop) (by fun_prop)
    simp; linarith
  -- for small δ the hypothesis of `deep_tail_aux` holds
  have hsmall : ∀ᶠ δ in 𝓝[>] (0 : ℝ), (W + δ) * (R + δ) * b < 1 := by
    have : ContinuousAt (fun δ : ℝ => (W + δ) * (R + δ) * b) 0 := by fun_prop
    have hv : (W + 0) * (R + 0) * b < 1 := by simpa using hb
    have h' : ∀ᶠ δ in 𝓝 (0 : ℝ), (W + δ) * (R + δ) * b < 1 :=
      this.eventually_lt continuousAt_const hv
    exact nhdsWithin_le_nhds h'
  have hle : ∀ᶠ δ in 𝓝[>] (0 : ℝ), f x ≤ (W + δ) / (1 - (W + δ) * (R + δ) * x) := by
    filter_upwards [hsmall, self_mem_nhdsWithin] with δ hδb hδ
    exact deep_tail_aux f f' R W b δ hR hW hδ hδb hf hf' hbound h0 x hx
  have hlim := hcont.tendsto
  simp only [add_zero] at hlim
  exact ge_of_tendsto hlim hle

/-! ## The comparison lemma -/

/-- The exponential supersolution `w(s) = η (e^{λ s} - 1)/λ ζ` of `w' = M w + e`
(Lemma 4.3 of the paper): if `M ζ ≤ λ ζ` and `e ≤ η ζ`, then `M w(s) + e ≤ w'(s)`. -/
theorem supersolution {ι : Type*} [Fintype ι] (M : ι → ι → ℝ) (ζ e : ι → ℝ) (lam η s : ℝ)
    (hlam : 0 < lam) (hη : 0 ≤ η) (hs : 0 ≤ s)
    (hMζ : ∀ i, ∑ j, M i j * ζ j ≤ lam * ζ i) (he : ∀ i, e i ≤ η * ζ i) (i : ι) :
    ∑ j, M i j * (η * (exp (lam * s) - 1) / lam * ζ j) + e i ≤ η * exp (lam * s) * ζ i := by
  set c := η * (exp (lam * s) - 1) / lam with hcdef
  have hc : 0 ≤ c := by
    have : 0 ≤ exp (lam * s) - 1 := by
      have := Real.add_one_le_exp (lam * s); nlinarith [mul_nonneg hlam.le hs]
    positivity
  have hsum : ∑ j, M i j * (c * ζ j) = c * ∑ j, M i j * ζ j := by
    rw [Finset.mul_sum]; congr 1; funext j; ring
  rw [hsum]
  have h1 : c * ∑ j, M i j * ζ j ≤ c * (lam * ζ i) := mul_le_mul_of_nonneg_left (hMζ i) hc
  have h2 : c * (lam * ζ i) + η * ζ i = η * exp (lam * s) * ζ i := by
    rw [hcdef]; field_simp; ring
  linarith [he i]

/-! ## The section time -/

/-- If `h' ≥ m > 0` on `[a, b]`, `h (s₂) = 0` and `|h (s₁)| ≤ ξ`, then `|s₂ - s₁| ≤ ξ / m`
(Lemma 4.4 of the paper). -/
theorem section_time (h : ℝ → ℝ) (a b m ξ s₁ s₂ : ℝ) (hm : 0 < m)
    (hcont : ContinuousOn h (Icc a b)) (hdiff : DifferentiableOn ℝ h (interior (Icc a b)))
    (hmono : ∀ x ∈ interior (Icc a b), m ≤ deriv h x)
    (hs₁ : s₁ ∈ Icc a b) (hs₂ : s₂ ∈ Icc a b) (hzero : h s₂ = 0) (hsmall : |h s₁| ≤ ξ) :
    |s₂ - s₁| ≤ ξ / m := by
  have grow := (convex_Icc a b).mul_sub_le_image_sub_of_le_deriv hcont hdiff hmono
  rw [le_div_iff₀ hm]
  have hs := abs_le.mp hsmall
  rcases le_total s₁ s₂ with h12 | h21
  · have := grow s₁ hs₁ s₂ hs₂ h12
    rw [abs_of_nonneg (by linarith)]
    nlinarith [hs.1, hs.2]
  · have := grow s₂ hs₂ s₁ hs₁ h21
    rw [abs_of_nonpos (by linarith)]
    nlinarith [hs.1, hs.2]

/-! ## The induction over steps -/

/-- The viscosity parameter contracts: `ε_{n+1} ≤ q ε_n` with `q ≤ 1` and `ε_n ≥ 0` gives
`ε_n ≤ ε_0`. -/
theorem eps_le (ε : ℕ → ℝ) (q : ℝ) (hq1 : q ≤ 1) (hε : ∀ n, 0 ≤ ε n)
    (hstep : ∀ n, ε (n + 1) ≤ q * ε n) : ∀ n, ε n ≤ ε 0 := by
  intro n
  induction n with
  | zero => exact le_refl _
  | succ n ih =>
    calc ε (n + 1) ≤ q * ε n := hstep n
      _ ≤ 1 * ε n := mul_le_mul_of_nonneg_right hq1 (hε n)
      _ = ε n := one_mul _
      _ ≤ ε 0 := ih

/-- `Λ^{2α} < Λ^{2/3} G` if and only if `α < 1/3 + log G / (2 log Λ)`. -/
theorem threshold_iff (Λ α G : ℝ) (hΛ : 1 < Λ) (hG : 0 < G) :
    Λ ^ (2 * α) < Λ ^ (2 / 3 : ℝ) * G ↔ α < 1 / 3 + Real.log G / (2 * Real.log Λ) := by
  have hΛ0 : 0 < Λ := by linarith
  have hlog : 0 < Real.log Λ := Real.log_pos hΛ
  rw [← Real.log_lt_log_iff (by positivity) (by positivity), Real.log_mul (by positivity) hG.ne',
    Real.log_rpow hΛ0, Real.log_rpow hΛ0]
  constructor
  · intro h
    rw [← sub_pos]
    have : 0 < (2 / 3 * Real.log Λ + Real.log G - 2 * α * Real.log Λ) / (2 * Real.log Λ) := by
      apply div_pos <;> linarith
    have e : (2 / 3 * Real.log Λ + Real.log G - 2 * α * Real.log Λ) / (2 * Real.log Λ)
        = 1 / 3 + Real.log G / (2 * Real.log Λ) - α := by
      field_simp
    linarith
  · intro h
    have e : 2 * α * Real.log Λ < 2 * Real.log Λ * (1 / 3 + Real.log G / (2 * Real.log Λ)) := by
      have := mul_lt_mul_of_pos_left h (by positivity : 0 < 2 * Real.log Λ)
      linarith
    have e2 : 2 * Real.log Λ * (1 / 3 + Real.log G / (2 * Real.log Λ))
        = 2 / 3 * Real.log Λ + Real.log G := by field_simp
    linarith

/-- The durations of the steps are summable when they decay geometrically. -/
theorem time_sum (Δ : ℕ → ℝ) (a q : ℝ) (hq0 : 0 ≤ q) (hq : q < 1)
    (hΔ : ∀ n, 0 ≤ Δ n) (hΔle : ∀ n, Δ n ≤ a * q ^ n) :
    Summable Δ ∧ ∑' n, Δ n ≤ a / (1 - q) := by
  have hg : Summable (fun n => a * q ^ n) := (summable_geometric_of_lt_one hq0 hq).mul_left a
  have hs : Summable Δ := Summable.of_nonneg_of_le hΔ hΔle hg
  refine ⟨hs, ?_⟩
  calc ∑' n, Δ n ≤ ∑' n, a * q ^ n := hs.tsum_le_tsum hΔle hg
    _ = a / (1 - q) := by rw [tsum_mul_left, tsum_geometric_of_lt_one hq0 hq, div_eq_mul_inv]

/-- If `Y_n ≥ C G^n` and `L G > 1`, then `L^n Y_n → ∞`. -/
theorem weighted_tendsto (L G C : ℝ) (Y : ℕ → ℝ) (hL : 0 < L) (hC : 0 < C)
    (hLG : 1 < L * G) (hY : ∀ n, C * G ^ n ≤ Y n) :
    Tendsto (fun n => L ^ n * Y n) atTop atTop := by
  have h1 : Tendsto (fun n : ℕ => C * (L * G) ^ n) atTop atTop :=
    (tendsto_pow_atTop_atTop_of_one_lt hLG).const_mul_atTop hC
  refine tendsto_atTop_mono (fun n => ?_) h1
  calc C * (L * G) ^ n = L ^ n * (C * G ^ n) := by rw [mul_pow]; ring
    _ ≤ L ^ n * Y n := mul_le_mul_of_nonneg_left (hY n) (by positivity)

end DyadicNS
