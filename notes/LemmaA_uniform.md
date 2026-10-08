# Lemma A — uniform version (n ≥ N_A = 10^5)

Status 2026-10-08. Everything below is proved with explicit constants, and every number comes from Arb (python-flint)
interval arithmetic. The one exception, η, is stated in §6: the lower bound on the step of L₂ is **not** proved here.
It is the exact inequality that Lemma B must supply. Notation is the same as LemmaA.md.

## 0. Result

**Theorem A-unif.** Let n ≥ N_A := 10⁵ and take k in the window W = {μ−2β < k < μ+2β}. Write θ for the saddle, K = k+1 and q = e^{−Kθ}. Then
  |L_k − L₂,k| ≤ ε_k ≤ Γ_max · θ q,   Γ_max = 0.4282.
For each ν' := Kθ + log θ the sharper cell-wise constant Γ(ν') is the one in ws_p12_regionI.out (θ ≥ 0.002) and
ws_p12_regionII.out (θ ≤ 0.002), and the larger of the two applies. Representative values: Γ ≤ 0.071 at ν' = −2.5, ≤ 0.140 at ν' = 0 (where the crossing is),
≤ 0.225 at ν' = 1, ≤ 0.428 at ν' = 2.1.
**Corollary (for the final logic).** Suppose Lemma B shows, for every k in W,
  (B-η) L₂,k − L₂,k+1 ≥ 2·Γ_max·max(θ_k q_k, θ_{k+1} q_{k+1})  (that is, the step is at least 0.857·θq), and
  (B-sep) |L₂,k| > Γ_max θ_k q_k for k outside the crossing pair.
Then sign L_k is determined for every k ∈ W except at most one.
*Proof.* Call k ∈ W undetermined if |L₂,k| ≤ ε_k. For every other k, |L_k − L₂,k| ≤ ε_k < |L₂,k|, so sign L_k = sign L₂,k.
By (B-sep), every undetermined k lies in the crossing pair {k*, k*+1}. Suppose both were undetermined. Then
  L₂,k* − L₂,k*+1 ≤ |L₂,k*| + |L₂,k*+1| ≤ ε_k* + ε_k*+1 ≤ Γ_max(θ_k* q_k* + θ_k*+1 q_k*+1) ≤ 2Γ_max·max(θ_k* q_k*, θ_k*+1 q_k*+1) = 0.8564·max(θq).
Lemma B (B1) gives L₂,k* − L₂,k*+1 ≥ Σ_min θ_k* q_k* ≥ (Σ_min/1.01)·max(θq) ≥ 0.9100·max(θq) (Σ_min = 0.91918, θ′q′ ≤ 1.01θq), a contradiction.
So at most one k ∈ W is undetermined. Lemma B's (B2) leaves both k* and k*+1 unseparated individually; uniqueness comes from (B1) as above
(REFEREE A.6). The determined signs, L₂ strictly decreasing, and (B3) at the window edges give the pattern +…+ ? −…−. □ The numbers in LemmaA.md (probe at n = 10⁵) give
step/(θq) between 1.02 and 1.21 across the window, so (B-η) holds there numerically with about 20% slack. The ν'-dependent Γ(ν') loosens it further.

## 1. Rescaling (exact)
Put b := Bθ³, s_r := κ_r θ^{r+1} (r = 2..6, s₂ = b), g_r := c_r θ^r = ν^r Li_{1−r}(e^{−ν}) (r ≥ 2) and ν := Kθ.
Then c₁ = 1 + ν e^{−ν'}/(1−q) and q = θ e^{−ν'}. Every quantity in lemA_bound.py depends on B only through w = φ√B.
With A_r := κ_r B^{−r/2} = s_r b^{−r/2} θ^{r/2−1}, C_r := c_r B^{−r/2} = g_r b^{−r/2} θ^{r/2}, C₁ = c₁(θ³/b)^{1/2}, w₀ = 0.4(b/θ)^{1/2},
the bound ε of LemmaA.md §6 equals the same formula with B = 1, (κ_r, c_r) → (A_r, C_r), φ₀ → w₀. The minor-arc term becomes
2π·sup|χ|·(b/(2πθ³))^{1/2}. `ws_p12_eps.core` implements this. It reproduces lemA_bound.lemma_A to 12 digits at
(10⁵, 1358) and (10⁵, 1851) (`python3 ws_p12_eps.py n k`).
The term |M_model − M₂| is replaced by an exact symbolic identity (ws_p12_sym.py): M_model − M₂ = Num/N₀, where Num is an
explicit polynomial with 11 monomials, all of order θ³ (or θ⁵·c₁²). This replaces Arb's cancellation-prone difference.

## 2. (P1) Ingredient bounds
**Lemma P1.1 (Riemann sum).** If f has total variation TV(f) on [0,∞), then |θΣ_{i≤k} f(iθ) − ∫₀^{kθ} f| ≤ θ·TV(f).
*Proof.* Each |θf(iθ) − ∫_{(i−1)θ}^{iθ} f| ≤ θ·Var_{[(i−1)θ,iθ]} f; sum the pieces. □
**Lemma P1.2.** Let f_r(u) = u^r Li_{1−r}(e^{−u}). Then
  TV(f₁..f₆) ≤ 1.000001, 1.000001, 2.000001, 6.086765, 25.691943, 138.425   (certified 2026-10-08 by ws_mech_tv.py).
*Proof (ws_mech_tv.py, independent, imports nothing from the project).* On u ≤ 2.5, f_r and f_r′ are evaluated from the Bernoulli series
(n ≤ 160, rigorous tail via |B_n|/n! ≤ 4/(2π)ⁿ); on u > 2.5 by Arb's polylog with f_r′ = (r f_r − f_{r+1})/u (identity E14).
On 6000 cells of [0, 60]: if Arb certifies a constant sign of f_r′ on the cell, the variation there is |f(b) − f(a)| (f monotone);
otherwise the cell is bisected (depth ≤ 14) and at the bottom Var ≤ width · sup|f_r′| = ∫|f′| bound. Both charges are valid
variation bounds whether or not f turns on the cell. On [60, ∞), f_r is decreasing (each l^{r−1}u^r e^{−lu} decreases for u > r), so Var = f_r(60). □
*History.* The earlier ws_p12_tv.py charged sup−inf (g(b)^r h(a) − g(a)^r h(b)) on sign-undetermined cells, which is not a variation bound
(REFEREE A.1); its constants 1.06265 … 1626.90 happen to exceed the true values but were never certified. ws_p12_tv.py is retired.
**Corollary P1.3.** s_r ∈ F_r(ν−θ) ± θ·TV_r, where F_r(v) = ∫₀^v f_r = r!(ζ(2) − Σ_{j=0}^r v^j/j!·Li_{2−j}(e^{−v})).
The closed form was checked against quadrature. Since kθ = ν−θ, this follows from P1.1. The quantities g_r and c₁ are exact functions of (θ, ν').
**Lemma W (window range).** Write β = √(6n)/π, ℓ = log β and δ = 5ℓ/β. For n ≥ 10⁵ and k ∈ W:
  (1−δ)/β < θ ≤ r/β with r = (1 − (ℓ+2+1/β)/(ζ(2)β))^{−1/2},
  so θ ≤ θ_A = 0.0040939, ν' ∈ [−2.5102, 2.0840] and ν ≥ 3.1158.
*Proof (ws_p12_window.py).* S(t) = Σ_{i≤k} i/(e^{it}−1) is decreasing. t²S(t) = tΣf₁(it), and f₁ is decreasing, so
F₁(Kt) − t ≤ t²S(t) ≤ ζ(2). The upper bound follows from θ²m ≤ ζ(2) and m ≥ n − μ − 2β. For the lower bound, S(θ₋) > m holds when
(2δ−δ²)ζ(2) > T(Kθ₋) + θ₋, with T(v) = ζ(2) − F₁(v) ≤ (v+1)e^{−v}/(1−e^{−v}). This reduces to G(β) > 0 for an explicit G that increases in β for β ≥ e².
Arb gives G(β_A) = 1.046 > 0. The ν' bounds then follow because Kθ + log θ is increasing in θ. The upper ν'-bound ((ℓ+2+1/β)r − ℓ + log r) and the
lower bound (−2 − δ(ℓ−2) + log(1−δ)) are monotone in β in the favourable direction for β ≥ β_A. (Upper: ℓ-coefficient r−1 = O(ℓ/β) decreasing. Lower:
δ(ℓ−2) decreasing for β ≥ e⁴.) □

## 3. (P2) Minor arcs
**Lemma P2.** Suppose θ ≤ 0.0042, ν = Kθ ≥ 3.0 and 0.4θ ≤ φ ≤ π. Then log|χ(φ)| ≤ −c_*/θ with **c_* = 0.15927**.
*Proof (ws_p12_minor.py).* Use −2log|χ| = S(φ) = Σ_{i≤k} log(1 + sin²(iφ/2)/sinh²(iθ/2)), with every summand ≥ 0. Write c = φ/θ.
(a) 0.4 ≤ c ≤ 30. Take T ∈ (0, π/2] and s = sinT/T. For i ≤ 2T/φ we have sin(iφ/2) ≥ s·iφ/2. The summand (u = iθ)
ℓ_c(u) = log(1 + s²c²(u/2)²/sinh²(u/2)) is decreasing in u, so S ≥ θ⁻¹∫_θ^{min(2T/c,ν)} ℓ_c ≥ θ⁻¹∫_{θ_A}^{…}.
ℓ_c is increasing in c and the upper limit is decreasing in c. On each of 300 geometric c-cells [c₁, c₂], evaluate with c₁ in the integrand
and c₂ in the limit. Take the best of 29 values of T. Right-endpoint Riemann sums of the decreasing integrand give lower bounds. Minimum: θS ≥ 0.3185.
(M) 30θ ≤ φ ≤ 0.6. Set d = φ/2 and α ∈ (0, π/2). The indices with dist(id, πℤ) ≥ α have density at least ρ = (π−2α−d)/π ≥ (π−2α−0.3)/π
in every initial segment beyond (1+c₀)/ρ, where c₀ = 2α/d + 1. This is block counting: each block of ⌈π/d⌉ consecutive i loses at most
2α/d + 1 bad indices. The summand ≥ log(1 + sin²α/sinh²(iθ/2)) is decreasing in i, so by Abel summation
S ≥ ρθ⁻¹∫_{u₁}^{ν−θ} log(1+sin²α/sinh²(u/2)) du with u₁ = ((2θ+4α/30)/ρ + θ). At α = 0.525: θS ≥ 1.057.
(b) 0.6 ≤ φ ≤ π. For each pair (2j−1, 2j), one of |sin((2j−1)φ/2)|, |sin(jφ)| is ≥ sin(φ/4) ≥ sin(0.15), since the two angles differ by
φ/2 ≤ π/2 and their sum is not near a multiple of π. Bound that index from below by the larger sinh, sinh(jθ). This gives
S ≥ θ⁻¹∫_{θ_A}^{(ν−θ)/2} log(1+sin²(0.15)/sinh²v) dv. Result: θS ≥ 0.4077.
c_* = min/2 = 0.15927. □ (The c_* ≈ 0.11 target is met. Region checks need c_* ≥ 3.5·θ, which holds with large slack.)

## 4. Region I: θ ∈ [0.002, θ_A] (finite box cover)
In ws_p12_boxes.py, θ runs over geometric cells (ratio 1.01) and ν' ∈ [−2.52, 2.10] runs over cells of width 0.01, clipped to ν ≥ 3.1158 (Lemma W).
On each box, Arb evaluates core() on balls containing every (A_r, C_r, w₀, minor) value allowed by P1.3, the exact g_r, and Lemma P2.
This bounds ε over the whole box, which is then divided by min θq = θ_lo²e^{−ν'_hi}. All asserts pass (δ₀ < 1, D₀ − E_D > 0, Δ < M₂).
Maximum Γ ≤ 0.42024 (at ν' ≈ 2.1; re-run 2026-10-08 with the certified TV and ν′_hi = 2.10 passed explicitly — the cell formula extends the cover to 2.11). Previously 0.42815 with the uncertified TV. Γ_max = 0.4282 is kept as the declared constant, now with slack. Runtime is 27 s. Output: ws_p12_regionI.out.

## 5. Region II: θ ≤ 0.002 (monotone majorant)
In ws_p12_asym.py, write t = √θ. By §1 each of A_r, C_r, C₁ is a monomial t^p q^{b_q} ν^j times a coefficient that does not depend on θ
except through b, s_r and y = e^{−ν}. Those coefficients are bounded uniformly by s_r ≤ r!ζ(2) + θ₁TV_r (the integrand is positive, so F_r ≤ r!ζ(2)),
b ≥ F₂(ν_min − θ₁) − θ₁TV₂, and y ≤ e^{−ν_min}. Every error polynomial has nonnegative coefficients, so it is a nonnegative combination of
t^p q^{b_q} ν^j. After division by θq = t²q, each monomial is t^P ν^j e^{−(b_q−1)ν'}. With ν = ν' + 2log(1/t), t^P ν^j is nondecreasing in t
on (0, √θ₁] whenever P > 0 and ν_min > 2j/P. The script **asserts** this for every monomial and then evaluates at θ₁. Gaussian-tail and minor-arc
terms carry the factors e^{−0.04b/θ}/θ² and e^{−c_*/θ}θ^{−3.5}. Both are increasing for θ < 0.02b and θ < c_*/3.5, which is asserted.
δ₀ is θ-free (w₀²A₄ = 0.16 s₄/b). D₀ ≥ 1 − (5/24)A₃². M₂ is bounded below termwise. Result: **Γ ≤ 0.17605 for all θ ≤ 0.002** (0.1778 before the TV re-certification) and all ν' ∈ [−2.52, 2.10]. Runtime is under 1 s.
(θ → 0 is therefore covered with no further argument. All terms are t^{P>0}, so Γ → 0 as θ → 0.)

## 6. η: what is and is not proved here
Not proved here: a lower bound on the L₂-step. The requirement handed to Lemma B, exactly as stated:
  **(B-η)** for n ≥ 10⁵ and consecutive k, k+1 ∈ W near the crossing: L₂,k − L₂,k+1 ≥ 0.857·max(θ_k q_k, θ_{k+1} q_{k+1});
  **(B-sep)** for every other k ∈ W: |L₂,k| > 0.4282·θ_k q_k. Alternatively, use the cell-wise Γ(ν'_k) from §4/§5 in place of 0.4282.
Heuristic (not proved): d/dk[−log(1−q)] ≈ −θq·(1 + O(νθ)), and the M₂ correction is O(θ²·…), so step ≈ θq. Pointwise Arb values at
n = 10⁵ give step/(θq) = 1.21, 1.10, 1.05, 1.03, 1.02 at x = −2, −1, 0, 1, 2, which is above 0.857 everywhere.

## 7. Scripts
ws_p12_common.py (f_r, F_r, Eulerian polynomials), ws_mech_tv.py (P1.2; ws_p12_tv.py retired), ws_p12_window.py (Lemma W), ws_p12_minor.py (P2),
ws_p12_sym.py (exact M-gap), ws_p12_eps.py (scaled ε, consistency check), ws_p12_boxes.py (Region I), ws_p12_asym.py (Region II).
Reproduce: `python3 ws_mech_tv.py; python3 ws_p12_window.py 100000 5; python3 ws_p12_minor.py 0.0042 3.0;
python3 ws_p12_boxes.py 0.002 0.0040939 -2.52 2.10 3.1158 0.15926 1.01 0.01; python3 ws_p12_asym.py 0.002 -2.52 2.10 0.15926 0.05`.
The TV constants are hard-coded in ws_p12_eps.py, rounded up from the Arb enclosures printed by ws_mech_tv.py.

## 8. Caveats (honest)
- The pointwise LemmaA.md form takes the minor-arc sup in float64. This uniform version replaces that with Lemma P2, which is fully Arb.
- The algebra of §§3–6 of LemmaA.md (T_ω, T_H, tails, combination, M₂ vs Ñ₀/D₀) is now machine-checked by ws_mech_lemA.py (MECH.md §6); ws_p12_eps.core is cross-checked against an independent implementation built from the verified polynomials.
- [Superseded by ws_mech_minor.py, MECH.md §1, which replaces (M)/(b) with the Dirichlet-identity proof.] In Lemma P2 regime (M), the block-counting density claim and the Abel-summation step are argued by hand above. They are the least-checked parts. A conservative
  fallback is regime (b) alone for φ ≥ 0.6, plus regime (a) extended (θS ≥ 0.13 at c = 100 numerically). This needs CA → φ/θ up to 0.6/θ, which is not done.
