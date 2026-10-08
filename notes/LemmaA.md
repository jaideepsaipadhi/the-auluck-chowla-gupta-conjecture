# Lemma A — explicit saddle-point remainder for the ratio R_k = p(n,k+1)/p(n,k)

> **SUPERSEDED (2026-10-08) as a proof component.** Only §§1–6 (representation, derivative bounds, model, error algebra, combination) are used,
> and they are machine-checked by ws_mech_lemA.py. The pointwise Theorem A, its float64 minor-arc supremum, and the table below (which fails at n = 10⁴)
> are **not** part of the proof; LemmaA_uniform.md replaces them with the uniform Arb version and Lemma P2.

Status (2026-10-08): **pointwise rigorous** (every (n,k) at which it is evaluated, in Arb ball
arithmetic; minor-arc sup in float64 with explicit safety margin). Holds near the crossing for
all sampled n ≥ 20 000. **Not yet uniform in n** — see §7.

## 1. Probabilistic representation
Fix 2 ≤ k ≤ n−2, m = n−k, K = k+1. Let θ > 0 be the unique solution of
  m = Σ_{i≤k} i/(e^{iθ}−1),   x := e^{−θ},  q := x^K.
Let Y_1,…,Y_k be independent, P(Y_i = y) = (1−x^i) x^{iy} (y ≥ 0), and X := Σ_{i≤k} i Y_i.
Then P(X = N) = f_k(N) x^N / F_k(x) (expand ∏(1−x^i)/(1−x^i z^i)), and E X = m (saddle equation).
Since F_{k+1} = F_k/(1−x^K) termwise,
  **R_k = a_{k+1}/a_k = Σ_{s≥0} x^{1+sK} · P(X = m−1−sK) / P(X = m).**   (1)
With χ(φ) := E e^{iφ(X−m)} and P(X = m−j) = (1/2π)∫_{−π}^{π} χ(φ) e^{ijφ} dφ,
  R_k = ∫χ G / ∫χ,  G(φ) := Σ_{s≥0} x^{1+sK} e^{i(1+sK)φ} = x e^{iφ}/(1 − q e^{iKφ}).
Put g0 := G(0) = x/(1−q) and H(φ) := log G(φ) − log g0 = iφ − log(1−qe^{iKφ}) + log(1−q). Then
  **R_k/g0 − 1 = Ñ/D,  D := ∫_{−π}^{π} χ,  Ñ := ∫_{−π}^{π} χ·(e^{H} − 1).**   (2)
(Form (2), rather than N/D, is essential: a relative error δ in D then costs only δ·|R/g0−1| ≈ δ·θ².)

## 2. Uniform derivative bounds (valid for ALL real φ)
Cumulants κ_r := Σ_{i≤k} i^r Li_{1−r}(x^i) (κ₂ = B). ψ := log χ (factorwise principal logs) satisfies
  ψ(0)=ψ′(0)=0,  ψ^{(r)}(0) = i^r κ_r,  **|ψ^{(r)}(φ)| ≤ κ_r for all φ ∈ ℝ, r ≥ 2.**
Proof: −log(1 − y e^{iaφ}) = Σ_l y^l e^{ilaφ}/l, so its r-th φ-derivative is (ia)^r Σ_l l^{r−1} y^l e^{ilaφ},
of modulus ≤ a^r Li_{1−r}(y); sum over factors (y = x^i, a = i). □
Likewise, with c₁ := 1 + Kq/(1−q), c_r := K^r Li_{1−r}(q) (r ≥ 2):
  H(0)=0, H′(0)=i c₁, H″(0)=−c₂, H‴(0)=−i c₃, H⁗(0)=c₄,  **|H′| ≤ c₁, |H^{(r)}| ≤ c_r (r≥2), Re H ≤ 0** (|G| ≤ g0).

## 3. Model and main term
On ℝ, with ω₃ = −iκ₃φ³/6, ω₄ = κ₄φ⁴/24, ω₅ = iκ₅φ⁵/120 and H₁ = ic₁φ, H₂ = −c₂φ²/2, H₃ = −ic₃φ³/6, H₄ = c₄φ⁴/24:
  Q := 1 + ω₃ + ω₄ + ω₅ + ω₃²/2 + ω₃ω₄ + ω₃³/6,   P_H − 1 := H₁+H₂+H₃+H₄ + H₁²/2,
  D₀ := ∫_ℝ e^{−Bφ²/2} Q,   Ñ₀ := ∫_ℝ e^{−Bφ²/2} Q (P_H − 1)   (exact, via Gaussian moments; Ñ₀/D₀ ∈ ℝ).
Closed-form main term used downstream (Lemma B must use THIS main term):
  N₀ = 1 + κ₄/(8B²) − 5κ₃²/(24B³)
  E₁ = [κ₃/(2B²) − κ₅/(8B³) + 105κ₃κ₄/(144B⁴) − 945κ₃³/(1296B⁵)] / N₀
  E₂ = (1/B)(1 + κ₄/(2B²) − 5κ₃²/(4B³))
  **M₂ = 1 + c₁E₁ − (c₂ + c₁²)E₂/2 − 5c₃κ₃/(12B³) + c₄/(8B²),   L₂,k := −θ − log(1−q) + log M₂.**
(Here h′ = −c₁, h″ = c₂, h‴ = −c₃, h⁗ = c₄ in the notation of the notes.)

## 4. Central arc |φ| ≤ φ₀ (φ₀ = 0.4θ)
Write χ = e^{−Bφ²/2} e^{ω}. From §2: Re ω ≤ κ₄φ⁴/24 + κ₆φ⁶/720 ≤ δ₀Bφ²/2 with
δ₀ := (κ₄φ₀²/12 + κ₆φ₀⁴/360)/B < 1. With a_r := κ_r/r! and t = |φ|:
  |e^{ω} − Q| ≤ T_ω(t) e^{δ₀Bt²/2},
  T_ω(t) = a₆t⁶ + (a₄²/2 + a₃a₅)t⁸ + a₄a₅t⁹ + (a₅²/2)t¹⁰ + ½(a₄t⁴+a₅t⁵)(a₃t³+a₄t⁴+a₅t⁵)² + (a₃t³+a₄t⁴+a₅t⁵)⁴/24
  (from e^{ω}−e^{ω̂} ≤ |ω−ω̂| e^{max Re}, |ω−ω̂| ≤ a₆t⁶; the omitted quadratic/cubic terms of e^{ω̂}; and the
  quartic Taylor remainder |e^z − Σ_{r<4} z^r/r!| ≤ |z|⁴/24 · max(1, e^{Re z})).
  |e^{H} − P_H| ≤ T_H(t) := c₅t⁵/120 + (c₁c₂/2 + c₁³/6) t³
  (|e^H−1−H−H²/2| ≤ |H|³/6 as Re H ≤ 0; |H − ΣH_r| ≤ c₅t⁵/120; |H²−H₁²|/2 ≤ (c₂t²/2)(2c₁t)/2).
  |e^{H} − 1| ≤ c₁ t,  |Q| ≤ Q_abs(t) (coefficient-wise moduli).
Hence (integrals over |φ| ≤ φ₀):
  |∫χ(e^H−1) − ∫e^{−Bφ²/2}Q(P_H−1)| ≤ ∫e^{−(1−δ₀)Bt²/2} c₁t T_ω(t) + ∫ e^{−Bt²/2} Q_abs(t) T_H(t),
  |∫χ − ∫e^{−Bφ²/2}Q| ≤ ∫ e^{−(1−δ₀)Bt²/2} T_ω(t),
evaluated exactly with E|Z|^p = (2/B′)^{p/2} Γ((p+1)/2)/√π.

## 5. Tails and minor arcs
Model tails: ∫_{|t|>φ₀} e^{−Bt²/2}P(t) ≤ e^{−Bφ₀²/4} ∫ e^{−Bt²/4} P(t).
Minor arcs φ₀ ≤ |φ| ≤ π: |χ|² = ∏_i (1 + 4y_i sin²(iφ/2)/(1−y_i)²)^{−1}; |∫χ(e^H−1)| ≤ 2·2π·M_χ, |∫χ| ≤ 2π·M_χ,
M_χ := sup_{φ₀≤φ≤π} |χ| computed by subdividing [φ₀, π] and using, per subinterval and per i, the exact
minimum of sin²(iφ/2) (0 if the interval i[a,b]/2 contains a multiple of π, else the endpoint minimum).

## 6. Theorem A (pointwise form)
With E_N, E_D the sums of the central, tail and minor-arc errors above (normalised by √(2π/B)):
  |R_k/g0 − (1 + Ñ₀/D₀)| ≤ Δ_model := (E_N + |Ñ₀/D₀| E_D)/(|D₀| − E_D),
  **|L_k − L₂,k| ≤ ε_k := −log(1 − (Δ_model + |1 + Ñ₀/D₀ − M₂|)/M₂).**
Implementation: lemA_bound.py (Arb/python-flint; saddle θ enclosed by a verified sign change).

### Verified values (rigorous upper bounds)
Near the crossing k* (η := per-step decrement of L₂ at k*):

| n | k* | η | ε (max over k*−1..k*+1) | ε/η |
|---|---|---|---|---|
| 10 000 | 340 | 1.88e-4 | 6.4e-4 | 3.40 (fails) |
| 20 000 | 519 | 9.10e-5 | 3.16e-5 | **0.347** |
| 40 000 | 788 | 4.44e-5 | 9.5e-6 | 0.214 |
| 60 000 | 1003 | 2.94e-5 | 4.8e-6 | 0.163 |
| 100 000 | 1358 | 1.74e-5 | 2.0e-6 | 0.117 |
| 300 000 | 2586 | 5.70e-6 | 3.3e-7 | 0.058 |
| 1 000 000 | 5191 | 1.69e-6 | 4.6e-8 | 0.027 |

Requirement for "≤ 1 undetermined k": ε < η/2. Across the window |x| ≤ 2: ε/|L₂| ≤ 0.016 (n = 20 000),
≤ 0.0044 (n = 60 000), largest near the crossing.

## 7. What remains to make Lemma A uniform in n ≥ N_A (N_A ≈ 2·10⁴)
(P1) **Ingredient bounds**: explicit two-sided bounds, valid for all n ≥ N_A and k in the window, on
  θ, ν := Kθ (q = e^{−ν}), and the scaled quantities Bθ³, κ_rθ^{r+1} (r = 3..6), c_rθ^r — via
  sum ↔ integral comparisons with explicit remainders. Every term of ε is then ≤ θ^{3.5}·poly(ν) and
  η ≥ θ²·(explicit), giving ε/η ≤ explicit function of θ, checked for θ ≤ θ_A.
(P2) **Analytic minor-arc lemma**: Re ψ(φ) ≤ −c/θ on [φ₀, π] with explicit c (numerically c ≈ 0.22 at φ₀ = 0.4θ;
  small φ via (u/2)²/sinh²(u/2) Riemann-sum bound; large φ via the consecutive-pair argument
  |sin(iφ/2)| or |sin((i+1)φ/2)| ≥ sin(φ/4)).
(P1) is shared with Lemma B. Lemma B must analyse L₂ (not L̃).
