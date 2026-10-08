# Lemmas C and D: the tails of the ACG row (n ≥ 10⁵)

Status 2026-10-08. Every inequality below has explicit constants. Every number comes from Arb (python-flint) ball
arithmetic, except ws_cd_sanity.py, which is labelled NON-RIGOROUS and is only a cross-check. Notation follows
INTERFACE.md and LemmaA_uniform.md: β = √(6n)/π, ℓ = log β, μ = βℓ, m = n−k, K = k+1, θ is the saddle, q = e^{−Kθ},
ν = Kθ, ν' = ν + log θ, and **v := kθ = ν − θ**. Also f₁(u) = u/(eᵘ−1), f₂(u) = (u/2)²/sinh²(u/2), with f_r and F_r as in
ws_p12_common.py, and s_r = θ^{r+1}κ_r = θΣ_{i≤k} f_r(iθ), b = s₂ = θ³B.

## 0. Result

**Theorem CD.** Let n ≥ N_A = 10⁵ and k_D(n) := ⌊1.7·n^{1/3}⌋ (≥ 78). Then:
- **(D)** R_k > 1 for 1 ≤ k ≤ k_D(n). This is elementary (§4).
- **(C-left)** R_k > 1 for k_D(n) < k ≤ μ − 2β. Lemma CL-a covers kθ ≤ 0.1 and Lemma CL-b covers kθ ≥ 0.1.
- **(C-right)** R_k < 1 for μ + 2β ≤ k ≤ ⌈n/2⌉. In fact log R_k ≤ −0.614·θ.
- **(E)** a_{k+1} ≤ a_k for ⌈n/2⌉ ≤ k ≤ n−1. This is exact (§5).

The fixed split of INTERFACE.md is used without change, with k_D = ⌊1.7n^{1/3}⌋, which is larger than ⌊n^{1/3}⌋. Lemma CL-a holds for
every k ≥ 47 with kθ ≤ 0.1 and every n, so CL and D overlap: D could stop anywhere between k = 47 and 1.7n^{1/3}.
Margins (re-run 2026-10-08 with the certified TV of ws_mech_tv.py): in CL-b, (q/θ)/(1 + X/(1−θX)) ≥ 5.82 (was 5.816), against the requirement > 1. In CL-a, ρ ≤ 0.161
(0.0673 with the mechanical N₃₄ of ws_mech_cla.py), against the requirement < 0.897. In CR, e^{−ν'}/(1−q) + X ≤ 0.3852, against the requirement < 1. In D, Ψ(s_A) = 0.305, against the requirement > 0.

## 1. Lemma C0: crude saddle bound (the zeroth-order version of Lemma A)

Start from the exact representation (LemmaA.md (2)): R_k/g₀ − 1 = Ñ/D, with g₀ = e^{−θ}/(1−q), D = ∫_{−π}^{π}χ and Ñ = ∫χ(e^H−1).
Here χ = e^ψ is the characteristic function of X−m and H = log G − log g₀. The bounds used are those of LemmaA.md §2:
|ψ^{(r)}| ≤ κ_r, |H'| ≤ c₁ = 1 + Kq/(1−q), |H''| ≤ c₂ = K²q/(1−q)², and Re H ≤ 0.

**Lemma C0.** Let φ₀ = c₀θ and suppose Re ψ(φ) ≤ −(1−δ)Bφ²/2 on |φ| ≤ φ₀ with δ < 1. Put
T₁ = c₁κ₃/(2B²), T₂ = (c₂+c₁²)/(2B), A₃ = κ₃B^{−3/2}, w₀ = φ₀√B and M := √(2πB)·sup_{φ₀≤|φ|≤π}|χ|. Then
  |R_k/g₀ − 1| ≤ ρ := (T₁/(1−δ)^{5/2} + T₂/(1−δ)^{3/2} + 2M) / (1 − erfc(w₀/√2) − √(2/π)A₃/(3(1−δ)²) − M),
provided the denominator is positive.

*Proof.* Normalise integrals by √(2π/B).
(i) **Numerator, central arc.** Since (e^H)'' = (H''+H'²)e^H and Re H ≤ 0, we have e^H − 1 = ic₁φ + E(φ) with |E| ≤ (c₂+c₁²)φ²/2.
Because χ(−φ) = conj χ(φ), ∫_{|φ|≤φ₀} χ·φ = 2i∫₀^{φ₀} φ Im χ. Here Im χ = e^{Re ψ} sin(Im ψ). Im ψ is odd (χ(−φ) = conj χ(φ)) and vanishes to second order at 0 (ψ'(0) = 0 as E X = m; ψ''(0) = −B is real),
and |(Im ψ)'''| ≤ κ₃, so |Im ψ| ≤ κ₃|φ|³/6. Hence |c₁∫χφ| ≤ c₁∫_ℝ e^{−(1−δ)Bφ²/2}κ₃φ⁴/6 = √(2π/B)·T₁/(1−δ)^{5/2}.
Similarly, ∫|χ||E| ≤ √(2π/B)·T₂/(1−δ)^{3/2}.
(ii) **Numerator, minor arcs.** Here |e^H − 1| ≤ 2, which contributes at most 2·2π·sup|χ| = √(2π/B)·2M.
(iii) **Denominator.** On |φ| ≤ φ₀, write ω = ψ + Bφ²/2. Since ω(0) = ω'(0) = ω''(0) = 0 and |ω'''| ≤ κ₃, we get |χ − e^{−Bφ²/2}| ≤ e^{−Bφ²/2}|e^ω−1|.
This is at most e^{−(1−δ)Bφ²/2}κ₃|φ|³/6, using |e^z−1| ≤ |z|max(1,e^{Re z}) and Re ω ≤ δBφ²/2. Integrating gives
√(2π/B)·√(2/π)A₃/(3(1−δ)²). Then ∫_{|φ|≤φ₀}e^{−Bφ²/2} = √(2π/B)(1−erfc(w₀/√2)). The minor arcs contribute at most √(2π/B)M. □

**Two admissible values of δ.**
(δ_T) Re ψ is even, with (Re ψ)''(0) = −B and |(Re ψ)''''| ≤ κ₄. So Re ψ ≤ −Bφ²/2 + κ₄φ⁴/24, and δ_T = κ₄φ₀²/(12B) = c₀²s₄/(12b).
(δ_C), when c₀v < 2: −2Re ψ = Σ_i log(1 + sin²(iφ/2)/sinh²(iθ/2)). For x = iφ/2 ≤ c₀v/2 we have sin²x ≥ x²(1−x²/3), and log(1+λa) ≥ λlog(1+a)
for λ ∈ [0,1]. With a_i = (iφ/2)²/sinh²(iθ/2) ≤ c₀², we have log(1+a) ≥ a·log(1+c₀²)/c₀² and Σa_i = Bφ². This gives
δ_C = 1 − (log(1+c₀²)/c₀²)(1 − (c₀v)²/12).

**Consequences** (with log(1−q) ≤ −q, log(1+ρ) ≤ ρ and log(1−y) ≥ −y/(1−y)):
  log R_k ≤ −θ + q/(1−q) + ρ,   log R_k ≥ −θ + q + log(1−ρ),   log R_k ≥ −θ − log ν + log(1−ρ)  (1−q ≤ ν).

## 2. Lemma C1: θ-uniform cell bound X ≥ ρ/θ  (ws_cd_common.X_cell)

**Lemma C1.** Fix 0 < v₁ < v₂ ≤ ∞, c₀ and a point θ̄. Suppose Ψ satisfies θS(φ) ≥ Ψ (S := −2log|χ|) for all φ ∈ [c₀θ, π], all θ ≤ θ̄
and all v ∈ [v₁, v₂], and suppose Ψ > 5θ̄. Then the number X = X_cell(θ̄, v₁, v₂, c₀) satisfies ρ ≤ θX for every saddle θ ≤ θ̄ with kθ ∈ [v₁, v₂].
*Proof.* The ingredient bounds are as follows.
- s₃ ≤ min(F₃(v₂) + θTV₃, v₂ sup f₃) and s₄ ≤ min(F₄(v₂) + θTV₄, v₂ sup f₄). These come from P1.3 of LemmaA_uniform.md together with s_r ≤ kθ·sup f_r. For v₂ = ∞, F_r(∞) = r!ζ(2).
- b ≥ max(F₂(v₁) − θ, v₁f₂(v₂)), since f₂ is decreasing and θΣf₂(iθ) ≥ ∫_θ^{(k+1)θ}f₂ ≥ F₂(kθ) − θ; also b ≥ kθ·f₂(kθ). And b ≤ F₂(v₂).
- c₁θ = θ + f₁(ν) ≤ θ + f₁(v₁), and c₂θ² = f₂(ν) ≤ f₂(v₁), since ν ≥ v ≥ v₁ and f₁, f₂ are decreasing.

In scaled form, T₁ = θ(c₁θ)s₃/(2b²), T₂ = θ(f₂(ν)+(c₁θ)²)/(2b), A₃ = s₃θ^{1/2}b^{−3/2}, w₀ = c₀(b/θ)^{1/2}, δ_T = c₀²s₄/(12b), and
M ≤ √(2πb)θ^{−3/2}e^{−Ψ/(2θ)}. After substituting the bounds above, T₁/θ, T₂/θ, A₃ and δ are nondecreasing in θ and w₀ is nonincreasing in θ.
M/θ = √(2πb)θ^{−5/2}e^{−Ψ/(2θ)} is increasing in θ when Ψ ≥ 5θ (derivative of −(5/2)logθ − Ψ/(2θ)). So ρ/θ is at most its value at θ̄. sup f₃ ≤ 2.0302009 and
sup f₄ ≤ 6.1253939 come from cell bounds f_r = g^r h_r ≤ g(b)^r h_r(a) on [a,b] ⊂ [0,60] and decrease on [60,∞) (ws_cd_common.sup_f, which is rechecked in __main__). □

X_best minimises over c₀ ∈ {0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0}. Any single admissible c₀ suffices.

### Lemma C-minor: the minor-arc constant Ψ  (ws_cd_common.psi_minor)
S(φ) = Σ_{i≤k} log(1 + sin²(iφ/2)/sinh²(iθ/2)), and every summand is ≥ 0. Write c = φ/θ. The bound is valid for all θ ≤ θ̄ and all k with kθ ≥ v₁.
For c_M ∈ {30, 4π/v₁, 8π/v₁}, Ψ := min(A, min(B1, B2)), where:
- **(A)** c ∈ [c₀, c_M], on 50 geometric cells [c_a, c_b]. For T ∈ (0, π/2] and i ≤ 2T/φ we have sin(iφ/2) ≥ (sinT/T)(iφ/2), so the summand is ≥ log(1+s²c_a²f₂(iθ)), which is
decreasing in iθ. A right Riemann sum then gives θS ≥ ∫_{θ̄}^{min(v₁, 2T/c_b)} log(1+s²c_a²f₂(u))du. The script takes the best of 16 values of T.
- **(B1)** φ ∈ [c_Mθ, 0.6], with d = φ/2. Indices with dist(id, πℤ) < α lie in windows around jπ, j ≤ (Nd+α)/π, and each window holds ≤ 2α/d+1 integers. So
#good(≤N) ≥ ρN − C₀, with ρ = 1 − (2α+0.3)/π and C₀ = (1+α/π)(2α/d+1) ≤ (1+α/π)(4α/(c_Mθ)+1). The j-th good index is ≤ (j+C₀)/ρ, and the
summand there is ≥ w(iθ) := log(1+sin²α/sinh²(iθ/2)), which is decreasing. So θS ≥ ρ∫_{θ(1+C₀)/ρ}^{kθ} w ≥ ρ∫_{(θ̄+θ̄C₀)/ρ+θ̄}^{min(v₁,40)} w. The script takes the best α ∈ {0.10, …, 0.75}.
- **(B2)** φ ∈ [0.6, π]. In each pair (2j−1, 2j), the angles (2j−1)φ/2 and jφ differ by φ/2 ∈ [0.3, π/2]. They cannot both lie within φ/4 of πℤ, so one of them has
|sin| ≥ sin(0.15). Its sinh is ≤ sinh(jθ), so θS ≥ ∫_{θ̄}^{min(v₁/2, 20)} log(1+sin²(0.15)/sinh²w)dw.

All integrals are right-endpoint Riemann sums of decreasing integrands evaluated in Arb, so they are rigorous lower bounds. The (A)/(B2) regimes are the
LemmaA_uniform P2 arguments, re-parametrised by v = kθ instead of ν. (B1) uses a more conservative count than LemmaA_uniform (C₀ carries the factor (1+α/π) and the lower limit includes the θ̄C₀ term). Both counts are valid — see LemmaB §7b, which shows the original limit already contains the θN_b term — so (B1) is a conservative variant, not a correction. Both are in any case superseded by the Dirichlet-identity bound of ws_mech_minor.py / ws_mech_cd.py (MECH.md §1).

## 3. Regions CL and CR

**Lemma CW (parameter ranges; ws_cd_window.py, n ≥ 10⁵).**
(i) θ_k is increasing in k: S(t;k) is increasing in k and decreasing in t, and m = n−k decreases.
(ii) θ²m ≤ ζ(2) (LemmaA_uniform, Lemma W).
CL: k ≤ μ−2β, so by (i) and Lemma W θ ≤ θ_A = r/β ≤ 0.0040939, and ν' ≤ ν_CL := (ℓ−2+1/β)r − ℓ + log r ≤ −1.95361.
CR: μ+2β ≤ k ≤ ⌈n/2⌉, so m ≥ (n−1)/2 and θ ≤ θ_CR := √(2ζ(2)/(n−1)) ≤ 0.0057358. The proof of Lemma W gives θ > (1−δ)/β, δ = 5ℓ/β.
That proof only needs m ≤ n and Kθ₋ ≥ (ℓ−2)(1−δ), and S(θ₋;k)−m only grows with k. Hence ν' ≥ ν_CR := (ℓ+2)(1−δ) − ℓ + log(1−δ) ≥ 1.04305.
Monotonicity in β: ν_CL is decreasing, θ_A and θ_CR are decreasing, and ν_CR = 2 − δ(ℓ+2) + log(1−δ) is increasing (δ and δℓ are decreasing for ℓ ≥ 2). So the
values at n = 10⁵ hold for all n ≥ 10⁵. □

**Lemma CR (ws_cd_CR.py).** For n ≥ 10⁵ and μ+2β ≤ k ≤ ⌈n/2⌉: v = ν' − logθ − θ ≥ ν_CR + log(1/θ_CR) − θ_CR ≥ 6.198, and
X = X_cell(θ_CR, 6.198, ∞, 0.4) ≤ 0.03219 (Ψ = 0.3264, M ≤ 4.6·10⁻⁹). Since q/θ = e^{−ν'} ≤ e^{−ν_CR} and q ≤ e^{−v},
  log R_k ≤ θ(e^{−ν'}/(1−q) + X − 1) ≤ θ(0.3853 − 1) < 0. □
Large k needs no separate treatment. For k up to n/2 we have v = kθ → ∞ while θ ≤ θ_CR, and the last cell [v₁, ∞) uses only F_r(∞) = r!ζ(2).
So nothing degenerates as q → 0, and the bound is uniform up to k = ⌈n/2⌉.

**Lemma CL-a (ws_cd_CLa.py; n-free).** If k ≥ K₀ = 47 and v = kθ ≤ V₀ = 0.1, then R_k > 1, for every n.
*Proof.* Here θ = v/k ≤ θ₀ := V₀/K₀, and the bounds b ≥ vf₂(V₀), s₃ ≤ v·sup f₃, s₄ ≤ v·sup f₄, c₁θ ≤ 1+θ₀, f₂ ≤ 1 hold. They give
T₁ ≤ (1+θ₀)sup f₃/(2K₀f₂(V₀)²), T₂ ≤ (1+(1+θ₀)²)/(2K₀f₂(V₀)), A₃ ≤ sup f₃·K₀^{−1/2}f₂(V₀)^{−3/2}, w₀ ≥ c₀(K₀f₂(V₀))^{1/2}, and δ = min(δ_T, δ_C).
The normalised minor-arc integral N := √(B/2π)∫_{φ₀≤|φ|≤π}|χ| replaces M in Lemma C0. Its proof only used integrals, so N ≤ N₁+N₂+N₃+N₄:
- N₁: c ∈ [c₀, 10] in 60 cells. Since kφ/2 ≤ c_{j+1}V₀/2 =: T_j < π/2, each summand is ≥ λ_j·(iφ/2)²/sinh²(iθ/2) with λ_j = log(1+s_j²c_{j+1}²)/c_{j+1}².
So |χ| ≤ e^{−λ_jBφ²/2}, and N₁ ≤ Σ_j erfc(c_j(λ_jK₀f₂(V₀)/2)^{1/2})/√λ_j, using θ²B = b/θ ≥ kf₂(v).
- N₂: φ ∈ [10θ, π/k]. Here sin(iφ/2) ≥ (2/π)(iφ/2), so |χ| ≤ ((2/π)²f₂(V₀))^{−k/2}c^{−k}, and N₂ ≤ √(2/π)√k((2/π)²f₂(V₀))^{−k/2}10^{1−k}/(k−1), using θ√B ≤ √k.
- N₃: φ ∈ [π/k, 0.6]. Block counting with d ≥ π/(2k) gives ≥ kγ−2 good indices, γ = 1 − (2α+0.3)/π − (1+α/π)4α/π, α = 1/4. Each has summand
≥ 2log(2sinα/(vσ₀)), with σ₀ = sinh(V₀/2)/(V₀/2). So N₃ ≤ √(2π)k^{3/2}v^{−1}(vσ₀/(2sinα))^{(kγ−2)/2}, using √B ≤ k^{3/2}/v.
(The exponent is halved, which is conservative because the base is < 1.)
- N₄: φ ∈ [0.6, π]. Pairs give N₄ ≤ √(2π)k^{3/2}v^{−1}(vσ₀/(2 sin 0.15))^{⌊k/2⌋/2}.

N₃ and N₄ are nondecreasing in v, since their exponents are ≥ 1 (asserted), so they are evaluated at V₀. N₁ already uses K₀. N₂, N₃ and N₄ are decreasing in k: the ratios
N(k+1)/N(k) are < 1 (asserted), and for N₄ the 2-step ratio is < 1, so the script takes max(N₄(K₀), N₄(K₀+1)).
Result at c₀ = 0.5: N ≤ 0.0420 and ρ ≤ 0.16088. Then log R_k ≥ −θ − log ν + log(1−ρ), with ν ≤ V₀(1+1/K₀), which is > 0 because 1−ρ > e^{θ₀}V₀(1+1/K₀) = 0.1024. □

**Lemma CL-b (ws_cd_CLb.py).** For n ≥ 10⁵ and k_D < k ≤ μ−2β with v = kθ ≥ 0.1, R_k > 1.
*Proof.* Cover v ∈ [0.1, ∞) by 33 geometric cells [v₁, v₂] (ratio 1.15) plus a final cell [8.76, ∞). On each cell, θ ≤ θ̄ := min(θ_A, v₂/78), since k > k_D ≥ 78.
The cell constant is X = X_best(θ̄, v₁, v₂). Then log R_k ≥ −θ + q − θX/(1−θX) > 0 iff q/θ > 1 + X/(1−θX). For q/θ we use the larger of
e^{−ν'} ≥ e^{−ν_CL} = 7.054 and e^{−ν}/θ ≥ e^{−(v₂+θ̄)}/θ̄. Every cell passes, and the minimum margin of (q/θ)/(1+X/(1−θ̄X)) is 5.816 (cell v ∈ [3.29, 3.79]). The full
table is in ws_cd_CLb.out. □

C-left therefore holds: every k ∈ (k_D, μ−2β] has k ≥ 79 ≥ 47, and it falls under CL-a (v ≤ 0.1) or CL-b (v ≥ 0.1).

## 4. Lemma D (elementary; ws_cd_D.py)

**D1.** C(n−1,k−1)/k! ≤ p(n,k) ≤ C(n−1+k(k−1)/2, k−1)/k!.
*Proof.* Lower bound: there are C(n−1,k−1) compositions of n into k positive parts, and each partition arises from at most k! of them.
Upper bound: λ ↦ (λ_i + k − i)_i maps partitions of n into k parts injectively to partitions of N = n + k(k−1)/2 into k *distinct* parts.
Each of those arises from exactly k! compositions of N, and there are C(N−1,k−1) such compositions. □
**D2.** If Φ(k,n) := log(n−k) − log(k(k+1)) − k(k−1)²/(2(n−k+1)) > 0, then p(n,k+1) > p(n,k).
*Proof.* By D1 it suffices that C(n−1,k)/(k+1)! > C(N−1,k−1)/k!. Now C(n−1,k)/C(n−1,k−1) = (n−k)/k, and
C(n−1,k−1)/C(N−1,k−1) = Π_{j=0}^{k−2}(n−1−j)/(n−1−j+D) with D = k(k−1)/2. That product is ≥ exp(−Σ D/(n−1−j)) ≥ exp(−(k−1)D/(n−k+1)). □
**D3.** For n ≥ 10⁵ and 1 ≤ k ≤ 1.7n^{1/3}, Φ > 0.
*Proof.* Φ is decreasing in k, since each term is. Put s = n^{1/3} and γ = 1.7. Then Φ(γs, s³) ≥ Ψ(s) := log((s²−γ)/(γ(γs+1))) − γ³/(2(1−γ/s²)), using k(k−1)² ≤ k³ and n−k+1 > n−k.
Ψ is increasing in s, because 2s(γs+1) > γ(s²−γ) and the last term is decreasing. (Exact form, PROOFS_FULL.md §7 / ws_close_D.py: Ψ′(s) = (γs²+2s+γ²)/((s²−γ)(γs+1)) + γ⁴s/(s²−γ)² > 0 for s > √γ.) Arb gives Ψ(10^{5/3}) = 0.3046 > 0. □
For comparison, the full binomial criterion is sharper: at n = 10⁵ it holds for k ≤ 82, and at n = 10⁹ for k ≤ 2202 ≈ 2.2n^{1/3}. Its true reach is c·n^{1/3}
with c slowly growing (√(2 log) behaviour is not attained; the k³/n term limits it). D3 already overlaps CL-a, which starts at k = 47.

## 5. Region E (exact)
For k ≥ n/2, p(n,k) = p(n−k). Subtracting 1 from every part is a bijection onto partitions of n−k with at most k parts, and that constraint is vacuous since n−k ≤ k.
Also p(j) ≤ p(j+1), because appending a part 1 is an injection. Hence for ⌈n/2⌉ ≤ k ≤ n−1, a_{k+1} = p(n−k−1) ≤ p(n−k) = a_k. The inequality is strict except at k = n−1,
where p(0) = p(1).

## 6. Sanity check (NON-RIGOROUS, float64; ws_cd_sanity.py)
Exact column recurrence for n = 10⁵ (k ≤ 4000) and n = 2·10⁵ (k ≤ 5000). There are no sign violations in the CL/CR ranges. The measured |R/g₀−1|/θ lies below the rigorous
cell X at every sampled k (k = 47 … 4999). For example, at n = 10⁵ the measured/rigorous values are 0.0010/254.8 at k = 47, 0.0255/0.187 at k = 864 = μ−2β, and 0.0009/0.0118 at
k = 1852 = μ+2β. The bounds are valid and conservative by a factor of about 7–10⁵. Outputs: ws_cd_sanity_2e5.out (n = 10⁵ was printed interactively; rerun with `python3 ws_cd_sanity.py 100000 4000`).

## 7. Scripts and reproduction (all < 1 min except ws_cd_sanity: about 15 s / 60 s)
ws_cd_common.py (Lemma C0/C1 evaluator, C-minor, sup f_r), ws_cd_window.py (Lemma CW), ws_cd_CR.py, ws_cd_CLa.py, ws_cd_CLb.py (49 s), ws_cd_D.py, and
ws_cd_sanity.py (non-rigorous). The .out files hold the recorded runs.
`python3 ws_cd_window.py 100000 5; python3 ws_cd_CR.py 100000; python3 ws_cd_CLa.py 47 0.1; python3 ws_cd_CLb.py 100000 0.1 78; python3 ws_cd_D.py`

## 8. Dependencies and caveats
- Imported from LemmaA_uniform.md: the TV constants (P1.2, now certified by ws_mech_tv.py), F_r, Lemma W's θ-lower-bound argument (used at n = 10⁵ with C = 5), and the exact representation (2) of LemmaA.md.
Lemma C0 re-derives its own error terms from scratch (§1) and does not use the T_ω/T_H algebra.
- [Superseded] The minor-arc block-counting (B1) and pair (B2) arguments were hand proofs; ws_mech_cd.py / ws_mech_cla.py replace them by the Dirichlet-identity bound (MECH.md §1) and reproduce PASS for CR and CL-b, so they are no longer load-bearing.
CR does not depend on B1 at all (v ≥ 6.2 makes c_M = 30 and B2/A dominate), but CL-b does in its middle cells.
- Floating-point grids are used only to choose cells. Every bound is evaluated on arb balls that contain the whole cell, and monotonicity in θ is proved in Lemma C1.
