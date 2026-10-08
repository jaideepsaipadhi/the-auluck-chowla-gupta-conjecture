# Lemma B — crossing structure of the main term L₂ on the window (n ≥ 10⁵)

Status 2026-10-08. Everything below is proved with explicit constants. Every number comes from Arb (python-flint) ball
arithmetic. Notation and conventions are those of LemmaA.md and LemmaA_uniform.md: W = {μ−2β < k < μ+2β}, θ = θ_k is the saddle
for (m, k) = (n−k, k), K = k+1, ν = Kθ, q = e^{−ν}, ν' = ν + log θ, s_r = θ^{r+1}κ_r, b = s₂, and
L₂,k = −θ − log(1−q) + log M₂, with M₂ exactly as in LemmaA.md §3 (B-scaled form in ws_b_core.M2_formula).
Primes denote the point k+1: θ' = θ_{k+1}, which is the saddle for (m−1, k+1).

## 0. Result

**Theorem B.** Let n ≥ 10⁵ and let k, k+1 ∈ W. Then
- **(B1)** L₂,k − L₂,k+1 ≥ Σ_min·θ_k q_k with Σ_min = **0.91918** (0.91905 before the 2026-10-08 TV re-certification), and θ_{k+1}q_{k+1} ≤ 1.01·θ_k q_k. So
  L₂,k − L₂,k+1 ≥ 0.9099·max(θ_k q_k, θ_{k+1}q_{k+1}) > 0.857 = 2Γ_max. L₂ is strictly decreasing on W.
- **(B2)** Let k* be the crossing (L₂,k* > 0 ≥ L₂,k*+1). For every k ∈ W \ {k*, k*+1}, |L₂,k| > 0.9099·θ_k q_k > Γ_max θ_k q_k.
- **(B3)** At the left end of W, L₂/θ ≥ 5.58. At the right end, L₂/θ ≤ −0.505. Hence there is exactly one sign change in W.

Combined with Theorem A-unif (|L_k − L₂,k| ≤ 0.4282·θ_k q_k), the hypotheses (B-η) and (B-sep) of the Corollary in
LemmaA_uniform.md §0 hold for every n ≥ 10⁵ with the uniform Γ_max. The per-cell Γ(ν') is not needed. The margin in (B1) is
0.9099/0.8564 ≈ 6%. The margin in (B2) is a factor of 2.1.

*Deduction of (B2) and the max in (B1) from Σ_min.* ν_{k+1} > ν_k (because θ' > θ, Lemma B.1), so q_{k+1} < q_k. Also θ' ≤ 1.01θ (Lemma B.1).
Hence θ'q' ≤ 1.01θq. For k ≤ k*−1: L₂,k ≥ L₂,k*+ (L₂,k − L₂,k+1) > Σ_min θ_k q_k. For k ≥ k*+2:
L₂,k ≤ L₂,k*+1 − (L₂,k−1 − L₂,k) ≤ −Σ_min θ_{k−1}q_{k−1} ≤ −(Σ_min/1.01)θ_k q_k. □

## 1. Decomposition of the step
L₂,k − L₂,k+1 = (θ' − θ) + log((1−q')/(1−q)) − (log M₂' − log M₂).
- θ' > θ (B.1), so the first term is ≥ 0. It is dropped.
- q' = e^{−(K+1)θ'} ≤ e^{−(K+1)θ} = q e^{−θ}, so the second term is ≥ T₂ := log(1 + q(1−e^{−θ})/(1−q)).
  This gives T₂ = θq(1 + O(θ + q)), which is the main term.
Hence **Σ := step/(θq) ≥ T₂/(θq) − |Δ log M₂|/(θq)**. Everything is about bounding Δ log M₂ = log M₂' − log M₂.

## 2. Controlling the move (θ, ν, s_r) → (θ', ν', s_r')   [the dθ/dk control]
Let S_k(t) = Σ_{i≤k} i/(e^{it}−1), so S_k(θ) = m and S_{k+1}(θ') = m−1. S_k is decreasing, with −S_k'(t) = B_k(t) = b(k,t)/t³.
By P1.1/P1.3, b(k,t) = t Σ_{i≤k} f₂(it) ≥ F₂(kt) − t·TV₂. Fix τ = 0.01 and put
  b_min := F₂(ν−θ) − θ(1+τ)TV₂ (a lower bound for b(k,t) on t ∈ [θ, θ(1+τ)], since kt ≥ ν−θ and F₂ ↑),   c₁ = 1 + Kq/(1−q).

**Lemma B.1 (step localisation).** If τ·b_min > (1+τ)³c₁θ², then θ < θ' ≤ (1+τ)θ.
*Proof.* S_{k+1}(θ) = m + K/(e^{Kθ}−1) = m + c₁ − 1 > m − 1, so θ' > θ. On [θ, (1+τ)θ] we have −S_{k+1}' ≥ B_k(t) ≥ b_min/((1+τ)θ)³. So
S_{k+1}((1+τ)θ) ≤ m + c₁ − 1 − τθ·b_min/((1+τ)³θ³) < m − 1. □

**Lemma B.2.** Under B.1, dθ := θ' − θ ≤ (1+τ)³c₁θ³/b_min.
*Proof.* MVT: c₁ = S_{k+1}(θ) − S_{k+1}(θ') = B_{k+1}(ξ)·dθ ≥ (b_min/((1+τ)θ)³)·dθ. □

**Lemma B.3.** dν := ν_{k+1} − ν_k = θ' + K·dθ ≤ (1+τ)θ + (ν/θ + 1)dθ. Also dν > 0.

**Lemma B.4.** For r = 2..5: |s_r' − s_r| ≤ dθ·((1+r)s̄_r + s̄_{r+1})/θ + (1+τ)θ·sup_{u∈[ν,ν+dν]} f_r(u),
where s̄_r := r!ζ(2) + (1+τ)θ·TV_r ≥ s_r(k,t) for all t ≤ (1+τ)θ.
*Proof.* s_r(k+1,θ') − s_r(k,θ) = [s_r(k,θ') − s_r(k,θ)] + θ' f_r(Kθ'), and Kθ' ∈ [ν, ν_{k+1}].
Since u f_r'(u) = r f_r − f_{r+1}, we get ∂_t s_r(k,t) = ((1+r)s_r − s_{r+1})/t. Bound this by its modulus with t ≥ θ. □

## 3. Region I: θ ∈ [0.002, θ_A] (box verification) — ws_b_step.py
The cover is the same as in LemmaA_uniform.md §4. θ runs over geometric cells of ratio 1.02 and ν' ∈ [−2.52, 2.10] over cells of width 0.02, clipped to ν ≥ 3.1158 (Lemma W).
This makes 8526 boxes. On each box:
1. Arb asserts the hypothesis of B.1 and computes dθ, dν and Δs_r as balls (B.2–B.4).
2. It forms a convex box H in the six variables (θ, ν, b, s₃, s₄, s₅) that contains both endpoints. The P1.3 enclosures F_r(·) ± θ(1+τ)TV_r are
   evaluated over [ν−θ, ν+dν], so H covers the point (k,θ) and the point (k+1,θ').
3. It evaluates ∇M₂ over H by forward-mode automatic differentiation in Arb (class D in ws_b_core.py), together with M₂(H) > 0.
4. By the mean value theorem on the segment ⊂ H: |Δ log M₂| ≤ Σ_j sup_H|∂_j M₂/M₂|·|Δv_j|.
5. T₂/(θq) is evaluated on the box.
**Result: Σ ≥ 0.91918 on every box** (re-run with certified TV; was 0.91905). The minimum is at ν' ≈ 2.1 and θ = θ_A. Runtime is 5 s.
The θ-increment dominates the M₂-variation (the ∂_ν M₂·dν term), and it grows with ν.

## 4. Region II: θ ≤ 0.002 (monotone majorant) — ws_b_asym.py
This uses the graded-monomial method of LemmaA_uniform.md §5 with t = √θ. Every quantity at both points k and k+1 is majorised by nonneg.
combinations of t^p q^{b_q} ν^j:
- θ'^e ≤ (1+τ)^e θ^e for e > 0 and ≤ θ^e for e < 0.
- q' ≤ q.
- ν' ≤ ρν with ρ = 1 + (1+τ)θ₁/ν_min + (1+τ)³(c₁θ²)_max(1+θ₁/ν_min)/b_min.
- f_r(u) ≤ (ρν)^r q P_{r−1}(y)/(1−y)^r with y = e^{−ν_min}.
- |∂_u g_r| ≤ (r g_r + g_{r+1})/u, and |h'(u)| ≤ (u+1)e^{−u}/(1−e^{−u})² for c₁ = 1 + h(ν)/θ.
The increments of A_r, C_r and C₁ come from the product rule with B.2–B.4. The increments of E₁, E₂ and N_n come from the quotient rule with N_n ≥ 1 − (5/24)A₃²_max.
After division by θq = t²q, every monomial has P > 0 and ν_min > 2j/P. The script asserts this, and the monomial is then evaluated at θ₁ = 0.002.
The same machinery bounds |log M₂|/θ ≤ |M₂−1|/(θ(1−|M₂−1|)), which is used in §5.
T₂/(θq) ≥ (1−θ₁/2)(1 − θ₁q₁/(2(1−q₁))), using x ≥ θq(1−θ/2), x ≤ θq/(1−q) and log(1+x) ≥ x(1−x/2).
**Result: Σ ≥ 0.96036 for all θ ≤ 0.002 and all ν' ∈ [−2.52, 2.10]; |log M₂|/θ ≤ 0.14236.** Runtime is under 2 s.

## 5. (B3) Window edges — ws_b_edges.py
With ε_β = (ℓ+2+1/β)/(ζ(2)β), r = (1−ε_β)^{−1/2} and δ = 5ℓ/β, Lemma W gives:
- Left end (K ≤ μ−2β+2, θ ≤ r/β): ν'_L ≤ −2 + (r−1)(ℓ−2) + 2r/β + log r ≤ −2 + ε_β(ℓ−1)/(1−ε_β) + 2/(β(1−ε_β)) ≤ **−1.906674**.
  (Corrected 2026-10-08: the printed majorant's supremum over all β ≥ β_A is −1.9066744, certified by ws_mech_window.py with no monotonicity
  argument; the former digits −1.9067 were 3·10⁻⁵ too small. The sharp middle expression is ≤ −1.94916 for all β ≥ β_A.)
- Right end (K > μ+2β, θ > (1−δ)/β): ν'_R ≥ 2 − δ(ℓ+2) + log(1−δ) ≥ **1.0431**. Certified for all β ≥ β_A by ws_mech_window.py (≥ 1.0430522), no monotonicity argument needed.
L₂/θ is then bounded as follows. For θ ≤ 0.002: L₂/θ = −1 + (−log(1−q))/θ + log M₂/θ with q/θ = e^{−ν'}, so left ≥ −1 + e^{1.906674} − 0.1424 = 5.588 and
right ≤ −1 + e^{−1.043}/(1−θ₁e^{−1.043}) + 0.1424 = −0.505. For θ ∈ [0.002, θ_A], Arb boxes give left ≥ 5.62 and right ≤ −0.643 (ws_mech_edges.py re-runs both regions with −1.906674).
Since L₂ is strictly decreasing on W by (B1), it has exactly one sign change in W. □

## 6. Consistency checks (non-essential) — ws_b_check.py
At n = 10⁵, M₂ from M2_formula agrees with ws_p12_eps.core to within 1e-13. The exact step/(θq) is 1.0767, 1.0519 and 1.0366 at k = 1210, 1358 (crossing) and 1500.
Σ_min = 0.919 is therefore a conservative lower bound (about 13% below the true value).

## 7. AUDIT of Lemma P2, regime (M): 30θ ≤ φ ≤ 0.6 — ws_b_minorM.py
**Re-derivation.** Put d = φ/2 ∈ [15θ, 0.3], take α ∈ (0, π/2) with 2α + d < π, and let w_i := log(1 + sin²α/sinh²(iθ/2)), which is decreasing in i.
Call i *good* if dist(id, πℤ) ≥ α. Every summand of S is ≥ 0, so S ≥ Σ_{good i ≤ k} w_i.
*(i) Counting.* The bad set {i : dist(id, πℤ) < α} splits into maximal runs, one per multiple of π approached. Each run lies in an open interval of length
2α/d, so it has at most N_b := ⌊2α/d⌋+1 elements. Consecutive runs are separated by a closed interval of length (π−2α)/d, which contains
N_g ≥ ⌊(π−2α)/d⌋ ≥ (π−2α)/d − 1 good integers. Split [1, I] into consecutive (bad run, good gap) pairs. Every complete pair, and every final partial good gap,
has good fraction ≥ ρ := min N_g/(N_g+N_b). Only a final bad run (≤ N_b elements) can go uncounted. (The run at 0, {1 ≤ i < α/d}, has ≤ N_b elements, so it pairs with the gap that follows it like any other run. There is no extra loss at the start.) Hence **G(I) ≥ ρ(I − N_b)**.
N_g/(N_g+N_b) increases in N_g and decreases in N_b. At the extremes, N_g+N_b = π/d, so
ρ ≥ ((π−2α)/d − 1)/(π/d) = (π−2α−d)/π ≥ (π−2α−0.3)/π. **This confirms the constant claimed in LemmaA_uniform.md.**
(Note: a first pass of this audit wrote (π−2α−d)/(π+d) by bounding N_g+N_b separately. That was an unnecessary loss and is withdrawn.)
*(ii) Abel summation.* With w_{k+1} := 0: Σ_{good} w_i = Σ_{I≤k}(w_I − w_{I+1})G(I) ≥ ρΣ_{I≤k}(w_I − w_{I+1})(I−N_b)₊ = ρΣ_{N_b<i≤k} w_i.
Each coefficient w_I − w_{I+1} ≥ 0. w is decreasing, so w_i ≥ θ⁻¹∫_{iθ}^{(i+1)θ} w, and the sum is ≥ θ⁻¹∫_{(N_b+1)θ}^{kθ}.
The lower limit satisfies (N_b+1)θ ≤ 2θ + 2αθ/d ≤ 2θ_A + 4α/30 (d ≥ 15θ). The upper limit satisfies kθ = ν−θ ≥ ν_A − θ_A.
Goodness gives sin²(iφ/2) = sin²(id) ≥ sin²α, so each good summand of S is ≥ w_i.
*(iii)* This gives θS ≥ ρ∫_{2θ_A+4α/30}^{ν_A−θ_A} log(1+sin²α/sinh²(u/2))du. The integrand is decreasing, so a right Riemann sum in Arb is rigorous.
**Verdict: regime (M) is correct as stated.** The original script's lower limit ((2θ_A+4α/30)/ρ + θ_A) is larger than necessary, which only makes it more conservative.
ws_b_minorM.py (Arb, N = 4000) gives θS ≥ **1.2304** with the audited limit (α = 0.6) and θS ≥ 1.0580 with the original limit. This reproduces the 1.057 in ws_p12_minor.out.
Either value is far above the binding regime (a) (0.3185), so **c* = 0.15927 stands**. The hand-argued step is now fully written out above and checked.


### 7b. Resolution of the LemmaCD.md (B1) objection (2026-10-08)
**Claim (LemmaCD.md §2 (B1), line 80).** The LemmaA_uniform count of bad indices is missing a factor (1+α/π), and its lower integration limit drops a θ̄C₀ term.
**Their derivation.** The bad indices i ≤ N lie in windows around jπ with j ≤ (Nd+α)/π. That makes ≤ (Nd+α)/π + 1 windows, each holding ≤ 2α/d+1 integers. Hence
#bad ≤ N(2α+d)/π + (1+α/π)(2α/d+1), so #good ≥ ρN − C₀ with C₀ = (1+α/π)(2α/d+1). This is a **valid** argument. Its window count is cruder, but nothing in it is wrong.
**Why the original count is also valid (and sharper).** The run/gap pairing of §7(i) gives G(I) ≥ ρ(I − N_b) with N_b = ⌊2α/d⌋+1 ≤ 2α/d+1. Checking the end cases:
- Group [1, I] as (run, gap) units. The first unit is the run at 0 (i < α/d), which has ≤ N_b elements.
- A complete unit has good fraction ⌊A⌋/(⌊A⌋+⌊2α/d⌋+1) ≥ (A−1)d/π = (π−2α−d)/π, with A = (π−2α)/d.
- If I ends inside a run, that partial run has ≤ N_b elements, and they are the only ones not covered by the ρ bound.
- If I ends inside a gap, the final run (≤ N_b) is charged to the additive constant, and the partial gap is all good, so its fraction is 1 ≥ ρ.
So the additive loss is N_b, not (1+α/π)N_b. The factor 1+α/π comes from counting windows (including a fractional one) and multiplying by the maximum window size. The pairing does not need it.
**Lower limit.** The original limit θ(N_b+1)/ρ + θ ≤ (2θ_A + 4α/30)/ρ + θ_A already contains the θN_b term (θN_b ≤ 2αθ/d + θ ≤ 4α/30 + θ). Nothing is dropped.
In the Abel form of §7(ii) the limit is simply θ(N_b+1), which is smaller still.
**Verdict.** The original LemmaA_uniform regime (M) is correct. LemmaCD's version is a correct but more conservative variant, not a correction.
**Re-certification** (ws_b_minorM.py, Arb, all three variants):
- Abel limit: θS ≥ 1.2304.
- Original limit: θS ≥ 1.0580.
- **LemmaCD count** (C₀ = (1+α/π)(4α/30+θ_A), lower limit (θ_A + (1+α/π)(4α/30+θ_A))/ρ + θ_A: θS ≥ 1.0137.
Every variant is far above the binding regime (a) value 0.3185, so **c* = 0.15927 stands** whichever count is adopted.

## 8. Scripts and reproduction
- ws_b_core.py: M₂ as a function of the six base variables (identical algebra to ws_p12_eps.core), the forward-mode AD class D, and Lemmas B.1–B.4 as box functions.
- `python3 ws_b_step.py 0.002 0.0040939 -2.52 2.10 3.1158 1.02 0.02` → ws_b_step.out (Region I, Σ ≥ 0.91918 after the TV re-certification; previously 0.91905).
- `python3 ws_b_asym.py 0.002 -2.52 2.10 0.05` → ws_b_asym.out (Region II, Σ ≥ 0.96039, |log M₂|/θ ≤ 0.14215; previously 0.96036 / 0.14236, still valid and still used as input to ws_b_edges).
- `python3 ws_b_edges.py 100000 0.002 0.0040939 0.14236` → ws_b_edges.out (B3).
- `python3 ws_b_minorM.py 0.0042 3.0` → ws_b_minorM.out (audit).
- `python3 ws_b_check.py 100000 1210 1358 1500` → ws_b_check.out (sanity, non-essential).
Every run takes under 10 s.

## 9. Caveats (honest)
- Lemma B inherits Lemma W and the P1 TV constants from LemmaA_uniform.md. TV is now certified by ws_mech_tv.py and Lemma W by ws_mech_window.py / ws_mech_lemA.py (part G).
- [Resolved] The edge bounds ν'_L and ν'_R (§5) are certified for all β ≥ β_A by ws_mech_window.py; no monotonicity argument is used.
- In Region I, the convex hull box H is built from P1.3 enclosures at both endpoints. This is rigorous but loose, and Σ_min is driven by the ∂_νM₂·dν term at ν' ≈ 2.1.
