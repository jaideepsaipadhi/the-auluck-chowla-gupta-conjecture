# Hostile referee report: computer-assisted proof of the Auluck–Chowla–Gupta conjecture

Date: 2026-10-08. Scope: ACG_attack.md, INTERFACE.md, LemmaA.md, LemmaA_uniform.md, LemmaB.md (including §7/§7b), LemmaCD.md,
fc/FiniteCheck.md, fc/FiniteCheckExact.md, and all scripts and logs. Re-runs are in ../rerun/.

**Verdict.** I found no FATAL error. There is one SERIOUS rigour gap. It sits in a constant that every analytic file uses
(the TV bounds), but the true values are far below the hard-coded ones, so the gap can be repaired and does not threaten
the result. The rest are MINOR issues: documentation, and monotonicity claims that are argued by hand. The logical assembly
closes, and every script reproduces its recorded output.

---

## A. Issues

### 1. SERIOUS: the TV certificate (Lemma P1.2, ws_p12_tv.py) is not rigorous as written
*Location:* LemmaA_uniform.md §2, Lemma P1.2; ws_p12_tv.py, the `else` branch; constants hard-coded in ws_p12_eps.py:12.

*Defect:* the Arb enclosure of f_r′ may fail to have a constant sign on a cell. The script then charges that cell
g(b)^r h_r(a) − g(a)^r h_r(b). That is (sup f − inf f) on the cell. It is the variation **only if f is monotone on the cell**.
When f really turns on the cell, the variation can reach 2·(range), and the correct product bound is
Var(GH) ≤ sup H·Var G + sup G·Var H = H(a)(G(b)−G(a)) + G(b)(H(a)−H(b)).
The enclosure fails to have a sign on many cells (64, 365, 1042, 2323, 2148 and 2631 cells for r = 1..6). These failures are
mostly overestimation in the ball arithmetic, but the script cannot tell which failures are real turning points.

*Effect:* I recomputed with the correct product bound (../rerun/tvprod.py) and got
TV ≤ 1.06268, 1.70560, 7.71293, 51.5566, 236.267, 1628.70.
**Every one of these exceeds the hard-coded value** (1.062650, 1.705253, 7.708993, 51.518986, 236.047594, 1626.897074).
So the constants actually used are **not certified** by the stated method.

*Why it is not fatal:* the true variations are far smaller. Numerically (mpmath, fine grid) TV(f₁..f₆) ≈ 1.00, 1.00, 2.00,
6.09, 25.7, 138. The hard-coded values overestimate by a factor of 1 to 12, so certifying them correctly (finer cells on the
cells without a sign, or the product bound plus a small bump in the constants) will not reduce the margins.
It still has to be done. Afterwards, re-run ws_p12_boxes, ws_p12_asym, ws_b_step, ws_b_asym and the ws_cd_* scripts,
since all of them import TV.

### 2. MINOR: INTERFACE.md and ACG_attack.md still define the old main term L̃
INTERFACE.md "Objects" and ACG_attack.md Step 3 give L̃ (h′, h″, κ₃ only) as "the explicit main term". The proof uses
L₂ = −θ − log(1−q) + log M₂ (LemmaA.md §3) throughout. I checked that L₂ is used consistently:
- ws_b_core.M2_formula agrees with ws_p12_eps.core to 1e-13 (ws_b_check.out, reproduced).
- The symbolic gap M_model − M₂ = Num/N₀ matches a direct evaluation to 6 significant figures at (10⁵, 900), (10⁵, 1358),
  (10⁵, 1850) and (3·10⁵, 2586).

The final-logic paragraph of INTERFACE.md should be updated to name L₂. As written it says something that is not what was proved.

### 3. MINOR: monotonicity in β is argued by hand only
Locations:
- Lemma W (θ_A, ν′_hi, ν′_lo, ν_min, G(β));
- LemmaB §5 (ν′_L, ν′_R);
- LemmaCD Lemma CW (ν_CL, ν_CR, θ_CR);
- Lemma D3 (Ψ(s) increasing).

These claims are what extend the n = 10⁵ values to all n ≥ 10⁵. The hand arguments are correct as far as I can follow them.
I also checked numerically (mpmath, β from β_A to about 10²⁰, geometric ratio 1.05) that each quantity is strictly monotone in
the claimed direction. Ψ(s) is increasing, with minimum 0.3042 at s = n^{1/3} for n = 10⁵. A short written derivative proof
for each would remove this caveat.

### 4. MINOR: Lemma P2 regime (M) has two versions of the count
LemmaCD.md §2 (B1) says its count "corrects" LemmaA_uniform. LemmaB.md §7b replies that both counts are valid.
I agree with §7b: the run/gap pairing gives G(I) ≥ ρ(I − N_b), and Abel summation then gives the lower limit (N_b+1)θ.
All three variants certify θS ≥ 1.01, far above the binding regime (a) value of 0.3185.
The wording in LemmaCD §2, line 80 should be changed so the documents do not contradict each other.

### 5. MINOR: the documented ν′ range does not match the command
LemmaA_uniform.md §7 runs ws_p12_boxes with ν′_hi = 2.09. The cell count `int((VHI−VLO)/h)+1` extends the cover to 2.10,
so ν′_hi = 2.0840 (Lemma W) is covered: the last output cell is [2.09, 2.10]. This works by a side effect of the formula.
Pass 2.10 explicitly.

### 6. MINOR: wording in Theorem B (B2)
(B2) leaves both k* and k*+1 unseparated. "At most one undetermined k" therefore comes from (B1), not from (B2):
if both |L₂,k*| ≤ ε_k* and |L₂,k*+1| ≤ ε_{k*+1}, then the step is ≤ ε_k* + ε_{k*+1} ≤ 2Γ max(θq) < 0.9099 max(θq).
That contradicts (B1). The argument is correct, but this step should be written out in LemmaA_uniform §0, where the
Corollary is stated without proof.

### 7. MINOR: notes on the pointwise LemmaA.md
The minor-arc supremum is taken in float64, and the table shows the bound fails at n = 10⁴. Neither matters, because the
uniform version (LemmaA_uniform) replaces this file and uses only Arb plus Lemma P2. LemmaA.md should carry a header line
saying it is superseded and not part of the proof.

### 8. MINOR: conservative slack (no defect)
- CL-b uses θ̄ = v₂/78, but k ≥ 79.
- Lemma C1 uses b ≥ F₂(v₁) − θ, which does not need TV.
- The (M) lower limit is larger than necessary.

All of these are sound and conservative.

---

## B. Verified OK

### Logical assembly (n ≥ 10⁵)
The regions chain with no gaps:
- D: 1 ≤ k ≤ ⌊1.7n^{1/3}⌋, at least 78.
- CL: k_D < k ≤ μ−2β. CL-a needs k ≥ 47; CL-b covers v ≥ 0.1.
- W: μ−2β < k < μ+2β, integers.
- CR: μ+2β ≤ k ≤ ⌈n/2⌉.
- E: k ≥ ⌈n/2⌉.

Each R_k gets its sign on its own, so boundary k need not lie in two regions. A two-sided check of the window ends:
- First k of W: K ≤ μ−2β+2. Last k of W: K ≥ μ+2β.
- Actual values at n = 10⁵: ν′ = −2.198 at k = 865 and +2.067 at k = 1851.
- At n = 2·10⁵ the values are −2.156 and +2.049. At n = 10⁶ they are −2.089 and +2.027.
- All lie inside [−2.5102, 2.0840].

The gluing works as follows:
- In W, |L − L₂| ≤ Γθq < |L₂| except possibly at one k, so the determined signs equal sign L₂.
- L₂ is strictly decreasing on W, which gives +…+ ? −…−.
- CL gives + on the left and CR/E give − on the right.
- The result is unimodal, ties included.

The finite range n ≤ 2·10⁵ overlaps the analytic range n ≥ 10⁵.

### Constants carried between files
- Γ_max = 0.42815, so 2Γ = 0.8563.
- Σ_min = 0.91905. Using θ′q′ ≤ 1.01θq gives step ≥ 0.9099·max(θq) > 0.8563. ws_b_step uses the requirement 2Γ(1+τ) = 0.8650 < 0.919.
- c_* = 0.15927 in ws_p12_minor and 0.15926 passed to the box scripts (conservative).
- θ_A = 0.0040939 and ν_min = 3.1158 are consistent in every file.
- P2 requires θ ≤ 0.0042 and ν ≥ 3.0. Both hold.
- The θ ≤ 0.002 / θ ≥ 0.002 split is the same in Lemma A and Lemma B.

### Mathematics re-derived
- Representation (1)/(2): R = Σ_s x^{1+sK} P(X = m−1−sK)/P(X = m). G(φ) = xe^{iφ}/(1 − qe^{iKφ}). R/g₀ − 1 = Ñ/D.
- |ψ^{(r)}| ≤ κ_r for all real φ.
- |H′| ≤ c₁, |H^{(r)}| ≤ c_r, and Re H ≤ 0 (because |G| ≤ g₀).
- T_ω, built from:
  - the omitted ω̂²/2 and ω̂³/6 terms, using 3a²r + 3ar² + r³ ≤ 3r(a+r)²;
  - the quartic remainder;
  - the factor e^{δ₀Bt²/2}, with δ₀ including the κ₆ term.
- T_H = c₅t⁵/120 + (c₁c₂/2 + c₁³/6)t³.
- Gaussian tails; the √2 normalisation; the ratio-perturbation formula; and ε = −log(1 − Δ/M₂).
- Lemma C0 (the T₁, T₂ and A₃ terms with their (1−δ) powers) and δ_T.
- Lemma W, including G and T(v) ≤ (v+1)e^{−v}/(1−e^{−v}).
- Lemmas B.1–B.4, including u f_r′ = r f_r − f_{r+1}.
- The T₂ lower bound and the convex-hull mean-value-theorem step.
- The monotone-majorant criterion ν_min > 2j/P, which the scripts assert for every monomial.
- The θ-monotonicity of the tail and minor-arc factors.
- Region II δ₀ is free of θ.
- P2 regimes (a), (M) and (b): Riemann-sum direction, the cell-endpoint choices and the pair argument.
- CL-a pieces N₁–N₄, including √B ≤ k^{3/2}/v.
- CR.
- Lemma D1–D3: the distinct-parts injection and the product bound.
- Region E.

### Scripts re-run (all under 1 min except the finite check)
- Identical output: ws_p12_regionII, ws_p12_minor, ws_p12_window, ws_b_asym, ws_b_edges, ws_b_minorM, ws_b_check,
  ws_cd_window, ws_cd_CR, ws_cd_CLa, ws_cd_CLb and ws_cd_D.
- ws_p12_regionI and ws_b_step differ only in the timing column.
- ws_p12_tv reproduces its constants. Note issue 1: the method behind them is flawed.
- ws_p12_sym and ws_p12_eps (10⁵, 1358) agree with lemA_bound to 6 digits.

### Finite check
- I rebuilt ws_fc_check and ran it for N = 10⁵ in 28 s. It found 0 undetermined pairs, and its SHA-256 matches log100k.txt.
- I rebuilt the **exact-integer** ws_fc_exact and ran it for N = 2·10⁵ (about 8 min). Result: PASS, 0 violations, and its
  SHA-256 equals log200k.txt (45ebb662…).
- log2k.txt is byte-identical to ref2k.txt.
- The column recurrence and the pair range (k−1 ≤ ⌈n/2⌉) are correct.
- The independent exact verifier makes the x87 rounding caveat irrelevant.

### Numerical spot-check against exact p(n,k) at n = 10⁵
I used the scan.py-style column recurrence and Arb M₂. Over k = 865 … 1851 (both ends of W and the crossing 1357–1359),
|L − L₂|/(θq) ≤ 0.0004. That is about 1000 times below Γ_max = 0.428, and well below the rigorous ε.
L₂/(θq) = +0.574 at k = 1358 and −0.480 at k = 1359, so both exceed Γ_max: no k is undetermined at n = 10⁵.
The exact mode is 1359, which agrees with the finite-check log.

## C. What must be done before this counts as a proof
1. Re-certify TV₁…TV₆ with a correct per-cell variation bound (issue 1), then re-run all dependent scripts.
2. Update INTERFACE.md and ACG_attack.md to state L₂ as the main term, and mark LemmaA.md as superseded (issues 2 and 7).
3. Optionally, write derivative proofs for the β-monotonicity claims (issue 3).
