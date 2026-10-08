# MECH: removing the hand-argued steps of the ACG proof (n ≥ 10⁵)

Status 2026-10-08 (updated: §3a audit + adoption, §6 new). Every new script runs in under 3 minutes. The longest is ws_mech_cd.py at 168 s; the others take under 45 s.
Interval arithmetic is Arb (python-flint 0.9.0). Exact identities are checked in sympy 1.14. The outputs are in ws_mech_*.out.

## 0. Shared elementary facts: ws_mech_elem.py → ws_mech_elem.out
Every one-line calculus fact used below is checked here, either as an exact sympy identity or as a sign that is a product of manifestly signed factors:
E1 u/(eᵘ−1)↓; E2 x/sinh x↓; E3 sin x/x↓ on (0,π); E4 log(1+Ax) ≥ x log(1+A) on [0,1] (concavity);
E5/E9 the Dirichlet identity Σ_{i≤I} sin²(id) = I/2 − sin(Id)cos((I+1)d)/(2 sin d) (telescoping + sum-to-product, symbolic; plus 500 Arb spot checks);
E6 lᵃe^{−l}↓ for l > a; E7 d/dt[t^P ν^j] = t^{P−1}ν^{j−1}(Pν−2j) with ν = ν'+2log(1/t); E8 1+csch² = coth², log coth↓;
E10 sin x ≥ x − x³/6; E11; E12 d/dθ log(θ^{−p}e^{−a/θ}) = (a−pθ)/θ²; E14 u f_r' = r f_r − f_{r+1} (r = 1..6, exact);
E15 h' for h = u/(eᵘ−1); E16 the Abel summation identity (k ≤ 10, exact); E17 ∫_v^∞ f₁ = −v log(1−e^{−v}) + Li₂(e^{−v}).

## 1. Minor arcs (LemmaA_uniform P2 regimes (M)/(b), LemmaCD C-minor (B1)/(B2), CL-a N₃/N₄)

**Old status.** Regime (M) used block counting plus a density claim plus Abel summation, and regime (b) used a pair argument. Both were hand proofs, audited in LemmaB §7/§7b.
LemmaCD (B1)/(B2) and CL-a's N₃/N₄ used the same kind of argument.

**New proof (ws_mech_minor.py).** This replaces all counting with one identity.
- **(a)** For c = φ/θ ∈ [0.4, c_a], the proof is unchanged in substance. It uses E2 and E3 and is evaluated on 400 geometric c-cells with 41 values of T.
- **(D)** For φ ∈ [c_aθ, π], put d = φ/2 ∈ (0, π/2] and w_i = 2 log coth(iθ/2) (decreasing, E8).
  1. E4 gives each summand ≥ sin²(id)·w_i.
  2. Abel summation (E16) applies, with partial sums G_I = Σ_{i≤I} sin²(id) ≥ (I − ⌈1/sin d⌉)₊/2 (E9).
  3. Together these give S ≥ ½Σ_{i>N₀} w_i ≥ θ⁻¹∫_{(N₀+1)θ}^{(k+1)θ} log coth(u/2) du.
  4. sin d ≥ (c_aθ/2)·sinc(c_aθ_max/2) by E3. This gives the θ-free lower limit L = 2θ_max + 2/(c_a sinc(c_aθ_max/2)).

No density, block or pair argument remains. There is also no θ-box cover and no θ → 0 asymptotics: θ enters only through θ ≤ θ_max and (k+1)θ ≥ u_top, so the bound holds uniformly as θ → 0.
Every integral is a right Riemann sum of a decreasing integrand, evaluated in Arb.

**Output** (`python3 ws_mech_minor.py 0.0042 3.0 0.4`):
- (a) gives θS ≥ 0.31946 on c ∈ [0.4, 1.331].
- (D) gives θS ≥ 0.34433 on φ ∈ [1.331θ, π].
- **c_* ≥ 0.15973.** The old value was 0.15927, and every downstream script uses 0.15926, so the input is still valid. The constant has improved slightly.

**LemmaCD (ws_mech_cd.py → .out).** Lemma C-minor's Ψ is replaced by the mechanical bound with u_top = v₁ (since (k+1)θ ≥ kθ ≥ v₁). CR and CL-b are then re-run with the project's own Lemma C0/C1 code.
- CR: Ψ = 0.33808 (old value 0.32647). e^{−ν'}/(1−q) + X ≤ 0.385274 < 1, unchanged.
- CL-b: every one of the 34 cells passes. The minimum margin is **5.8165**, unchanged; the binding cell [3.29, 3.79] does not depend on the minor-arc bound. Small-v cells move slightly, for example 13.16 → 12.23 at v ∈ [0.1, 0.115], and all remain ≥ 11.

**CL-a (ws_mech_cla.py → .out).** N₃ + N₄ (block counting + pairs) is replaced by one term N₃₄.
- For φ ∈ [π/k, π], E4 and E9 give Σ sin²(id) ≥ k/2 − 1/(2 sin(π/2k)).
- Hence |χ| ≤ (vσ₀/2)^{2gk} with g = (½ − 1/(π sinc(π/2K₀)))/2.
- N₃₄ is increasing in v and decreasing in k. The monotonicity conditions are asserted.
- **Result: N₃₄ = 6.3·10⁻⁸** (old N₃ + N₄ = 0.0388). ρ ≤ 0.0673 (old 0.1609), against the requirement < 0.8977.

## 2. Lemma B edge bounds: monotonicity in β
**Old status.** ν'_L ≤ −1.9067 and ν'_R ≥ 1.0431 were computed at β_A and extended to all β by a hand monotonicity argument. The same pattern was used for θ_A, ν'_hi/lo, ν_min, ν_CL and ν_CR in Lemma W and Lemma CW.

**New proof (ws_mech_window.py).** No monotonicity is needed.
- Write x = 1/β, ℓ = log β, and express each bound in the atoms x, ℓx, ℓ²x, 1/ℓ. (r−1) is computed stably as eps·ℓ/(s(1+s)) with s = √(1−eps).
- Run an adaptive Arb bisection over ℓ ∈ [ℓ_A, 60]: 5486 boxes, each claim certified on every box.
- Close [60, ∞) with one tail box: x ∈ [0, e^{−60}], ℓx ∈ [0, 60e^{−60}], ℓ²x ∈ [0, 3600e^{−60}] (E6), 1/ℓ ∈ [0, 1/60].

**Output.** Every claim is certified for all n ≥ 10⁵. Worst box bound and stated constant:

| Quantity | Worst box bound | Stated constant |
|---|---|---|
| G | ≥ 1.0092 | > 0 |
| θ_A | ≤ 0.00409389 | 0.0040939 |
| ν'_hi | ≤ 2.08400 | 2.0840 |
| ν'_lo | ≥ −2.510198 | −2.5102 |
| ν_min | ≥ 3.115818 | 3.1158 |
| ν'_L (sharp form) | ≤ −1.94916 | −1.9067 |
| ν'_L (printed majorant) | ≤ −1.9066744 | −1.9067 |
| ν_CL | ≤ −1.9536125 | −1.95361 |
| ν_CR | ≥ 1.0430522 | 1.04305 |

**Constant change.** The printed bound ν'_L ≤ −1.9067 is only an Arb-rounded upper bound of its value at β_A, which is −1.906674922. Taken over all β ≥ β_A, the supremum of the *printed majorant* is −1.9066744. That is above −1.9067 by 3·10⁻⁵, so the stated digits are marginally wrong.
- The sharp expression −2 + (r−1)(ℓ−2) + 2r/β + log r is ≤ −1.94916 for all β, which is far below the requirement.
- (B3) was re-run with both −1.906674 and −1.94915 (ws_mech_edges.py → .out). Left: L₂/θ ≥ 5.588 (Region II) and 5.621 (Region I); with the sharp value, 5.880 and 5.911. Right: ≤ −0.505 and ≤ −0.643.
- **Fix needed in LemmaB.md §5:** replace "−1.9067" by "−1.906674", or by the sharp −1.9491. No conclusion changes.

The other constants hold: the ν_min, ν_CL and ν_CR values printed in the docs are below/above the certified worst-box bounds in the needed direction.
Independent sanity check (non-rigorous): actual saddles at both window ends for n = 10⁵, 10⁶ and 10⁷ lie in [(1−δ)/β, r/β]. At n = 10⁵ they give ν' = −2.198 and +2.067.

## 3. Independent recomputation of the P1 constants
**(3a) TV_r (ws_mech_tv.py → .out).** This was written from scratch and imports nothing from the project.
- On u ≤ 2.5, f_r and f_r' come from the Bernoulli series to n = 160, with a rigorous tail bound (|B_n|/n! ≤ 4/(2π)ⁿ).
- On u > 2.5, f_r uses Arb's polylog, and f_r' = (r f_r − f_{r+1})/u (E14).
- The variation is |Δf| on cells where Arb certifies a constant sign of f', with adaptive bisection elsewhere. The tail [60, ∞) is monotone by E6.

| r | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| new TV_r ≤ | 1.000000 | 1.000000 | 2.000000 | 6.086765 | 25.69195 | 138.4250 |
| project | 1.06265 | 1.70526 | 7.70900 | 51.5190 | 236.048 | 1626.90 |

TV₁ = TV₂ = 1 exactly, which matches f₁ and f₂ being monotone with sup 1.

**Audit of ws_mech_tv.py for the REFEREE A.1 flaw (2026-10-08): sound.** On a cell, the charge is |f(b)−f(a)| **only** when Arb certifies
f_r′ > 0 or f_r′ < 0 on the whole cell (then f is monotone there and |Δf| is the variation). Otherwise the cell is bisected, and at depth 14 the
charge is width·sup|f_r′| ≥ ∫|f_r′| = Var — valid whether or not f turns. No sup−inf charge is used anywhere. Other checks: the f′ enclosure is
taken on the ball covering the whole cell (`a.union(b)`); the series branch is used only when the ball's upper end is ≤ 2.5 < 2π (radius of
convergence), with a tail bound 2·(first omitted term) justified by the asserted term ratio < ½ and |B_n|/n! = 2ζ(n)/(2π)ⁿ ≤ 4/(2π)ⁿ and
|(n−1)…(n−r+1)| ≤ n^{r−1}; the n = 0 coefficient gives f_r(0) = (r−1)!; the u > 2.5 branch uses Arb polylog and the exact identity E14;
[60, ∞) is monotone by E6 with f_r → 0. The two branches agree with Arb polylog at 200 random points. The TV used in P1.1 is over [0, kθ] ⊂ [0, ∞).

**Adopted.** ws_p12_eps.TV now holds 1.000001, 1.000001, 2.000001, 6.086765, 25.691943, 138.425 (the printed Arb upper ends rounded up;
note 25.69 would be too small — the enclosure is 25.6919422). ws_p12_tv.py is marked RETIRED. All importers were re-run (old outputs in old_out/):

| Script | Old | New (certified TV) |
|---|---|---|
| ws_p12_boxes (Region I, now with ν′_hi = 2.10 explicit) | Γ ≤ 0.42815 | Γ ≤ 0.42024 |
| ws_p12_asym (Region II) | Γ ≤ 0.17782 | Γ ≤ 0.17605 |
| ws_b_step (Region I) | Σ ≥ 0.91905 | Σ ≥ 0.91918 (need 0.8650: OK) |
| ws_b_asym (Region II) | Σ ≥ 0.96036, |log M₂|/θ ≤ 0.14236 | Σ ≥ 0.96039, ≤ 0.14215 |
| ws_b_edges / ws_mech_edges | left 5.588/5.621, right −0.505/−0.643 | unchanged to 4 digits |
| ws_cd_CR | 0.385275 | 0.385132 |
| ws_cd_CLb, ws_mech_cd | margin 5.8165 | 5.8206 |
| ws_cd_CLa, ws_cd_D, ws_cd_window, ws_p12_window, ws_p12_minor, ws_b_check | — | byte-identical |
| ws_mech_mono | asserts pass | asserts pass (quotes the new Γ/Σ) |

The declared constants Γ_max = 0.4282 and 2Γ_max = 0.857 are kept; they now hold with slack (0.42024).

**(3b) Lemma W.** The reduction chain was re-derived from the LemmaA_uniform text.
- The T(v) tail bound is checked by Arb spot checks; the proof is that Li₁(y) and Li₂(y) are each ≤ y/(1−y) termwise.
- G > 0 and the θ-range constants are now certified on all β ≥ β_A, as in §2, not just at β_A. The worst G is 1.0092, attained near ℓ_A, against the old pointwise 1.0459.
- **Caveat:** the *algebraic* reduction (S(θ₋) > m ⇐ G > 0) is the project's derivation, re-read and re-encoded independently but not machine-derived.

## 4. "Monotone majorant in θ" claims (ws_mech_mono.py → .out)
- **LemmaCD C1.** The θ-derivatives of T₁/θ, T₂/θ, A₃ and δ_T, and the negative θ-derivative of w₀², are computed symbolically in sympy. Each is a polynomial with all coefficients ≥ 0 in positive symbols, divided by a positive denominator.
  M and M/θ are increasing for θ < Ψ/5 (E12), which X_cell asserts. Hence ρ/θ ≤ its value at θ̄. **Verified.**
- **LemmaB Region II.** c₁θ² is increasing because d/dθ[θ²ν] = θ(2ν−1) and ν > 3 (E13 identity re-checked inline). The (1+τ)-step majorants follow from Lemma B.1, which is asserted.
- **LemmaA Region II and LemmaB Region II (graded monomials).** The criterion "t^Pν^j nondecreasing iff Pν ≥ 2j" is now the symbolic identity E7. Both scripts are re-executed with live assertions on every monomial: Γ ≤ 0.17782 and Σ ≥ 0.96036, both unchanged. The tail factors e^{−a/θ}θ^{−p} are covered by E12 together with the scripts' asserts.
- **LemmaCD cells (CL-b, CR).** These rely only on C1, which is verified above.
- **LemmaB asym.** Same as the graded-monomial item; verified.

## 5. Downstream margins (with the certified TV, 2026-10-08)

| Requirement | Value | Status |
|---|---|---|
| ε ≤ 0.4282·θq (Γ_max) | 0.42024 (Region I) / 0.17605 (Region II) | Holds with more slack (was 0.42815 / 0.17782) |
| step ≥ 0.857·max θq | Σ_min = 0.91918 → 0.9100·max | Holds (was 0.91905) |
| c_* ≥ requirement (≥ 3.5θ and the 0.15926 used) | 0.15973 | Holds |
| CL-b margin > 1 | 5.8206 | Holds (was 5.8165) |
| CR < 1 | 0.3851 | Holds |
| CL-a ρ < 0.8977 | 0.0673 (mech) / 0.161 (orig) | Holds |
| (B3) left > 0 / right < 0 | 5.588 / −0.505 | Holds, with ν'_L = −1.906674 |

## 6. LemmaA.md §§3–6 error algebra and the Lemma W reduction (ws_mech_lemA.py → .out, 5 s)
sympy exact checks (all [OK]):
- **A** Gaussian moments (p ≤ 20, by sympy integration), the √(B/B′) and √2 normalisation factors, the tail-shift inequality.
- **B** ψ-factor and H derivatives at 0 (r ≤ 6 / ≤ 5) equal i^rκ_r and (ic₁, −c₂, −ic₃, c₄, ic₅) with the Eulerian closed forms of Li_{1−r};
  |1−qe^{ix}|² − (1−q)² = 2q(1−cos x), giving |G| ≤ g₀ and Re H ≤ 0.
- **C** D₀ = N₀ exactly; Ñ₀ real; doc-M₂ ≡ code-M₂; **M_model − M₂ = Num/N₀ with Num identical to the 11-monomial polynomial hard-coded in
  ws_p12_eps.core**; every piece is invariant under the scaling κ_r→λ^rκ_r, c_r→λ^rc_r, B→λ²B (justifies the B = 1 scaled core).
- **D** T_ω: the polynomial part of e^ω − Q (cubic Taylor of e^{ω̂} minus Q) is dominated coefficient-wise by the stated terms; the integral
  Taylor remainder identity for e^z; Re ω̂ = κ₄t⁴/24; δ₀Bt²/2 − (κ₄t⁴/24 + κ₆t⁶/720) factors as t²[…] ≥ 0 on |t| ≤ φ₀.
- **E** e^H − P_H = (e^H−1−H−H²/2) + (H−ΣH_r) + (H−H₁)(H+H₁)/2 exactly; T_H is the sum of the piece bounds; Q_abs = |Q| coefficient-wise;
  PH1abs dominates |P_H − 1|.
- **F** integrand identity χ(e^H−1) − gQ(P_H−1) = g[(e^ω−Q)(e^H−1) + Q(e^H−P_H)]; ratio-perturbation identity; the log bound.
  **Independent Arb implementation of ε from the verified sympy polynomials agrees with ws_p12_eps.core to 1.7·10⁻¹³ (relative) on 25 random inputs**,
  so the code that produces Γ computes exactly the formula proved.
- **G** Lemma W: (2δ−δ²)ζ₂ identity; T(v) = Σ e^{−lv}(v/l + 1/l²) (sympy integral) ≤ (v+1)e^{−v}/(1−e^{−v}); e^{−v} in atoms; v+1 = (ℓ−1) − (ℓ−2)δ;
  division by ℓx yields exactly G. G > 0 on all β ≥ β_A is ws_mech_window.py.
Lines tagged [ELEM, by hand, one line] are elementary steps stated but not machine-derived: Re ω ≤ κ₄t⁴/24 + κ₆t⁶/720 (Taylor with |ψ^{(6)}| ≤ κ₆
and odd orders imaginary), the three |·| bounds in T_H (|H| ≤ c₁t, |e^z−1−z−z²/2| ≤ |z|³/6 for Re z ≤ 0, Taylor of H), L − L₂ = log(R/(g₀M₂)),
the minor-arc length bound, and t²S(t) ≥ F₁(Kt) − t (left Riemann sum of decreasing f₁ ≤ 1).

## 7. What remains not machine-checked
- The [ELEM] one-liners of §6 and the "|ψ^{(r)}| ≤ κ_r, |H^{(r)}| ≤ c_r for all real φ" termwise-series argument (LemmaA.md §2).
- The exact representation (1)/(2) of LemmaA.md §1 (probabilistic/generating-function identity) — checked numerically by the referee, not symbolically.
- Lemma P1.1 (Riemann sum ≤ θ·TV) and the closed form of F_r (checked against quadrature only).
- Lemma B.1–B.4 and the convex-hull mean-value step of Lemma B are implemented in ws_b_core but their derivations are by hand (referee-audited).
- Lemma C0 error terms and Lemma D1–D3 are hand derivations (referee-audited); their monotonicity-in-θ claims are mechanized (§4).
- The finite check relies on two compiled C programs (independent; exact-integer one has no floating point).
- Software trust base: python-flint/Arb 0.9.0, sympy 1.14, gcc 13.3.

Reproduce: `python3 ws_mech_elem.py; python3 ws_mech_tv.py; python3 ws_mech_window.py 100000 5; python3 ws_mech_minor.py 0.0042 3.0 0.4;
python3 ws_mech_cd.py 100000; python3 ws_mech_cla.py 47 0.1; python3 ws_mech_edges.py; python3 ws_mech_mono.py; python3 ws_mech_lemA.py`
