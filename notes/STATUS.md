# STATUS — computer-assisted proof of the Auluck–Chowla–Gupta conjecture (2026-10-08)

## Theorem
For every n ≥ 1, the sequence p(n,1), p(n,2), …, p(n,n) is weakly unimodal. Here p(n,k) is the number of partitions of n into exactly k parts.

Proof split:
- n ≤ 2·10⁵: exhaustive exact computation.
- n ≥ 10⁵: analytic argument with explicit constants, every number certified in Arb.

The two ranges overlap.

## Notation and region cover for n ≥ 10⁵
β = √(6n)/π, ℓ = log β, μ = βℓ, R_k = p(n,k+1)/p(n,k), L_k = log R_k. The regions are:

| Region | k range | Claim | Source |
|---|---|---|---|
| D | 1 ≤ k ≤ ⌊1.7n^{1/3}⌋ (≥ 78) | R_k > 1 | LemmaCD §4 (D1–D3) |
| CL | k_D < k ≤ μ−2β | R_k > 1 | LemmaCD CL-a (k ≥ 47, kθ ≤ 0.1) + CL-b (kθ ≥ 0.1) |
| W | μ−2β < k < μ+2β | sign L_k = sign L₂,k except ≤ 1 k; L₂ has one sign change | LemmaA_uniform (Thm A-unif + Corollary §0) + LemmaB (B1–B3) |
| CR | μ+2β ≤ k ≤ ⌈n/2⌉ | log R_k ≤ −0.614θ | LemmaCD C-right |
| E | ⌈n/2⌉ ≤ k ≤ n−1 | a_{k+1} ≤ a_k | LemmaCD §5 (p(n,k) = p(n−k), exact) |

Each R_k gets its sign independently. The regions only need to cover every k; they need not overlap.

Window logic:
1. Lemma A gives |L_k − L₂,k| ≤ ε_k ≤ Γ_max θq with Γ_max = 0.4282. The certified values are 0.42024 in Region I and 0.17605 in Region II.
2. Lemma B (B1) gives a step L₂,k − L₂,k+1 ≥ 0.9100·max(θq) > 2Γ_max·max(θq) = 0.8564·max(θq).
3. Lemma B (B2) gives |L₂,k| > Γ_max θq off the crossing pair.
4. Lemma B (B3) gives L₂/θ ≥ 5.588 at the left edge and ≤ −0.505 at the right edge.
5. Hence the sign pattern on W is +…+ ? −…−. Together with the other regions, the whole row has the pattern + … + (?) − … −, which is unimodal.

Ties are allowed: the claim is weak unimodality.

## Lemma chain and key constants

**Lemma W / CW (parameter ranges)**
- Constants: θ ≤ 0.0040939, ν′ ∈ [−2.5102, 2.0840], ν ≥ 3.1158, ν_CL ≤ −1.95361, ν_CR ≥ 1.04305.
- Certified for all β ≥ β_A with no monotonicity argument (ws_mech_window.py).
- The reduction to G > 0 is mechanized (ws_mech_lemA.py part G).

**P1.2 (TV of f_r)**
- Constants: TV ≤ 1.000001, 1.000001, 2.000001, 6.086765, 25.691943, 138.425.
- Certified by ws_mech_tv.py, using a valid variation bound on every cell.

**P2 (minor arcs)**
- Constant: c_* ≥ 0.15973. The downstream scripts use 0.15926.
- Proof: ws_mech_minor.py, by the Dirichlet identity. It needs no block or pair counting.

**Lemma A error algebra (LemmaA.md §§1–6)**
- The ε formula is machine-checked (ws_mech_lemA.py).
- The code (ws_p12_eps.core) equals an independent implementation of the proved formula to 1.7·10⁻¹³.

**Thm A-unif**
- Region I: Γ ≤ 0.42024 (ws_p12_boxes).
- Region II: Γ ≤ 0.17605 (ws_p12_asym).

**Lemma B**
- (B1): Σ ≥ 0.91918 (Region I, ws_b_step) and 0.96039 (Region II, ws_b_asym).
- (B3): ws_b_edges and ws_mech_edges, with ν′_L ≤ −1.906674.

**Lemma CD**

| Part | Value | Requirement | Script |
|---|---|---|---|
| CL-b margin | 5.8206 | > 1 | ws_cd_CLb, ws_mech_cd |
| CL-a ρ | 0.0673 (ws_mech_cla) / 0.161 | < 0.8977 | ws_mech_cla, ws_cd_CLa |
| CR | 0.3851 | < 1 | — |
| D, Ψ(s_A) | 0.30456 | > 0, Ψ increasing (exact Ψ′) | ws_cd_D, ws_close_D |

**Finite check**
- PASS for every n ≤ 200 000, with 0 violations.
- Two independent verifiers whose logs are byte-identical:
  - x87 interval arithmetic: fc/FiniteCheck.md, ws_fc_check.c.
  - exact integers: fc/FiniteCheckExact.md, ws_fc_exact.c.
- SHA-256 values are in fc/SHA256SUMS.

## Reproduction (in this folder; python-flint 0.9.0, sympy 1.14, gcc 13.3)
```
# analytic part, n >= 1e5 (total ~6 min)
python3 ws_mech_elem.py                 > ws_mech_elem.out
python3 ws_mech_tv.py                   > ws_mech_tv.out      # P1.2 TV constants
python3 ws_mech_lemA.py                 > ws_mech_lemA.out    # LemmaA §§3-6 algebra, Lemma W reduction
python3 ws_mech_window.py 100000 5      > ws_mech_window.out  # Lemma W / CW / B-edge constants, all beta
python3 ws_mech_minor.py 0.0042 3.0 0.4 > ws_mech_minor.out   # P2, c_* >= 0.15973
python3 ws_p12_window.py 100000 5       > ws_p12_window.out
python3 ws_p12_minor.py 0.0042 3.0      > ws_p12_minor.out
python3 ws_p12_boxes.py 0.002 0.0040939 -2.52 2.10 3.1158 0.15926 1.01 0.01 > ws_p12_regionI.out
python3 ws_p12_asym.py 0.002 -2.52 2.10 0.15926 0.05                    > ws_p12_regionII.out
python3 ws_b_step.py 0.002 0.0040939 -2.52 2.10 3.1158 1.02 0.02        > ws_b_step.out
python3 ws_b_asym.py 0.002 -2.52 2.10 0.05                              > ws_b_asym.out
python3 ws_b_edges.py 100000 0.002 0.0040939 0.14236                    > ws_b_edges.out
python3 ws_mech_edges.py > ws_mech_edges.out
python3 ws_cd_window.py 100000 5 > ws_cd_window.out
python3 ws_cd_CR.py 100000 > ws_cd_CR.out
python3 ws_cd_CLa.py 47 0.1 > ws_cd_CLa.out
python3 ws_cd_CLb.py 100000 0.1 78 > ws_cd_CLb.out
python3 ws_cd_D.py > ws_cd_D.out
python3 ws_mech_cd.py 100000 > ws_mech_cd.out
python3 ws_mech_cla.py 47 0.1 > ws_mech_cla.out
python3 ws_mech_mono.py > ws_mech_mono.out
# closure pass (PROOFS_FULL.md; total ~1 min)
python3 ws_close_rep.py > ws_close_rep.out; python3 ws_close_D.py > ws_close_D.out
python3 ws_close_elem.py > ws_close_elem.out; python3 ws_close_thresh.py expected ../notes/FiniteCheckExact.md > ws_close_thresh.out   # in the repo: reads expected/ (run_all.sh passes actual/)
# audits / sanity (non-essential)
python3 ws_b_minorM.py 0.0042 3.0; python3 ws_b_check.py 100000 1210 1358 1500; python3 ws_p12_eps.py 100000 1358
# finite part (in fc/; ~8 min + ~10 min)
gcc -O3 -march=native -Wall -o ws_fc_exact ws_fc_exact.c && ./ws_fc_exact 200000 exlog200k.txt && cmp exlog200k.txt log200k.txt
gcc -O2 -fopenmp -frounding-math -o ws_fc_check ws_fc_check.c -lm && ./ws_fc_check 200000 und200k.txt log200k.txt
```
ws_p12_tv.py is retired: its TV method was invalid (REFEREE A.1). lemA_*.py, scan.py and the other early scripts are exploratory and are not part of the proof.

## Trust list (updated 2026-10-08, after the PROOFS_FULL.md pass)
Items 1–8 of the previous list are closed. Each now has a complete written proof in **PROOFS_FULL.md** (no "clearly"; every inequality
justified), and each has the machine checks listed below. The identities are exact sympy or rigorous Arb; the numerical checks are sanity tests of
statements that are proved in writing.

| # | Item | Written proof | Machine check |
|---|---|---|---|
| 1 | Representation (1)/(2) | PROOFS_FULL §1 | ws_close_rep R1–R4: (1) exact with rational x on 4959 (n,k,x) triples; (2) by quadrature to 1e-20 |
| 2 | \|ψ^(r)\| ≤ κ_r, \|H^(r)\| ≤ c_r, Re H ≤ 0 | §2 (termwise series, M-test) | ws_close_rep R5 (random φ); identities at 0 in ws_mech_lemA B |
| 3 | [ELEM] one-liners | §3 | ws_close_elem X1/X2 (sympy remainder identities), X3 (numeric) |
| 4 | P1.1 and closed form of F_r | §4 | ws_close_elem F1/F2 (sympy telescoping + limits, r ≤ 8), F3 (Arb rigorous quadrature) |
| 5 | Lemma B.1–B.4, T₂, hull MVT, Region II majorants | §5 | ws_close_elem B1/B4/B5 (sympy), B6 (B.1–B.4 + hull at real saddles, n = 1e5 and 1e6) |
| 6 | Lemma C0/C1 error terms, δ_T/δ_C, CL-a N₁/N₂/N₃₄ | §6 | ws_close_elem C0a–d, C1a–b, CLa (sympy) |
| 7 | Lemma D1–D3, **Ψ(s) increasing** | §7 | ws_close_D: **Ψ′ in exact closed form (> 0 for s > √γ)**, Arb Ψ(s_A) = 0.30456, Arb Ψ′-sweep, exact a_{k+1} > a_k at n = 1e5, 2e5 |
| 8 | Region E | §8 | ws_close_rep E1/E2 (exact, n ≤ 400) |

**Threshold audit** (PROOFS_FULL §9, ws_close_thresh.py). Every lemma's hypotheses hold for **all** n ≥ 10⁵, and D ∪ CL ∪ W ∪ CR ∪ E = [1, n−1]:
- k_D ≥ 78;
- μ + 2β ≤ n/2;
- CL, W and CR lie in 2 ≤ k ≤ n−2;
- θ_CR is decreasing and v ≥ 6.198;
- the box covers ⊇ the certified ranges.

No genuine gap was found. Side observation: (B3) is not load-bearing for unimodality (§9.3).

**Remaining trust base (only):**
- python-flint/Arb 0.9.0 (ball arithmetic, polylog, zeta, erfc, acb.integral);
- sympy 1.14;
- gcc 13.3 and the two independent C finite-check programs (the exact-integer one uses no floating point);
- the classical theorems quoted in PROOFS_FULL §11: Taylor's theorem with integral remainder, the mean value theorem, the M-test,
  additivity of total variation, and the identity theorem.

mpmath appears only in sanity tests.
