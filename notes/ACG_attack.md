# Auluck–Chowla–Gupta conjecture (1942) — attack notes

**Conjecture.** For every n ≥ 1, the row p(n,1), …, p(n,n) is (weakly) unimodal,
where p(n,k) = #partitions of n into exactly k parts.
Posed by Auluck–Chowla–Gupta (J. Indian Math. Soc. 6, 1942); Erdős (Bull. AMS 1946):
"I am unable to prove or disprove this conjecture." Szekeres (1953) proved it for
n > N₀ with N₀ not effective. No proof for all n known (Canfield 1997; Canfield–Corteel–Savage 1998).

## Notation
f_k(m) = p_{≤k}(m) = [x^m] F_k(x),  F_k(x) = ∏_{j≤k} (1−x^j)^{−1}.
a_k = a_k(n) = p(n,k) = f_k(n−k).   m := n−k.
β = √(6n)/π,  μ = β log β,  x = (k−μ)/β   (Erdős–Lehner Gumbel scaling).

## Step 1 — exact reduction to a single generating function
Since F_{k+1} = F_k · (1−x^{k+1})^{−1},
  a_{k+1} = f_{k+1}(n−k−1) = Σ_{s≥1} f_k(n − s(k+1)),
so, with g(x) := x/(1 − x^{k+1}),
  **R_k := a_{k+1}/a_k = [x^m](g·F_k) / [x^m] F_k.**
Numerator and denominator are coefficient extractions of the SAME F_k at the SAME m,
so all exponential main terms cancel exactly; only the smooth multiplier g remains.

## Step 2 — unimodality from first differences only (no log-concavity needed)
Let L_k := log R_k. Suppose an explicit L̃_k satisfies |L_k − L̃_k| ≤ ε_k, and
  (a) on the crossing window W, L̃ is strictly decreasing with steps ≥ η_min, and 2ε < η_min;
  (b) outside W, |L̃_k| > ε_k.
Then sign(L_k) is determined for every k except at most ONE k, so the sign pattern of
(a_{k+1} − a_k) is (+…+, ?, −…−): unimodal. Region k ≥ n/2 is exact: a_k = p(n−k), nonincreasing.

> **Note (2026-10-08).** The proof does NOT use L̃ below. It uses the second-order main term L₂ = −θ − log(1−q) + log M₂ of LemmaA.md §3
> (see INTERFACE.md and STATUS.md). L̃ and the table are kept as the historical motivation.

## Step 3 — the approximation (saddle point of F_k at m)
θ > 0 solves m = Σ_{j≤k} j/(e^{jθ}−1);  B = Σ j² e^{jθ}/(e^{jθ}−1)²;  κ₃ = Σ j³ e^{jθ}(e^{jθ}+1)/(e^{jθ}−1)³.
K = k+1, q = e^{−Kθ}, h(u) = −u − log(1−e^{−Ku}):
  h′ = −1 − Kq/(1−q),   h″ = K²q/(1−q)².
  **L̃_k = −θ − log(1−q) + log(1 − (h′² + h″)/(2B) − h′κ₃/(2B²)).**
(Derivation: Gaussian expansion of the ratio of Cauchy integrals around u = θ, with the
weight exp(−Bφ²/2 + iκ₃φ³/6 + …); E_w[φ] ≈ iκ₃/(2B²).)

## Numerical evidence (floating point; not rigorous)
- Unimodality: no violation for n ≤ 100 000 (exact integers for n ≤ 1500).
- Required margins scale as Gumbel predicts: crossing step η ≈ θ² ≈ 1/β².
- Error of L̃ relative to η (max over the window |x| ≤ 3–4):

| n | crossing k | η | max|L−L₀|/η (no corr.) | max|L−L̃|/η |
|---|---|---|---|---|
| 1 000 | 80 | 2.13e-3 | 0.39 | 0.011 |
| 4 000 | 193 | 4.88e-4 | 1.09 | 0.011 |
| 16 000 | 454 | 1.14e-4 | 2.60 | 0.0105 |
| 60 000 | 1004 | 2.93e-5 | 5.48 | 0.0105 |

Requirement is < 0.5η ⇒ ~50× slack for the rigorous remainder bound. The residual looks like
O(θ³ polylog β) (≈ constant fraction of η over this range, decaying slowly).

## Remaining proof obligations
- **Lemma A (main).** Explicit bound |L_k − L̃_k| ≤ ε for n ≥ N₀, k in the window |x| ≤ x₀:
  Taylor/Laplace remainder for the ratio of two Cauchy integrals sharing the weight F_k,
  explicit cumulant bounds (B, κ₃, κ₄ via sums ↔ integrals), central arc |φ| ≤ φ₀,
  minor arcs |φ| > φ₀ (relative contribution e^{−c√n}, as for p(n)).
- **Lemma B.** L̃_k strictly decreasing on the window with step ≥ η_min; locate crossing.
  (Explicit asymptotics of θ(m,k), q, B via Euler–Maclaurin.)
- **Lemma C.** Left (x ≤ −x₀) and right (x ≥ x₀, k ≤ n/2) regions: crude bounds, margins ~θ ≫ η.
- **Lemma D.** Small k (k ≲ n^{1/3}): elementary p(n,k) bounds; k ≥ n/2: exact.
- **Finite check** n < N₀ (done: fc/FiniteCheck.md and the exact-integer fc/FiniteCheckExact.md, n ≤ 2·10⁵): column recurrence, O(N₀²) time, O(N₀) memory; all-positive sums
  ⇒ rigorous floating error bounds (exact integers only at ties).

Scripts: scan.py (unimodality sweep), margins.py, saddle.py, ratio.py.
