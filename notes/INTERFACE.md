# ACG proof — shared interface (all workstreams must use these conventions)

Read ACG_attack.md first (same folder) for background and numerics.

## Objects
n ≥ 1. a_k = p(n,k) = f_k(n−k), f_k(m) = [x^m] F_k(x), F_k(x) = ∏_{j=1}^{k} (1−x^j)^{−1}.
m := n − k.  R_k := a_{k+1}/a_k (defined for 1 ≤ k ≤ n−1).
Exact identity: R_k = [x^m](g_k F_k)/[x^m]F_k with g_k(x) = x/(1−x^{k+1}).

Cauchy form, u = θ + iφ, |φ| ≤ π:
  f_k(m) = (1/2π) ∫_{−π}^{π} e^{m u} F_k(e^{−u}) dφ,
  [x^m](g_k F_k) = (1/2π) ∫ e^{m u} F_k(e^{−u}) g_k(e^{−u}) dφ.

Saddle θ = θ(m,k) > 0:  m = Σ_{j≤k} j/(e^{jθ}−1).
  B  = Σ_{j≤k} j² e^{jθ}/(e^{jθ}−1)²
  κ₃ = Σ_{j≤k} j³ e^{jθ}(e^{jθ}+1)/(e^{jθ}−1)³
  κ₄ = Σ_{j≤k} j⁴ e^{jθ}(e^{2jθ}+4e^{jθ}+1)/(e^{jθ}−1)⁴
Weight: log[e^{mu}F_k(e^{−u})] = const − Bφ²/2 + iκ₃φ³/6 + κ₄φ⁴/24 + …  (u = θ+iφ)
K = k+1, q = e^{−Kθ}, h(u) = −u − log(1 − e^{−Ku}), h′ = −1 − Kq/(1−q), h″ = K²q/(1−q)².
Explicit main term used by the proof (2026-10-08):  L₂,k := −θ − log(1−q) + log M₂, with M₂ the closed form of LemmaA.md §3
(cumulants κ₃..κ₅, c₁..c₄). The older L̃_k := −θ − log(1−q) + log(1 − (h′²+h″)/(2B) − h′κ₃/(2B²)) of ACG_attack.md is historical only.
L_k := log R_k.

Scaling: β = √(6n)/π, μ = β log β, x = (k−μ)/β. θ ≈ 1/β near the crossing.

## Regions (k ranges, for fixed n) — FIXED split, workstreams cover these:
  D  (small):  1 ≤ k ≤ k_D(n) := ⌊n^{1/3}⌋        → show R_k > 1      [Workstream 3]
  CL (left):   k_D < k ≤ μ − 2β                   → show R_k > 1      [Workstream 1]
  W  (window): μ − 2β < k < μ + 2β                → |L_k − L₂,k| ≤ ε_k  [Workstream 1]
                                                    and L₂ crossing structure [Workstream 2]
  CR (right):  μ + 2β ≤ k ≤ ⌈n/2⌉                  → show R_k < 1      [Workstream 1]
  E  (end):    k ≥ n/2: a_k = p(n−k), nonincreasing (exact; trivial)
If a workstream needs a different split boundary, it must say so explicitly and give
the boundary as a formula in n.

## Final logic
For n ≥ N₀: unimodality follows if (D, CL give R>1), (CR, E give R≤1), and in W:
  |L_k − L₂,k| ≤ ε_k ≤ Γ_max θq  (LemmaA_uniform), L₂ strictly decreasing on W with step ≥ 0.9100·max(θq) > 2Γ_max·max(θq),
  |L₂,k| > Γ_max θq off the crossing pair, L₂ > ε at the left end of W and < −ε at the right end (LemmaB)
  ⇒ one crossing, ≤ 1 undetermined k (proof: LemmaA_uniform §0 Corollary).
As finally proved: N₀ = 10⁵, and D's boundary is k_D = ⌊1.7n^{1/3}⌋ (≥ 78 for n ≥ 10⁵; LemmaCD).
For n ≤ 2·10⁵: rigorous computer check [Workstream 4], fc/FiniteCheck.md (x87 interval) and the independent exact-integer
verifier fc/FiniteCheckExact.md (ws_fc_exact.c, PASS for all n ≤ 200 000). See STATUS.md for the assembled theorem.

## Rigour standard
Every inequality must have explicit constants. Numerical evaluation of bounds must be
rigorous (interval arithmetic via mpmath.iv or python-flint arb, or exact rationals), or
clearly flagged as non-rigorous. Every claimed threshold must be a formula or a verified
number, never "for n sufficiently large".
Deliverable per workstream: a markdown file WS<i>_*.md in this folder with statements,
complete proofs, constants, the resulting threshold N_i, and scripts used.
