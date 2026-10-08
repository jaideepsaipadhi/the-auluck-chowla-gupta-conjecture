# Referee report: "A proof of the Auluck–Chowla–Gupta conjecture on partitions into k parts"

Standard applied: Ramanujan J. / Annals of Combinatorics, hostile. Paper `paper/acg.tex` (commit 1167afc), checked against the
repository code and against a fresh run of `run_all.sh --quick` (ALL CHECKS PASSED) plus independent re-runs of
`finite/ws_fc_exact` at N = 10^5 (log SHA-256 matches `675740a5…`) and N = 2·10^5.

## Verdict

**No FATAL issue found.** I re-derived every step listed below by hand. Two stated numbers in the text were false as written
(they did not affect validity), and one appendix formula did not match the code. All of these, and the minor and presentation
items, are fixed in this revision. No certified constant was changed.

What I re-derived independently:
- Lemma 2.1, Lemma 2.2, Proposition 2.3 (both forms of the representation, including the step from the dominated exchange to
  `e^H = G/g_0`), and (2.6) for |χ|.
- The closed form of F_r (the telescoping and both limits), the scaled identities (2.8), and Lemmas 3.1–3.2 and Corollary 3.3.
- Lemma 4.1 (M_model − M_2 = Num/N_0). I checked this **independently in SymPy**, with D_0 = N_0 and Ñ_0 real, using my own
  script rather than the repository's.
- Lemma 4.2 term by term: the three-way split of e^ω − Q, the domination by T_ω, the decomposition of e^H − P_H, and the integral
  identity.
- The tails (4.4), the minor arcs (4.5), and Proposition 4.3, including the algebra for Ñ/D − Ñ_0/D_0.
- Lemma 4.4 (a) and (b): the Dirichlet identity, the Abel summation, the bound on (N_d+1)θ, and the sinc step.
- Lemma 5.1: the integral sandwich, the bound on T(v), the reduction to G (including the factor β/ℓ), and all four ν/ν′ bounds.
- Lemmas 5.2–5.4 and Corollary 5.5.
- The decomposition in Theorem B; Lemma 6.2 (i)–(iv), including ∂_t s_r = ((1+r)s_r − s_{r+1})/t; and Corollary 6.3, including
  the case analysis for k ≤ k*−1 and k ≥ k*+2.
- Lemma 7.1 (both pieces of the numerator, and the denominator constant √(2/π)A_3/(3(1−δ)²)), Lemma 7.2, and (7.1).
- Lemma 7.3: monotonicity of ρ, and d log(M/θ)/dθ.
- Proposition 7.4 numerically: θ_CR = 0.0057358, v ≥ 6.198, and 0.35310 + 0.0321 ≤ 0.3852.
- Proposition 7.5: T_1, T_2, A_3 and w_0 bounds, and N_1, N_2 and N_34, including θ√B ≤ √k and B ≤ k³/v².
- Lemma 7.7 and Proposition 7.8 (Steps 1–4, ∂_kΦ, Ψ_D′, and Ψ_D(s_A) ≈ 0.3046), and Lemma 7.9.
- The cover argument in Section 8: β(ℓ+2)/(n/2) = 0.0370 at n = 10^5.
- The finite algorithm in `ws_fc_exact.c`: the pair range 2(k−1) ≤ n+1 ⇔ k−1 ≤ ⌈n/2⌉; the in-place comparison order, including
  the m < k pre-loop starting at k−3; the limb budget (27 limbs at 2·10^5); and the abort on overflow.

Numbers in the paper versus script output: every certified value in Table 2 and in the text matches `analytic/expected/` and
my fresh run. Checked: G ≥ 1.0092 on 5486 boxes; θ_A, ν′ range, ν_min, ν_CL and ν_CR; TV_1…TV_6; c_* = 0.15973; Γ = 0.42024 and
0.17605; Σ = 0.91918 on 8526 boxes and 0.96039; |log M_2|/θ ≤ 0.14215; 5.588 and −0.505; X ≤ 0.0321; Ψ = 0.33808; 0.3852;
N_1, N_2, N_34 and ρ for CL-a; the CL-b margin 5.82; Ψ_D(s_A) = 0.30456; the finite counts; and mode 2042 at n = 2·10^5.
Exceptions are items 1, 2 and 6 below.

What is certified versus what the paper claims: the central inequalities are all Arb, with live `assert`s on arb comparisons
(an arb `<` is true only when it is certain). Float is used only to choose the grid. The hypotheses δ_0 < 1, N_0 − E_D > 0 and
Δ < M_2 are asserted in `ws_p12_eps.core` on every box. The monotonicity of the monomials is asserted in `ws_p12_asym` and
`ws_b_asym`. One weakness, now fixed: the asserts were disabled silently if `PYTHONOPTIMIZE` was set in the environment (item 8).

## Issues

### FATAL
None.

### SERIOUS (statements false or inconsistent as written; validity unaffected; all fixed)

1. **Corollary 4.5 (c_*), wrong c_a and cell count.** The text said c_a = 1.331, with L ≤ 1.5105 and "400 geometric cells of
   [0.4, 1.331]". With c_a = 1.331 one gets L = 1.51103 > 1.5105. The code (`ws_mech_minor.py`) actually uses
   c_a = 0.4·150^{96/400} = 1.33144…, giving L = 1.5105. Its 400 cells cover [0.4, 60] at ratio 150^{1/400}, so only 96 of them
   lie in [0.4, c_a]. The certified Ψ is unaffected (with L = 1.5110, ∫_L^3 log coth(u/2) du = 0.34419, which is still above
   0.31946). *Fixed:* the text now states c_a exactly, says "96 geometric cells of ratio 150^{1/400}", and uses the name L_*.
2. **Appendix A, tail terms, extra factor √2.** The text wrote √2·e^{−w_0²/4}·I_{1/2}[·], but I_A as defined already integrates
   against e^{−Aw²/2}, so I_{1/2} already includes the √2. The code (`ws_p12_eps.core`: `tailfac = e^{-w0^2/4}·√2` times
   `poly_abs_integral(·, 1/2)`, which is an expectation under N(0, 2)) evaluates exactly e^{−w_0²/4}·I_{1/2}[·]. So the appendix,
   which claims to record ε_k "as it is evaluated", did not match the code. The appendix value was larger, so it would still have
   been a valid bound. *Fixed:* the formula is corrected, and the parenthetical now explains how I_A relates to the
   expectation-form moments the code computes.

### MINOR (fixed unless stated otherwise)

3. **Corollary 6.3: k\* undefined when no k ∈ W has L_{2,k} > 0.** *Fixed:* added the convention k* := min(W∩Z) − 1.
4. **Proposition 7.5 (CL-a): missing the hypothesis k ≤ n−2** that Proposition 2.3 needs. *Fixed:* the statement is now
   "47 ≤ k ≤ n−2". The assembly already guarantees this.
5. **Proposition 7.6 (CL-b): wrong cell count.** The text said "33 geometric cells and a final cell [8.76, ∞)". The script
   uses 32 cells starting at 0.1 and the final cell [0.1·1.15^{32}, ∞) = [8.7565…, ∞). Writing 8.76 suggested a gap
   [8.7565, 8.76]. *Fixed.*
6. **Proposition 7.4 (CR): source of Ψ.** Ψ = 0.33808 comes from `ws_mech_cd` (the mechanical minor arcs). `ws_cd_CR` uses an
   older Ψ route and prints 0.32647, with the same conclusion. The text now names `ws_mech_cd` and states θ_max = θ_CR.
7. **Region I of Theorem A: 0.15926 was used without explanation** in place of c_* = 0.15973. *Fixed:* the text now says it is
   a rounded-down value.
8. **The asserts are the certificate.** `python -O` or `PYTHONOPTIMIZE` would remove them silently and still print "PROVED".
   *Fixed:* `run_all.sh` now runs `unset PYTHONOPTIMIZE`. *Remaining recommendation:* replace the bare `assert`s with explicit
   `if not …: sys.exit(1)` checks.
9. **Region II monomial condition.** The text said the program asserts "P > 0 and ν_min > 2j/P" for every monomial. The code
   also allows constant monomials (P = j = 0). *Fixed:* the wording in Section 5.3 and in Appendix B.
10. **Section 9.2 overstated the driver.** `run_all.sh` does not run the interval program F2, and it reaches 2·10^5 only with
    `--full` (the default is 10^5). *Fixed.*
11. **Trust base listed "the identity theorem for entire functions",** which no step of the paper uses. *Removed.*
12. **Not fixed, needs the author: the repository URL** (`\repo`, line 11) is a placeholder marked TODO. It must point to a
    real public archive before submission, ideally a Zenodo DOI pinned to the commit that is refereed.

### PRESENTATION

13. **Symbol clashes in Lemma 4.4.** N_0 meant both the Gaussian normalisation and ⌈1/sin d⌉, and L meant both log R and the
    lower limit. *Fixed:* renamed to N_d and L_*.
14. **Δ was overloaded.** It meant both the error quantity of Proposition 4.3 and the increment from k to k+1 in Section 6.
    *Fixed:* Section 6 now writes the increments explicitly as x′ − x.
15. **Unreduced fractions in (2.9):** 105/144 and 945/1296 are now written as 35/48. Same values; the code is unchanged.
16. **Dead notation:** D_0^* in Proposition 4.3 is removed.
17. **Bibliography.** `% VERIFY` flags remain on ACG42, EL41, Arb (page range), Sz51, Sz53 and Sz90. I could not confirm these
    against the originals online. They agree with how these works are standardly cited, but the author must check them against
    a library copy before submission.
    - Confirmed, and flags removed: Erdős 1946 (Bull. AMS 52, 185–188). The PDF on Rényi's site also confirms the intro claim
      that Erdős was "unable to prove or disprove" the conjecture.
    - Confirmed, and flag removed: Canfield, EJC 4(2), R6. Note that this paper is about P(n, k), partitions into *at most* k
      parts. The intro wording "a different derivation of Szekeres' formula" is accurate.
    - The novelty sentence is suitably hedged ("To our knowledge…"). Its `% VERIFY` flag should stay until the author has
      checked MathSciNet/zbMATH for any explicit-N_0 proof after 1990.
18. **README** used labels that conflict with the paper: "Lemma A/B", and "B2" meaning separation rather than Region II.
    *Fixed:* the window logic now uses the paper's Theorem A, Theorem B, Corollary and Remark, and the README names the paper
    as the primary account.
19. **Readability (no change needed).** The paper is readable as a journal paper: definitions come before use, and every
    computer-assisted step is isolated in Table 2. Script names in the body are appropriate for a computer-assisted proof.
    The remaining lab-notebook residue lives only in `notes/`, not in the paper.

## Reproduction

- `./run_all.sh --quick` (before and after the edits): ALL CHECKS PASSED, about 6 min on this machine.
- `finite/ws_fc_exact 100000`: PASS, with +85 086 401, −2 414 963 590 and 8 ties; the log SHA-256 matches the README.
- `finite/ws_fc_exact 200000`: PASS (27 limbs), +256 705 214, −9 743 394 777, 8 ties; the log SHA-256 `45ebb662…` matches `finite/logs/log200k.txt`.
- The paper compiles cleanly: no undefined references and no overfull boxes.
