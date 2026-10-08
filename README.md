# Unimodality of p(n,k): a computer-assisted proof of the Auluck–Chowla–Gupta conjecture

**Theorem.** For every n ≥ 1 the sequence p(n,1), p(n,2), …, p(n,n) is weakly unimodal,
where p(n,k) is the number of partitions of n into exactly k parts.

This repository contains everything needed to re-check the proof: an exact finite verification for
n ≤ 2·10⁵ and an analytic argument for n ≥ 10⁵ in which every numerical constant is certified with
ball arithmetic (Arb, via python-flint). The two ranges overlap.

## Proof architecture

Notation: β = √(6n)/π, ℓ = log β, μ = βℓ, R_k = p(n,k+1)/p(n,k), L_k = log R_k, θ the saddle point,
q = e^{−(k+1)θ}, and L₂,k the explicit main term of `notes/LemmaA.md` §3.

For n ≥ 10⁵ the index range is covered by five regions, and each R_k gets its sign independently:

| Region | k range | Claim | Lemma (notes/) | Scripts (analytic/) |
|---|---|---|---|---|
| D | 1 ≤ k ≤ ⌊1.7 n^{1/3}⌋ | R_k > 1 | LemmaCD §4 (D1–D3) | ws_cd_D |
| CL | k_D < k ≤ μ − 2β | R_k > 1 | LemmaCD CL-a, CL-b | ws_cd_CLa, ws_mech_cla, ws_cd_CLb, ws_mech_cd |
| W | μ − 2β < k < μ + 2β | sign L_k = sign L₂,k except for ≤ 1 k; L₂ has one sign change | LemmaA_uniform, LemmaB | ws_p12_*, ws_b_*, ws_mech_lemA, ws_mech_window, ws_mech_edges |
| CR | μ + 2β ≤ k ≤ ⌈n/2⌉ | log R_k ≤ −0.614θ | LemmaCD C-right | ws_cd_CR, ws_mech_cd |
| E | ⌈n/2⌉ ≤ k ≤ n−1 | p(n,k+1) ≤ p(n,k) | LemmaCD §5 (p(n,k) = p(n−k)) | — (exact bijection) |

Window logic (Region W):
1. Lemma A: |L_k − L₂,k| ≤ Γ_max θq with Γ_max = 0.4282 (certified 0.42024 in Region I, 0.17605 in Region II).
2. Lemma B (B1): L₂,k − L₂,k+1 ≥ 0.9100·max(θq) > 2Γ_max·max(θq).
3. Lemma B (B2): |L₂,k| > Γ_max θq off the crossing pair.
4. Lemma B (B3): L₂/θ ≥ 5.588 at the left edge of W and ≤ −0.505 at the right edge.

So the signs on W are + … + ? − … −, and the whole row is + … + (?) − … −, which is weakly unimodal.

Shared ingredients: parameter ranges (Lemma W/CW, `ws_mech_window`), total-variation constants of f_r
(P1.2, `ws_mech_tv`), the minor-arc constant c_* ≥ 0.15973 (P2, `ws_mech_minor`; downstream uses 0.15926),
the Lemma A error algebra (`ws_mech_lemA`), θ-monotonicity claims (`ws_mech_mono`) and elementary
calculus facts (`ws_mech_elem`). The `ws_close_*` scripts mechanize or test the items that were
previously hand-proved only (`ws_close_rep`, `ws_close_D`, `ws_close_elem`; written proofs in `notes/PROOFS_FULL.md`) and audit
all thresholds and region coverage for every n ≥ 10⁵ (`ws_close_thresh`). Details: `notes/STATUS.md`, `notes/MECH.md`, `notes/INTERFACE.md`.

## Layout

```
README.md  LICENSE  requirements.txt  run_all.sh
finite/    ws_fc_exact.c (exact integers), ws_fc_check.c (x87 interval arithmetic, independent),
           ws_fc_ref.py (Python big-int reference), Makefile, logs/ (reference logs + SHA256SUMS)
analytic/  proof scripts and their library modules; expected/ holds the reference outputs
notes/     lemma documents, full written proofs (PROOFS_FULL.md) and the internal referee report
paper/     the paper (acg.tex, acg.pdf)
```

Library modules in `analytic/` (imported, not run directly): `ws_p12_common.py`, `ws_p12_eps.py`
(Lemma A ε core; its `TV` table is the ws_mech_tv output), `ws_b_core.py`, `ws_cd_common.py`, and
`lemA_bound.py` (only its saddle-point/cumulant helpers are used, via `ws_p12_eps`, `ws_p12_asym`, `ws_b_check`;
its own pointwise bound is not part of the proof). All scripts are run with `analytic/` as working directory.

## Requirements

* Python 3 with the pinned packages in `requirements.txt` (`pip install -r requirements.txt`):
  python-flint 0.9.0, sympy 1.14.0, numpy 2.5.3, mpmath 1.3.0.
* gcc (tested with 13.3) with OpenMP support (for `ws_fc_check`), GNU make, bash, coreutils.

## How to reproduce

```
./run_all.sh --quick   # analytic part only, ~4.5 min
./run_all.sh           # analytic part + exact finite check to N = 100 000, ~5.5 min
./run_all.sh --full    # analytic part + exact finite check to N = 200 000 (the full proof range), ~13 min
```
`run_all.sh` builds the C programs, runs every analytic script, and diffs each output against
`analytic/expected/` after stripping timing tokens (`time X.Xs` lines and `(Ns)` tags). It then runs the
exact finite check, compares it byte-for-byte with the reference log in `finite/logs/`, and exits
nonzero on any failure or mismatch. Actual outputs go to `analytic/actual/` and `finite/out/`.

The independent interval verifier is not run by default (it is slower):
```
cd finite && make && ./ws_fc_check 200000 und200k.txt log200k.txt   # und200k.txt must be empty
```

## Expected outputs

Key analytic values (full text in `analytic/expected/`):

| Script | Certified value | Needed |
|---|---|---|
| ws_mech_minor | c_* ≥ 0.15973 | ≥ 0.15926 |
| ws_p12_boxes / ws_p12_asym | Γ ≤ 0.42024 / 0.17605 | ≤ 0.4282 |
| ws_b_step / ws_b_asym | Σ ≥ 0.91918 / 0.96039 | > 0.8564 |
| ws_b_edges / ws_mech_edges | left L₂/θ ≥ 5.588, right ≤ −0.505 | > ±Γ_max·q |
| ws_cd_CLb / ws_mech_cd | CL-b margin 5.8206 | > 1 |
| ws_mech_cla / ws_cd_CLa | ρ ≤ 0.0673 / 0.161 | < 0.8977 |
| ws_cd_CR | 0.3851 | < 1 |

Finite check (exact integers; sign counts over pairs (k,k+1), k ≤ ⌈n/2⌉):

| N | '+' | '−' | ties | violations |
|---|---|---|---|---|
| 2 000 | 151 295 | 849 696 | 8 | 0 |
| 100 000 | 85 086 401 | 2 414 963 590 | 8 | 0 |
| 200 000 | 256 705 214 | 9 743 394 777 | 8 | 0 |

SHA-256 (also in `finite/logs/SHA256SUMS`); the exact and interval verifiers produce byte-identical logs:
```
675740a5c82f111a986517bef17ce56ed9b9129665fbc98d6fdff696c4d89e2d  log100k.txt
45ebb662ffbcd5a46af32fe1b1133c78a30e9e5167080a9a148d67c3685e3132  log200k.txt
6f0a4aad6d1c7b01f3a3a80fa174afd292e4eebd8012517e3e32d953a0815e68  ws_fc_exact.c
67e92fa044d5bcc472b87d94cfd213d040ba9a7641cb77b4c4d9986060b4936d  ws_fc_check.c
```

## Runtimes (2-core x86_64, single-threaded, measured 2026-10-08)

| Step | Time |
|---|---|
| ws_mech_cd | 129 s |
| ws_cd_CLb | 48 s |
| ws_mech_minor | 28 s |
| ws_p12_boxes | 16 s |
| ws_p12_minor | 12 s |
| all other analytic scripts | ≤ 5 s each |
| `run_all.sh --quick` total | 4 m 20 s |
| finite exact, N = 100 000 | 63 s (`run_all.sh` total 5 m 24 s) |
| finite exact, N = 200 000 | ≈ 8 m 18 s (from notes/FiniteCheckExact.md; not re-run here) |
| finite interval (ws_fc_check), N = 200 000 | ≈ 10 min |

## Trust base

* Arb ball arithmetic via python-flint 0.9.0 (including polylog, zeta, gamma), sympy 1.14 for exact
  symbolic identities, gcc 13.3 integer code generation, and the two independent C finite-check programs
  (one exact-integer with no floating point, one x87 directed-rounding intervals).
* Every step that is not itself a certified computation has a complete written proof in `notes/PROOFS_FULL.md`
  (representation of R_k, derivative bounds, the [ELEM] one-liners, P1.1/F_r, Lemma B.1–B.4, Lemma C0/C1, Lemma D1–D3,
  Region E), and each is backed by the `ws_close_*` checks (exact sympy identities, Arb enclosures, or exact-integer
  instance tests). Beyond the software above, the proof uses only classical theorems: Taylor's theorem with integral
  remainder, the mean value theorem, the Weierstrass M-test, additivity of total variation, and the identity theorem.
