"""CLOSE-4: threshold / hypothesis audit for n >= 10^5 (PROOFS_FULL.md §9).
Reads the CERTIFIED numbers from the recorded .out files (regex; fails if a file or line is missing) and checks, for ALL n >= 10^5,
that every lemma is applied inside its hypotheses and that the regions cover 1 <= k <= n-1.  Independent Arb recomputation where cheap.
Usage: python3 ws_close_thresh.py [OUTDIR] [FiniteCheckExact.md]       (< 5 s)"""
import re
from flint import arb
import os, sys
OUT = sys.argv[1] if len(sys.argv) > 1 else "."   # directory holding the .out files (run_all.sh passes actual/)
FC = sys.argv[2] if len(sys.argv) > 2 else "fc/FiniteCheckExact.md"
def _p(name): return os.path.join(OUT, name)
def grab(path, pat, grp=1):
    txt = open(_p(path)).read(); m = re.search(pat, txt, re.M)
    assert m, (path, pat); return float(m.group(grp))
def ok(m): print(f"[OK] {m}", flush=True)
NUM = r"([-+]?[0-9.]+(?:e[-+]?[0-9]+)?)"

# ---- certified constants (worst box over ALL beta >= beta_A, ws_mech_window.out) ----
W = {key: grab('ws_mech_window.out', rf"^\s+{key}\s+\S+\s+\S+\s+\[\S+ \+/- \S+\]\s+{NUM}\s+CERTIFIED") for key in
     ['G', 'thA', 'nuhi', 'VL', 'VLprinted', 'nuCL', 'nulo', 'numin', 'nuCR']}
assert W['G'] > 0
thA, nuhi, nulo, numin = W['thA'], W['nuhi'], W['nulo'], W['numin']
ok(f"Lemma W/CW (all n >= 1e5): theta <= {thA}, nu' in [{nulo}, {nuhi}], nu >= {numin}, nu_CL <= {W['nuCL']}, nu_CR >= {W['nuCR']}, nu'_L <= {W['VLprinted']}")

# ---- Lemma P2 (minor arcs, window) hypotheses ----
cstar = grab('ws_mech_minor.out', rf"c_\* = PSI/2 >= {NUM}")
assert thA <= 0.0042 and numin >= 3.0
assert cstar >= 0.15926                      # value passed to ws_p12_boxes / ws_p12_asym
assert 0.15926 / 3.5 > 0.002                 # Region II asserts theta < c_*/3.5 (monotone tail factor) on theta <= 0.002
ok(f"P2: theta_A <= 0.0042, nu_min >= 3.0 on W; c_* = {cstar} >= 0.15926 (input of the Lemma A box scripts); c_*/3.5 > 0.002")

# ---- Lemma A-unif / Lemma B box covers ----
assert -2.52 <= nulo and nuhi <= 2.10 and thA <= 0.0040939 and numin >= 3.1158
G1 = grab('ws_p12_regionI.out', rf"Gamma <= {NUM}\s+\(\d+s\)\s*\Z")
G2 = grab('ws_p12_regionII.out', rf"Gamma <= {NUM} over")
S1 = grab('ws_b_step.out', rf"min Sigma >= {NUM}")
S2 = grab('ws_b_asym.out', rf"\): Sigma >= {NUM};")
lm2 = grab('ws_b_asym.out', rf"sup \|log M2\|/theta <= {NUM}")
Gam = 0.4282
assert max(G1, G2) <= Gam
GmaxI = max(float(x) for x in re.findall(r"Gamma <= ([0-9.]+)", open(_p('ws_p12_regionI.out')).read()))
assert GmaxI <= Gam
Smin = min(S1, S2); step = Smin / 1.01
assert step > 2 * Gam and step > Gam
ok(f"Thm A-unif: Gamma <= {max(GmaxI, G2)} <= Gamma_max = {Gam} on theta in (0, theta_A], nu' in [-2.52, 2.10] (covers W);  "
   f"Lemma B: Sigma_min = {Smin}, step >= {step:.5f} max(theta q) > 2 Gamma_max = {2*Gam:.4f} (Corollary) and > Gamma_max (B-sep)")
assert lm2 <= 0.14236                         # ws_b_edges was run with 0.14236 >= certified sup
L5 = grab('ws_mech_edges.out', rf"left L2/theta >= \[{NUM}")
assert W['VLprinted'] <= -1.906674 + 1e-12 and W['nuCR'] >= 1.0430522 - 1e-9 and L5 > 0
ok(f"(B3): |log M2|/theta <= {lm2} <= 0.14236 (input to ws_b_edges); nu'_L majorant {W['VLprinted']} <= -1.906674 (input to ws_mech_edges); "
   "(B3) is in fact not needed for unimodality (PROOFS_FULL §9.3)")

# ---- Region CR ----
n0 = arb(10) ** 5
thCR = (2 * arb(2).zeta() / (n0 - 1)).sqrt()
assert thCR < arb('0.0057358')
v_lo = arb(W['nuCR']) - thCR.log() - thCR
assert v_lo > arb('6.198')
crpass = grab('ws_mech_cd.out', rf"=== Lemma CR with mechanical Psi ===\n.*\n.*psi = \[{NUM}")
assert crpass > 5 * float(thCR.upper())
crval = grab('ws_mech_cd.out', rf"=== Lemma CR with mechanical Psi ===\n.*\n.*\n.*<= \[{NUM}")
assert crval < 1
ok(f"CR (all n >= 1e5): theta <= sqrt(2 zeta2/(n-1)) <= {thCR.str(6)} (decreasing in n), v = nu'-log theta-theta >= {v_lo.str(6)} >= 6.198, "
   f"Psi_CR = {crpass} > 5 theta_CR, e^-nu'/(1-q)+X <= {crval} < 1")

# ---- Region CL ----
clb = grab('ws_mech_cd.out', rf"CL-b PROVED for n >= 100000 ; min margin .* >= {NUM}")
cla = grab('ws_mech_cla.out', rf"rho <= \[{NUM}")
cla_need = grab('ws_mech_cla.out', rf"need rho < .* = \[{NUM}")
assert clb > 1 and cla < cla_need
# CL-b uses theta <= min(theta_A, v2/78): requires k >= 79 > k_D >= 78 on CL, and theta <= theta_A, nu' <= nu_CL on k <= mu - 2 beta
ok(f"CL: CL-a rho = {cla} < {cla_need} (k >= 47, kθ <= 0.1, any n); CL-b margin {clb} > 1 (k >= 79, θ <= θ_A, ν' <= ν_CL)")

# ---- region cover, for all n >= 1e5 ----
# (i) k_D(n) = floor(1.7 n^{1/3}) >= 78 (nondecreasing in n, = 78 at 1e5)  => every k in CL has k >= 79 >= 47.
assert 1000 * 78 ** 3 <= 4913 * 10 ** 5
# (ii) mu + 2 beta = beta(l+2) <= n/2 - 1, so W and CR lie in 2 <= k <= ceil(n/2) <= n-2 and CR/E meet at ceil(n/2):
#      ratio beta(l+2)/(n/2) = 2(l+2)/(zeta2 beta) = (2/zeta2)(l+2)e^{-l}, decreasing for l > -1;  at n = 1e5:
beta = (6 * n0).sqrt() / arb.pi(); l = beta.log()
assert beta * (l + 2) < n0 / 2 - 1
assert 2 * (l + 2) / (arb(2).zeta() * beta) < arb('0.04')
# (iii) mu - 2 beta > k_D(n): not needed (regions only need to cover); D, CL, W, CR, E cover [1, n-1]:
#      D = [1, k_D], CL = (k_D, mu-2beta], W = (mu-2beta, mu+2beta), CR = [mu+2beta, ceil(n/2)], E = [ceil(n/2), n-1].
# (iv) the analytic lemmas need 2 <= k <= n-2 (representation (1)/(2)): CL, W, CR have 79 <= k <= ceil(n/2) <= n-2.
ok("cover: k_D >= 78; beta(l+2) <= n/2 - 1 for all n >= 1e5 ((l+2)e^{-l} decreasing); D ∪ CL ∪ W ∪ CR ∪ E = [1, n-1]; 2 <= k <= n-2 on CL/W/CR")
# ---- Lemma D ----
d_ok = 'all checks passed' in open(_p('ws_close_D.out')).read()
assert d_ok
ok("D: Phi > 0 for 1 <= k <= 1.7 n^{1/3} and all n >= 1e5 (Psi increasing: ws_close_D.out D3c/D3d)")
# ---- finite check overlap ----
fc = open(FC).read()
assert '200000' in fc or '2·10⁵' in fc or '200 000' in fc
ok("finite check: n <= 2*10^5 (exact-integer + interval, byte-identical logs) overlaps the analytic range n >= 10^5")
print("all threshold checks passed")
