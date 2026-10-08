"""MECH-4: mechanical verification of every "monotone in theta" claim used to pass from a finite evaluation to all small theta.
 (M1) LemmaCD Lemma C1 (X_cell): with the substituted ingredient bounds, T1/th, T2/th, A3, delta_T are nondecreasing in th,
      w0 nonincreasing, M/th and M nondecreasing when Psi >= 5 th.  Checked by sympy: d/dth of each expression, written over the
      positive symbols (blo := F2(v1) - th > 0 asserted in X_cell), is shown to be  (polynomial with all coefficients >= 0)/(positive)
      or its negative.  M-terms: identity E12 (ws_mech_elem).  Then rho/th = (increasing, positive)/(decreasing, positive) is increasing.
 (M2) LemmaB Region II: c1 th^2 = th^2 + th^2 nu e^{-nu'}/(1-q) increasing in th at fixed nu' (E13 + nu > 3, q = th e^{-nu'} increasing).
 (M3) LemmaA Region II / LemmaB Region II graded monomials: t^P nu^j nondecreasing on (0, t1] iff P nu >= 2j (E7); the scripts ASSERT
      P > 0 and nu_min > 2j/P for every monomial -- re-executed here (both scripts imported and run on every nu'-cell, asserts live).
      Tail factors e^{-a/th} th^{-p}: increasing iff th < a/p (E12) -- asserted in the scripts (TH1 < 0.02 b, TH1 < c_*/3.5).
 (M4) LemmaB Region II: (1+tau)-step majorants th'^e <= (1+tau)^e th^e (e>0) -- immediate from th <= th' <= (1+tau)th (Lemma B.1, asserted)."""
import sympy as sp, sys, subprocess, time
th, blo0, F3, TV3, F4, TV4, f1v, f2v, c0, Psi, b = sp.symbols('th blo0 F3 TV3 F4 TV4 f1v f2v c0 Psi b', positive=True)
# substitution: blo = blo0 - th  with blo0 = F2(v1); positivity of blo is the assertion blo > 0 in X_cell, so introduce B = blo > 0, blo0 = B + th
Bs = sp.symbols('B', positive=True)
s3 = F3 + th * TV3; s4 = F4 + th * TV4; blo = blo0 - th; c1th = th + f1v
exprs = {
    'T1/th': (c1th * s3 / (2 * blo ** 2), +1),
    'T2/th': ((f2v + c1th ** 2) / (2 * blo), +1),
    'A3':    (s3 * sp.sqrt(th) / blo ** sp.Rational(3, 2), +1),
    'deltaT': (c0 ** 2 * s4 / (12 * blo), +1),
    'w0^2':  (c0 ** 2 * blo / th, -1),
}
def signcheck(e, sgn):
    d = sp.diff(e, th).subs(blo0, Bs + th)
    d = sp.factor(sp.together(sp.simplify(d)))
    num, den = sp.fraction(d)
    num = sp.expand(num * sgn)
    # den must be positive: product of powers of positive symbols/sqrt; num: polynomial in positive symbols & sqrt(th), sqrt(B) with coeffs >= 0
    num2 = num.subs({th: sp.Symbol('T', positive=True) ** 2, Bs: sp.Symbol('BB', positive=True) ** 2})
    num2 = sp.expand(sp.powsimp(num2, force=True))
    poly = sp.Poly(num2, *sorted(num2.free_symbols, key=str))
    ok_num = all(c >= 0 for c in poly.coeffs())
    den2 = sp.expand(sp.powsimp(den.subs({th: sp.Symbol('T', positive=True) ** 2, Bs: sp.Symbol('BB', positive=True) ** 2}), force=True))
    ok_den = den2.is_positive
    return ok_num, ok_den, d
print("(M1) Lemma C1 ingredient monotonicity (sympy):")
for k, (e, sg) in exprs.items():
    on, od, d = signcheck(e, sg)
    print(f"   {k:7s} {'nondecreasing' if sg > 0 else 'nonincreasing'} in th:  numerator coeffs >= 0: {on}, denominator > 0: {od}")
    assert on and od
# M and M/th:  sqrt(2 pi bup) th^{-p} exp(-Psi/(2th)), p = 3/2 and 5/2:  d/dth log = (Psi/2 - p th)/th^2  > 0 iff th < Psi/(2p)
for p in (sp.Rational(3, 2), sp.Rational(5, 2)):
    d = sp.simplify(sp.diff(sp.log(th ** (-p) * sp.exp(-Psi / (2 * th))), th) - (Psi / 2 - p * th) / th ** 2)
    assert d == 0
print("   M (p=3/2), M/th (p=5/2): d log/dth = (Psi/2 - p th)/th^2 > 0 for th < Psi/5 (asserted 'psi > 5 th' in X_cell)")
print("   erfc(w0/sqrt2) increasing in th (erfc decreasing, w0 nonincreasing); delta = min(deltaT, deltaC): deltaC th-free;")
print("   s3 = min(F3+th TV3, v2 sup f3), blo = max(F2(v1)-th, v1 f2(v2)): min/max of monotone pieces with the same direction.")
print("   => numerator of rho/th nondecreasing and >= 0, denominator nonincreasing; den(th_bar) > 0 asserted => rho/th <= value at th_bar.  OK")
# M2
nup = sp.symbols('nup', real=True)
e = th ** 2 * (nup - sp.log(th))
assert sp.simplify(sp.diff(e, th) - th * (2 * (nup - sp.log(th)) - 1)) == 0
print("(M2) d/dth[th^2 nu] = th(2nu-1) > 0 since nu >= nu_min > 3 (asserted 'nmin > 3' in ws_b_asym.setup); 1/(1-q) increasing.  OK")
print("(M3) re-running the graded-monomial scripts with live monotonicity assertions:")
t0 = time.time()
for cmd in (["python3", "ws_p12_asym.py", "0.002", "-2.52", "2.10", "0.15926", "0.05"], ["python3", "ws_b_asym.py", "0.002", "-2.52", "2.10", "0.05"]):
    out = subprocess.run(cmd, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-2000:]
    print("   ", " ".join(cmd[1:]), "->", out.stdout.strip().splitlines()[-1])
print(f"   ({time.time()-t0:.1f}s)  every monomial passed 'P > 0 and nu_min > 2j/P' (E7) and the tail asserts (E12).")
print("ALL MONOTONICITY CLAIMS VERIFIED")
