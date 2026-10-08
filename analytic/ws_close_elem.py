"""CLOSE-3: mechanization of the remaining derivations (PROOFS_FULL.md §§3-6): the [ELEM] one-liners, P1.1/F_r, Lemma B.1-B.4 +
T2 + Region II ingredient majorants, Lemma C0/C1 error terms, CL-a minor-arc normalisations.
Usage: python3 ws_close_elem.py      (~2-4 min)
Every [OK] line is (S) an exact sympy identity / sign-of-manifestly-signed-factors check, (A) a rigorous Arb enclosure, or
(N) a numerical sanity test of an inequality whose proof is written in PROOFS_FULL.md (labelled as such)."""
import time, random, sympy as sp
from sympy import I, Rational as R
import mpmath as mp
from flint import arb, acb
t0 = time.time()
def ok(m): print(f"[OK] {m}  ({time.time()-t0:.0f}s)", flush=True)

# ============ item 3: [ELEM] one-liners ============
z = sp.symbols('z', positive=True)   # both sides entire in z: identity on (0,inf) => on C (identity theorem)
s_, u, phi = sp.symbols('s u phi')
rem = z ** 3 * sp.integrate((1 - s_) ** 2 / 2 * sp.exp(s_ * z), (s_, 0, 1))
assert sp.simplify(rem - (sp.exp(z) - 1 - z - z ** 2 / 2)) == 0
rem1 = z * sp.integrate(sp.exp(s_ * z), (s_, 0, 1))
assert sp.simplify(rem1 - (sp.exp(z) - 1)) == 0
assert sp.integrate((1 - s_) ** 2 / 2, (s_, 0, 1)) == R(1, 6)
ok("(S) X1  e^z-1-z-z^2/2 = z^3 int_0^1 (1-s)^2/2 e^{sz} ds and e^z-1 = z int_0^1 e^{sz} ds; int (1-s)^2/2 = 1/6  => |.| <= |z|^3/6 if Re z <= 0")
pp = sp.symbols('p', positive=True)
for nn in range(0, 7):
    assert sp.simplify(sp.integrate((pp - u) ** nn / sp.factorial(nn), (u, 0, pp)) - pp ** (nn + 1) / sp.factorial(nn + 1)) == 0
# Taylor with integral remainder, checked symbolically on a generic exponential-polynomial test function
a_ = sp.symbols('a', positive=True)
ftest = sp.exp(I * a_ * u) * u ** 2
for nn in range(0, 5):
    T = sum(sp.diff(ftest, u, r).subs(u, 0) * pp ** r / sp.factorial(r) for r in range(nn + 1))
    Rm = sp.integrate(sp.diff(ftest, u, nn + 1) * (pp - u) ** nn / sp.factorial(nn), (u, 0, pp))
    assert sp.simplify(ftest.subs(u, pp) - T - Rm) == 0
ok("(S) X2  int_0^p (p-u)^n/n! du = p^{n+1}/(n+1)!  (n<=6); Taylor integral-remainder formula re-verified on a test function (n<=4)")

# numerical sanity of the [ELEM] inequalities at small (n, k) with the real saddle
mp.mp.dps = 30
def saddle(m, k):
    S = lambda t: mp.fsum(i / mp.expm1(i * t) for i in range(1, k + 1)) - m
    lo, hi = mp.mpf(10) ** -6, mp.mpf(20)
    for _ in range(110):
        mid = (lo + hi) / 2
        if S(mid) > 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2
random.seed(7)
for (m, k) in [(60, 8), (200, 15), (500, 25)]:
    th = saddle(m, k); x = mp.e ** (-th); K = k + 1; q = x ** K
    kap = {r: mp.fsum(mp.mpf(i) ** r * mp.polylog(1 - r, x ** i) for i in range(1, k + 1)) for r in range(2, 7)}
    cc = {1: 1 + K * q / (1 - q)}; cc.update({r: mp.mpf(K) ** r * mp.polylog(1 - r, q) for r in range(2, 6)})
    B = kap[2]
    psi = lambda p: -1j * m * p + mp.fsum(mp.log(1 - x ** i) - mp.log(1 - x ** i * mp.e ** (1j * i * p)) for i in range(1, k + 1))
    H = lambda p: 1j * p - mp.log(1 - q * mp.e ** (1j * K * p)) + mp.log(1 - q)
    for _ in range(25):
        tt = mp.mpf(random.uniform(0, 3.0)) * th
        om = psi(tt) + B * tt ** 2 / 2
        assert mp.re(om) <= kap[4] * tt ** 4 / 24 + kap[6] * tt ** 6 / 720 + mp.mpf(10) ** -25
        h = H(tt); h1 = 1j * cc[1] * tt
        assert abs(h) <= cc[1] * tt * (1 + mp.mpf(10) ** -20)
        assert abs(mp.e ** h - 1 - h - h ** 2 / 2) <= abs(h) ** 3 / 6 * (1 + mp.mpf(10) ** -20)
        T4 = h1 - cc[2] * tt ** 2 / 2 - 1j * cc[3] * tt ** 3 / 6 + cc[4] * tt ** 4 / 24
        assert abs(h - T4) <= cc[5] * tt ** 5 / 120 * (1 + mp.mpf(10) ** -15) + mp.mpf(10) ** -25
        assert abs(h - h1) <= cc[2] * tt ** 2 / 2 * (1 + mp.mpf(10) ** -15)
        Sx = mp.fsum(i / mp.expm1(i * tt) for i in range(1, k + 1)) if tt > 0 else None
        F1 = lambda v: mp.quad(lambda w: w / mp.expm1(w), [0, v])
        if tt > 0: assert tt ** 2 * Sx >= F1(K * tt) - tt
ok("(N) X3  at the real saddle, (m,k) in {(60,8),(200,15),(500,25)}, 75 random t: Re w <= k4t^4/24+k6t^6/720, |H|<=c1t, "
   "|e^H-1-H-H^2/2|<=|H|^3/6, |H - T4(H)| <= c5t^5/120, |H-H1| <= c2t^2/2, t^2 S(t) >= F1(Kt) - t")

# ============ item 4: P1.1 and F_r ============
v = sp.symbols('v', positive=True)
L = sp.Function('L')            # L(s) stands for Li_s(e^{-v});  d/dv L(s) = -L(s-1)
for r in range(1, 9):
    G = sp.factorial(r) * (sp.zeta(2) - sum(v ** jj / sp.factorial(jj) * sp.Symbol(f'L{2-jj}') for jj in range(r + 1)))
    # derivative with the chain rule d/dv L_s = -L_{s-1}
    dG = 0
    for jj in range(r + 1):
        Ls = sp.Symbol(f'L{2-jj}'); Lm = sp.Symbol(f'L{1-jj}')
        dG += -sp.factorial(r) * (sp.diff(v ** jj / sp.factorial(jj), v) * Ls + v ** jj / sp.factorial(jj) * (-Lm))
    assert sp.expand(dG - v ** r * sp.Symbol(f'L{1-r}')) == 0
ok("(S) F1  d/dv [ r!(zeta2 - sum_{j<=r} v^j/j! Li_{2-j}(e^{-v})) ] = v^r Li_{1-r}(e^{-v})  (telescoping, r = 1..8, using d/dv Li_s(e^{-v}) = -Li_{s-1}(e^{-v}))")
yv = sp.symbols('y', positive=True)
assert sp.simplify(sp.diff(sp.polylog(3, sp.exp(-v)), v) + sp.polylog(2, sp.exp(-v))) == 0
for jj in range(1, 9):
    if jj == 1: term = v * (-sp.log(1 - sp.exp(-v)))
    else: term = v ** jj * sp.expand_func(sp.polylog(2 - jj, sp.exp(-v)))
    assert sp.limit(term, v, 0, '+') == 0, jj
assert sp.limit(sp.polylog(2, sp.exp(-v)), v, 0, '+') == sp.zeta(2) or sp.N(sp.limit(sp.polylog(2, sp.exp(-v)), v, 0, '+') - sp.zeta(2)) == 0
ok("(S) F2  chain rule d/dv Li_s(e^{-v}) = -Li_{s-1}(e^{-v}); limits v^j Li_{2-j}(e^{-v}) -> 0 (j=1..8) and Li_2(e^{-v}) -> zeta(2) as v->0+  => F_r(0)=0")
import sys; sys.path.insert(0, '.')
from ws_p12_common import F_fun, f_fun, EUL
from ws_p12_eps import TV
def f_acb(r):
    def f(xx, an):
        y = (-xx).exp(); Pv = acb(0)
        for c in reversed(EUL[r - 1]): Pv = Pv * y + c
        return xx ** r * y * Pv / (1 - y) ** r
    return f
for r in range(1, 7):
    for (a, b) in [(0.25, 1), (0.25, 3), (1, 9), (3, 30)]:
        Iab = acb.integral(f_acb(r), a, b).real
        Fd = F_fun(r, arb(b)) - F_fun(r, arb(a))
        assert Iab.overlaps(Fd), (r, a, b, Iab, Fd)
ok("(A) F3  rigorous Arb quadrature int_a^b f_r == F_r(b) - F_r(a) (project F_fun) for r=1..6 on 4 intervals each")
# P1.1 numerical sanity (proof in PROOFS_FULL §4): |theta sum f_r(i theta) - F_r(k theta)| <= theta TV_r
for (th, k) in [(arb('0.004'), 900), (arb('0.0012'), 4000), (arb('0.05'), 60)]:
    for r in range(1, 6):
        ssum = th * sum((f_fun(r, th * i) for i in range(1, k + 1)), arb(0))
        assert abs(ssum - F_fun(r, th * k)) < th * TV[r]
ok("(A) P1  theta sum_{i<=k} f_r(i theta) within theta*TV_r of F_r(k theta) at 3 (theta,k), r=1..5 (instance check; proof is §4)")

# ============ item 5: Lemma B ============
tq, Ks = sp.symbols('t K', positive=True)
for kk in range(1, 7):
    Sk = sum(i / (sp.exp(i * tq) - 1) for i in range(1, kk + 1))
    Sk1 = sum(i / (sp.exp(i * tq) - 1) for i in range(1, kk + 2))
    assert sp.simplify(Sk1 - Sk - (kk + 1) / (sp.exp((kk + 1) * tq) - 1)) == 0
    Bk = -sp.diff(Sk, tq)
    f2 = lambda w: w ** 2 * sp.exp(w) / (sp.exp(w) - 1) ** 2
    assert sp.simplify(Bk - sum(f2(i * tq) for i in range(1, kk + 1)) / tq ** 2) == 0
qq = sp.exp(-Ks * tq)
assert sp.simplify(Ks * qq / (1 - qq) - Ks / (sp.exp(Ks * tq) - 1)) == 0
ok("(S) B1  S_{k+1}(t) - S_k(t) = K/(e^{Kt}-1) = c1 - 1 at t=theta; -S_k'(t) = t^{-2} sum f2(it) = b(k,t)/t^3 (k=1..6); f2(u) = u^2 e^u/(e^u-1)^2")
def fr_sym(r, w):
    y = sp.exp(-w); Pv = sum(c * y ** e for e, c in enumerate(EUL[r - 1])); return w ** r * y * Pv / (1 - y) ** r
w = sp.symbols('w', positive=True)
for r in range(2, 6):
    assert sp.simplify(w * sp.diff(fr_sym(r, w), w) - (r * fr_sym(r, w) - fr_sym(r + 1, w))) == 0
    kk = 3
    sr = tq * sum(fr_sym(r, i * tq) for i in range(1, kk + 1)); sr1 = tq * sum(fr_sym(r + 1, i * tq) for i in range(1, kk + 1))
    assert sp.simplify(sp.diff(sr, tq) - ((1 + r) * sr - sr1) / tq) == 0
    srk1 = tq * sum(fr_sym(r, i * tq) for i in range(1, kk + 2))
    assert sp.simplify(srk1 - sr - tq * fr_sym(r, (kk + 1) * tq)) == 0
ok("(S) B4  u f_r' = r f_r - f_{r+1}; d/dt s_r(k,t) = ((1+r) s_r - s_{r+1})/t; s_r(k+1,t) = s_r(k,t) + t f_r((k+1)t)   (r=2..5, k=3)")
th_, nu_, dth_ = sp.symbols('theta nu dtheta', positive=True)
# B.3: nu_{k+1} - nu_k = (K+1)(theta+dth) - K theta = theta + dth + K dth, K = nu/theta
assert sp.simplify((nu_ / th_ + 1) * (th_ + dth_) - nu_ - (th_ + dth_ + (nu_ / th_) * dth_)) == 0
# T2: (1 - q e^{-theta})/(1-q) = 1 + q(1-e^{-theta})/(1-q)
qs = sp.symbols('q', positive=True)
assert sp.simplify((1 - qs * sp.exp(-th_)) / (1 - qs) - (1 + qs * (1 - sp.exp(-th_)) / (1 - qs))) == 0
# Region II: |h'(u)| <= (u+1)e^{-u}/(1-e^{-u})^2,  h = u e^{-u}/(1-e^{-u})
h = u / (sp.exp(u) - 1); hp = sp.diff(h, u)
assert sp.simplify(hp - (sp.exp(u) - 1 - u * sp.exp(u)) / (sp.exp(u) - 1) ** 2) == 0
Nn = sp.exp(u) - 1 - u * sp.exp(u); assert Nn.subs(u, 0) == 0 and sp.simplify(sp.diff(Nn, u) + u * sp.exp(u)) == 0   # Nn <= 0
assert sp.simplify(((u + 1) * sp.exp(u) - (-Nn)) - (2 * sp.exp(u) - 1)) == 0                # (u+1)e^u - |Nn| = 2e^u - 1 > 0
assert sp.simplify((u + 1) * sp.exp(u) / (sp.exp(u) - 1) ** 2 - (u + 1) * sp.exp(-u) / (1 - sp.exp(-u)) ** 2) == 0
# c1 theta^2 = theta^2 + theta^2 nu e^{-nu'}/(1-q) with q = theta e^{-nu'}, nu = nu' + log(1/theta):  d/dtheta [theta^2 nu] = theta(2 nu - 1)
nup = sp.symbols('nup', real=True)
assert sp.simplify(sp.diff(th_ ** 2 * (nup - sp.log(th_)), th_) - th_ * (2 * (nup - sp.log(th_)) - 1)) == 0
# f_r Eulerian form: f_r(u) = u^r e^{-u} P_{r-1}(e^{-u})/(1-e^{-u})^r = u^r Li_{1-r}(e^{-u})
for r in range(1, 7):
    assert sp.simplify(sp.expand_func(sp.polylog(1 - r, yv)) - yv * sum(c * yv ** e for e, c in enumerate(EUL[r - 1])) / (1 - yv) ** r) == 0
ok("(S) B5  B.3 algebra; T2 identity; Region II: h' formula, h' numerator <= 0, |h'| <= (u+1)e^{-u}/(1-e^{-u})^2; d(theta^2 nu)/dtheta = theta(2nu-1); Eulerian Li_{1-r}")

# Numerical validation of B.1-B.4 and of the hull at REAL saddles (Arb bisection for theta_k, theta_{k+1})
from ws_b_core import deltas_box, hull_point, dlogM2_box, M2_formula, TAU
def saddle_arb(m, k, prec_steps=80):
    lo, hi = arb('1e-6'), arb(1)
    S = lambda t: sum((arb(i) / (arb(i) * t).expm1() for i in range(1, k + 1)), arb(0)) - m
    # invariant: S(lo) > 0 > S(hi) certified, so the root lies in [lo, hi] (S strictly decreasing)
    for _ in range(prec_steps):
        mid = (lo + hi) / 2
        val = S(mid)
        if val > 0: lo = mid
        elif val < 0: hi = mid
        else: break                               # sign undetermined at this precision: stop, enclosure still valid
    return lo.union(hi)
def s_vals(k, th): return {r: th * sum((f_fun(r, th * i) for i in range(1, k + 1)), arb(0)) for r in (2, 3, 4, 5)}
cases = [(10 ** 5, 870), (10 ** 5, 1358), (10 ** 5, 1850), (10 ** 6, 3700), (10 ** 6, 5191), (10 ** 6, 6700)]
for (n, k) in cases:
    m = n - k
    th = saddle_arb(m, k); th1 = saddle_arb(m - 1, k + 1)
    nu = (k + 1) * th; nu1 = (k + 2) * th1; nup = nu + th.log()
    dl = deltas_box(th, nup)
    assert th1 > th and th1 < (1 + TAU) * th                      # B.1
    assert (th1 - th) < dl['dth'] and (nu1 - nu) < dl['dnu'] and nu1 > nu      # B.2, B.3
    s0, s1 = s_vals(k, th), s_vals(k + 1, th1)
    for r in (2, 3, 4, 5): assert abs(s1[r] - s0[r]) < dl['ds'][r]   # B.4
    th_h, nu_h, sh = hull_point(th, nup, dl)
    for (val, box) in [(th, th_h), (th1, th_h), (nu, nu_h), (nu1, nu_h)] + [(s0[r], sh[r]) for r in (2, 3, 4, 5)] + [(s1[r], sh[r]) for r in (2, 3, 4, 5)]:
        assert box.contains(val), (n, k)
    tot, _, _ = dlogM2_box(th, nup)
    M0, _ = M2_formula(th, nu, s0[2], s0[3], s0[4], s0[5]); M1, _ = M2_formula(th1, nu1, s1[2], s1[3], s1[4], s1[5])
    assert abs(M1.log() - M0.log()) < tot
ok(f"(N) B6  at real saddles {cases}: B.1 (theta<theta'<(1+tau)theta), B.2-B.4 increments, both endpoints inside the hull box, |dlog M2| < MVT bound")

# ============ item 6: Lemma C0/C1 and CL-a normalisations ============
ph, aa, Bs, dl_, c1s, k3s, c2s, lam, A0 = sp.symbols('phi a B delta c1 kappa3 c2 lambda A0', positive=True)
omd = sp.symbols('omd', positive=True)        # omd = 1 - delta > 0
nrm = sp.sqrt(2 * sp.pi / Bs)
I4 = sp.integrate(ph ** 4 * sp.exp(-omd * Bs * ph ** 2 / 2), (ph, -sp.oo, sp.oo))
assert sp.simplify(c1s * k3s / 6 * I4 / nrm - (c1s * k3s / (2 * Bs ** 2)) / omd ** R(5, 2)) == 0
I2 = sp.integrate(ph ** 2 * sp.exp(-omd * Bs * ph ** 2 / 2), (ph, -sp.oo, sp.oo))
assert sp.simplify((c2s + c1s ** 2) / 2 * I2 / nrm - ((c2s + c1s ** 2) / (2 * Bs)) / omd ** R(3, 2)) == 0
I3 = 2 * sp.integrate(ph ** 3 * sp.exp(-omd * Bs * ph ** 2 / 2), (ph, 0, sp.oo))
assert sp.simplify(k3s / 6 * I3 / nrm - sp.sqrt(2 / sp.pi) * (k3s * Bs ** R(-3, 2)) / (3 * omd ** 2)) == 0
ok("(S) C0a  Gaussian integrals: c1 k3/6 int phi^4 g = T1/(1-d)^{5/2}, (c2+c1^2)/2 int phi^2 g = T2/(1-d)^{3/2}, k3/6 int|phi|^3 g = sqrt(2/pi) A3/(3(1-d)^2) (units sqrt(2pi/B))")
w0 = sp.symbols('w0', positive=True)
assert sp.simplify(sp.integrate(sp.exp(-Bs * ph ** 2 / 2), (ph, -w0 / sp.sqrt(Bs), w0 / sp.sqrt(Bs))) / nrm - sp.erf(w0 / sp.sqrt(2))) == 0
ok("(S) C0b  int_{|phi|<=phi0} e^{-B phi^2/2} = sqrt(2pi/B)(1 - erfc(w0/sqrt2)), w0 = phi0 sqrt B")
# Im psi parity & second-order vanishing, minor arc 2*2pi*sup = sqrt(2pi/B)*2M
Ms = sp.symbols('Msup', positive=True)
assert sp.simplify(2 * 2 * sp.pi * Ms / nrm - 2 * sp.sqrt(2 * sp.pi * Bs) * Ms) == 0
ok("(S) C0c  minor arcs: 2*2pi*sup|chi| = sqrt(2pi/B) * 2M,  M = sqrt(2 pi B) sup|chi|")
# delta_C: (x - x^3/6)^2 - x^2(1 - x^2/3) = x^6/36 ; sum a_i = B phi^2 ; |1-ye^{ia}|^2/(1-y)^2 = 1 + sin^2(a/2)/sinh^2(u/2), y=e^{-u}
xx = sp.symbols('x', positive=True)
assert sp.expand((xx - xx ** 3 / 6) ** 2 - xx ** 2 * (1 - xx ** 3 / 3 / xx)) == xx ** 6 / 36
al = sp.symbols('alpha', real=True)
yy = sp.exp(-u)
lhs = ((1 - yy * sp.cos(al)) ** 2 + (yy * sp.sin(al)) ** 2) / (1 - yy) ** 2
assert sp.simplify(sp.expand_trig((lhs - 1 - sp.sin(al / 2) ** 2 / sp.sinh(u / 2) ** 2).rewrite(sp.exp))) == 0
ii = sp.symbols('i', positive=True)
assert sp.simplify(((ii * ph / 2) ** 2 / sp.sinh(ii * th_ / 2) ** 2 - ph ** 2 * ii ** 2 * sp.exp(ii * th_) / (sp.exp(ii * th_) - 1) ** 2).rewrite(sp.exp)) == 0
ok("(S) C0d  (x - x^3/6)^2 = x^2(1-x^2/3) + x^6/36; |1-ye^{ia}|^2/(1-y)^2 = 1 + sin^2(a/2)/sinh^2(u/2); (i phi/2)^2/sinh^2(i theta/2) = phi^2 i^2 e^{i theta}/(e^{i theta}-1)^2")
# C1 scaled forms
bS, s3S, s4S, c0S, thS, kap3, kap4 = sp.symbols('b s3 s4 c0 theta kappa3 kappa4', positive=True)
Bv = bS / thS ** 3; k3v = s3S / thS ** 4; k4v = s4S / thS ** 5; c1v = sp.symbols('c1v', positive=True)
assert sp.simplify(c1v * k3v / (2 * Bv ** 2) - thS * (c1v * thS) * s3S / (2 * bS ** 2)) == 0
assert sp.simplify((c2s + c1v ** 2) / (2 * Bv) - thS * (c2s * thS ** 2 + (c1v * thS) ** 2) / (2 * bS)) == 0
assert sp.simplify(k3v * Bv ** R(-3, 2) - s3S * thS ** R(1, 2) * bS ** R(-3, 2)) == 0
assert sp.simplify(c0S * thS * sp.sqrt(Bv) - c0S * sp.sqrt(bS / thS)) == 0
assert sp.simplify(k4v * (c0S * thS) ** 2 / (12 * Bv) - c0S ** 2 * s4S / (12 * bS)) == 0
nuS = sp.symbols('nu', positive=True); Kv = nuS / thS; qv = sp.exp(-nuS)
assert sp.simplify((1 + Kv * qv / (1 - qv)) * thS - (thS + nuS / (sp.exp(nuS) - 1))) == 0
assert sp.simplify(Kv ** 2 * qv / (1 - qv) ** 2 * thS ** 2 - nuS ** 2 * sp.exp(nuS) / (sp.exp(nuS) - 1) ** 2) == 0
Psi_ = sp.symbols('Psi', positive=True)
assert sp.simplify(sp.diff(-R(5, 2) * sp.log(thS) - Psi_ / (2 * thS), thS) - (Psi_ - 5 * thS) / (2 * thS ** 2)) == 0
ok("(S) C1a  scaled T1 = theta(c1 theta)s3/(2b^2), T2 = theta(c2theta^2+(c1theta)^2)/(2b), A3 = s3 theta^{1/2} b^{-3/2}, w0 = c0(b/theta)^{1/2}, "
   "delta_T = c0^2 s4/(12b); c1 theta = theta + f1(nu), c2 theta^2 = f2(nu); d/dtheta log(M/theta) = (Psi-5theta)/(2theta^2)")
# rho monotone in its arguments (partial derivatives signed): rho = (T1 (1-d)^{-5/2} + T2 (1-d)^{-3/2} + 2M)/(1 - E - c A3 (1-d)^{-2} - M)
T1s, T2s, Mm, Es, A3s, cc_ = sp.symbols('T1 T2 M E A3 c', positive=True)
num_ = T1s * (1 - dl_) ** R(-5, 2) + T2s * (1 - dl_) ** R(-3, 2) + 2 * Mm
den_ = 1 - Es - cc_ * A3s * (1 - dl_) ** -2 - Mm
rho_ = num_ / den_
e_ = sp.symbols('e', positive=True)              # delta = e/(1+e) parametrises 0 < delta < 1
for var in [T1s, T2s, Mm, Es, A3s, dl_]:
    dn_ = sp.simplify(sp.diff(num_, var).subs(dl_, e_ / (1 + e_)))
    dd_ = sp.simplify(sp.diff(den_, var).subs(dl_, e_ / (1 + e_)))
    assert (dn_ == 0 or dn_.is_positive) and (dd_ == 0 or dd_.is_negative), (var, dn_, dd_)
    assert sp.simplify(sp.diff(rho_, var) - (sp.diff(num_, var) * den_ - num_ * sp.diff(den_, var)) / den_ ** 2) == 0
ok("(S) C1b  rho is increasing in T1, T2, M, delta, A3 and in erfc(w0/sqrt2) whenever the denominator is positive (direct: numerator "
   "increasing, denominator decreasing in each)")
# CL-a normalisations
lB = sp.symbols('lB', positive=True)
assert sp.simplify((2 * sp.sqrt(Bs / (2 * sp.pi)) * sp.integrate(sp.exp(-lam * Bs * ph ** 2 / 2), (ph, A0, sp.oo))).rewrite(sp.erf) - (sp.erfc(A0 * sp.sqrt(lam * Bs / 2)) / sp.sqrt(lam)).rewrite(sp.erf)) == 0
kq = sp.symbols('k', positive=True)
assert sp.simplify(sp.integrate(ph ** (-kq), (ph, 10, sp.oo), conds='none') - 10 ** (1 - kq) / (kq - 1)) == 0
assert sp.simplify(2 / sp.sqrt(2 * sp.pi) - sp.sqrt(2 / sp.pi)) == 0
assert sp.simplify(sp.sqrt(Bs / (2 * sp.pi)) * 2 * sp.pi - sp.sqrt(2 * sp.pi * Bs)) == 0
ok("(S) CLa  2 sqrt(B/2pi) int_a^inf e^{-lam B phi^2/2} = erfc(a sqrt(lam B/2))/sqrt(lam);  int_10^inf c^{-k} = 10^{1-k}/(k-1) (k>1); N34 prefactor sqrt(2 pi B)")
print(f"all checks passed ({time.time()-t0:.0f}s)")
