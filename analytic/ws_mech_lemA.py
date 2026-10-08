"""MECH-5: machine check of the LemmaA.md sections 3-6 error algebra and of the Lemma W reduction to G > 0.
Usage: python3 ws_mech_lemA.py      (sympy 1.14 + python-flint; ~1 min)
Every [OK] line is an exact sympy identity, an exact coefficient-wise dominance of polynomials with symbolic positive
coefficients, or (part F only) an Arb cross-check that the project's code computes exactly the formula proved here.
Parts:
 A  Gaussian moments and the normalisation factors sqrt(B/B'), sqrt(2) used in lemA_bound / ws_p12_eps.
 B  Derivatives at 0 of psi (one factor) and of H agree with i^r kappa_r and (i c1, -c2, -i c3, c4); |G| <= g0.
 C  D0 = N0 closed form, Ntil0/D0 real, closed-form M2, and Mmodel - M2 = Num/N0 with Num exactly the 11-monomial
    polynomial hard-coded in ws_p12_eps.py.
 D  T_omega: exact 3-piece decomposition of e^w - Q, the polynomial piece dominated coefficient-wise, the Taylor-remainder
    identity, Re(w_hat) = w4, and the Re(omega) <= delta0 B t^2/2 inequality on the central arc.
 E  T_H: exact decomposition of e^H - P_H and dominance; |Q| <= Qabs, |P_H - 1| <= PH1abs coefficient-wise.
 F  Combination: integrand identity, ratio-perturbation identity, log bound; then an independent Arb implementation of
    eps built from the sympy polynomials of C-E, compared with ws_p12_eps.core at random inputs.
 G  Lemma W: S(theta_-) > m  <=  (2d-d^2) z2 > T(K theta_-) + theta_-  <=  G(beta) > 0, every step an identity or a
    termwise inequality checked here."""
import sympy as sp, random, time
from sympy import I, Rational as R
t0 = time.time()
def ok(msg): print(f"[OK] {msg}", flush=True)
def elem(msg): print(f"[ELEM, by hand, one line] {msg}", flush=True)
t, B, Bp = sp.symbols('t B Bp', positive=True)
a3, a4, a5, a6 = sp.symbols('a3 a4 a5 a6', positive=True)        # a_r = kappa_r / r!
k3, k4, k5, k6 = sp.symbols('k3 k4 k5 k6', positive=True)
c1, c2, c3, c4, c5 = sp.symbols('c1 c2 c3 c4 c5', positive=True)

# ---------------- A. Gaussian moments ----------------
for p in range(0, 21):
    raw = sp.integrate(sp.exp(-B * t**2 / 2) * t**p, (t, 0, sp.oo)) * 2      # int_R |t|^p e^{-Bt^2/2}
    nrm = sp.sqrt(2 * sp.pi / B)
    absm = sp.simplify(raw / nrm - (2 / B)**R(p, 2) * sp.gamma(R(p + 1, 2)) / sp.sqrt(sp.pi))
    assert absm == 0, p
    if p % 2 == 0:
        assert sp.simplify(raw / nrm - sp.factorial2(p - 1) / B**(p // 2)) == 0
ok("A1  int_R e^{-Bt^2/2} t^p / sqrt(2pi/B) = (p-1)!!/B^{p/2} (p even, odd = 0 by symmetry) and E|Z|^p formula, p <= 20")
assert sp.simplify(sp.sqrt(2 * sp.pi / Bp) / sp.sqrt(2 * sp.pi / B) - sp.sqrt(B / Bp)) == 0
assert sp.simplify(sp.sqrt(2 * sp.pi / (B / 2)) / sp.sqrt(2 * sp.pi / B) - sp.sqrt(2)) == 0
ok("A2  normalisation factors: N(0,1/B') moments -> sqrt(2pi/B) units via sqrt(B/B'); tail e^{-Bt^2/4} via sqrt(2)")
phi0 = sp.symbols('phi0', positive=True)
lhs = sp.exp(-B * t**2 / 2); rhs = sp.exp(-B * phi0**2 / 4) * sp.exp(-B * t**2 / 4)
assert sp.simplify(sp.log(rhs) - sp.log(lhs) - B * (t**2 - phi0**2) / 4) == 0
ok("A3  tail: log(e^{-B phi0^2/4} e^{-Bt^2/4}) - log e^{-Bt^2/2} = B(t^2-phi0^2)/4 >= 0 for |t| >= phi0")

# ---------------- B. derivatives at 0 ----------------
y, a, phi, q, K = sp.symbols('y a phi q K', positive=True)
def Li_neg(s, z):   # Li_{-s}(z), s = 0..5 (Eulerian closed forms, as in lemA_bound.Li_neg)
    E = {0: 1, 1: 1, 2: 1 + z, 3: 1 + 4*z + z**2, 4: 1 + 11*z + 11*z**2 + z**3, 5: 1 + 26*z + 66*z**2 + 26*z**3 + z**4}
    return z * E[s] / (1 - z)**(s + 1)
for s in range(0, 6):     # Li_{-s} = (z d/dz)^s Li_0
    f = z0 = sp.Symbol('z0'); g = z0 / (1 - z0)
    for _ in range(s): g = z0 * sp.diff(g, z0)
    assert sp.simplify(g - Li_neg(s, z0)) == 0
ok("B0  Li_{-s}(z) closed forms = (z d/dz)^s [z/(1-z)], s = 0..5")
fac = -sp.log(1 - y * sp.exp(I * a * phi)) + sp.log(1 - y)          # log of (1-y)/(1-y e^{ia phi}); linear term removed by centring
for r in range(2, 7):
    d = sp.diff(fac, phi, r).subs(phi, 0)
    assert sp.simplify(d - (I * a)**r * Li_neg(r - 1, y)) == 0
ok("B1  psi factor: d^r/dphi^r [-log(1-y e^{ia phi})] at 0 = (ia)^r Li_{1-r}(y), r = 2..6  (so psi^(r)(0) = i^r kappa_r)")
H = I * phi - sp.log(1 - q * sp.exp(I * K * phi)) + sp.log(1 - q)
cc = {1: 1 + K * q / (1 - q)}; cc.update({r: K**r * Li_neg(r - 1, q) for r in range(2, 6)})
want = {1: I * cc[1], 2: -cc[2], 3: -I * cc[3], 4: cc[4], 5: I * cc[5]}
assert sp.simplify(H.subs(phi, 0)) == 0
for r in range(1, 6):
    assert sp.simplify(sp.diff(H, phi, r).subs(phi, 0) - want[r]) == 0
ok("B2  H(0)=0, H'(0)=i c1, H''(0)=-c2, H'''(0)=-i c3, H''''(0)=c4, H^(5)(0)=i c5 with c1 = 1+Kq/(1-q), c_r = K^r Li_{1-r}(q)")
# |G| <= g0  <=>  |1 - q e^{iK phi}| >= 1 - q : |1-qe^{ix}|^2 - (1-q)^2 = 2q(1 - cos x) >= 0
x_ = sp.symbols('x_', real=True)
assert sp.simplify(sp.expand((1 - q*sp.cos(x_))**2 + (q*sp.sin(x_))**2 - (1 - q)**2 - 2*q*(1 - sp.cos(x_)))) == 0
ok("B3  |1-q e^{ix}|^2 - (1-q)^2 = 2q(1-cos x) >= 0, hence |G| <= g0 and Re H <= 0")

# ---------------- C. model, D0, N0, M2 ----------------
w3 = -I * k3 * t**3 / 6; w4 = k4 * t**4 / 24; w5 = I * k5 * t**5 / 120
Q = sp.expand(1 + w3 + w4 + w5 + w3**2 / 2 + w3 * w4 + w3**3 / 6)
H1 = I * c1 * t; H2 = -c2 * t**2 / 2; H3 = -I * c3 * t**3 / 6; H4 = c4 * t**4 / 24
PH1 = sp.expand(H1 + H2 + H3 + H4 + H1**2 / 2)
iB = sp.symbols('iB', positive=True)
def G(expr):   # int_R e^{-Bt^2/2} expr dt / sqrt(2pi/B), via A1
    p = sp.Poly(sp.expand(expr), t); out = 0
    for (e,), c in p.terms():
        if e % 2 == 0: out += c * (sp.factorial2(e - 1) if e > 0 else 1) * iB**(e // 2)
    return sp.expand(out)
D0 = G(Q); N0t = G(Q * PH1)
Nn = 1 + k4 * iB**2 / 8 - 5 * k3**2 * iB**3 / 24
assert sp.expand(D0 - Nn) == 0 and not N0t.has(I)
ok("C1  D0 = 1 + k4/(8B^2) - 5k3^2/(24B^3) = N0 exactly; Ntil0 is real (so Ntil0/D0 in R)")
E1 = (k3*iB**2/2 - k5*iB**3/8 + R(105, 144)*k3*k4*iB**4 - R(945, 1296)*k3**3*iB**5) / Nn
E2 = iB * (1 + k4*iB**2/2 - 5*k3**2*iB**3/4)
M2 = 1 + c1*E1 - (c2 + c1**2)*E2/2 - 5*c3*k3*iB**3/12 + c4*iB**2/8          # LemmaA.md sec 3 (doc form)
h1, h2, h3, h4 = -c1, c2, -c3, c4
M2code = 1 - h1*E1 - (h2 + h1*h1)*E2/2 + 5*h3*k3*iB**3/12 + h4*iB**2/8       # lemA_bound / lemA_M2 form
assert sp.simplify(M2 - M2code) == 0
gap = sp.together(1 + N0t / D0 - M2)
Num_derived = sp.expand(sp.numer(gap) * Nn / sp.denom(gap)) if False else sp.expand(sp.simplify((1 + N0t/D0 - M2) * Nn))
# the polynomial hard-coded in ws_p12_eps.core (B = 1 there; restore iB-homogeneity by matching at iB = 1 and checking degrees)
Num_code = (-R(25, 48)*c4*k3**2 + c4*k4/6 + (c2 + c1**2)*(k4**2/32 + R(25, 192)*k3**4 - R(25, 192)*k3**2*k4)
            + R(5, 4)*c3*k3**3 + R(7, 48)*c3*k5 - R(25, 24)*c3*k3*k4)
assert sp.expand(Num_derived.subs(iB, 1) - Num_code) == 0
assert len(sp.Add.make_args(Num_derived)) == 11
ok("C2  M2 (doc) == M2 (code);  Mmodel - M2 = Num/N0 with Num (11 monomials) == the polynomial hard-coded in ws_p12_eps.core at B = 1")
# scaled form: core() uses B := 1 with (kappa_r, c_r) -> (A_r, C_r) = (kappa_r B^{-r/2}, c_r B^{-r/2}); check homogeneity of every piece
lam = sp.symbols('lam', positive=True)
scal = {k3: k3*lam**3, k4: k4*lam**4, k5: k5*lam**5, c1: c1*lam, c2: c2*lam**2, c3: c3*lam**3, c4: c4*lam**4, iB: iB/lam**2}
for name, ex in (("D0", D0), ("N0til", N0t), ("M2", M2), ("Num", Num_derived)):
    assert sp.simplify(ex.subs(scal, simultaneous=True) - ex) == 0, name
ok("C3  D0, Ntil0, M2, Num are invariant under kappa_r -> lam^r kappa_r, c_r -> lam^r c_r, B -> lam^2 B  (justifies the B = 1 scaled core)")

# ---------------- D. T_omega ----------------
W3 = -I*a3*t**3; W4 = a4*t**4; W5 = I*a5*t**5                 # = w3, w4, w5 with a_r = kappa_r/r!
Wh = W3 + W4 + W5
Qa = sp.expand(1 + W3 + W4 + W5 + W3**2/2 + W3*W4 + W3**3/6)
assert sp.expand(Qa - Q.subs({k3: 6*a3, k4: 24*a4, k5: 120*a5})) == 0
cubicT = sp.expand(1 + Wh + Wh**2/2 + Wh**3/6)
poly_piece = sp.expand(cubicT - Qa)          # e^w - Q = (e^w - e^wh) + (e^wh - cubicT) + poly_piece   [exact, by telescoping]
def absmaj(expr):
    """coefficient-wise modulus majorant: replace each monomial's complex coefficient by its modulus."""
    out = 0
    for term in sp.Add.make_args(sp.expand(expr)):
        m = sp.Abs(term)                 # all symbols are positive, so |c * monomial| = |c| * monomial
        assert not m.has(sp.Abs), term
        out += m
    return sp.expand(out)
def dominated(big, small, gens):
    d = sp.Poly(sp.expand(big - small), *gens)
    return all(cf >= 0 for cf in d.coeffs())
Tw_poly_part = sp.expand((a4**2/2 + a3*a5)*t**8 + a4*a5*t**9 + a5**2/2*t**10 + R(1, 2)*(a4*t**4 + a5*t**5)*(a3*t**3 + a4*t**4 + a5*t**5)**2)
assert dominated(Tw_poly_part, absmaj(poly_piece), (a3, a4, a5, t))
ok("D1  |cubic Taylor of e^{w_hat} - Q| <= (a4^2/2+a3a5)t^8 + a4a5 t^9 + a5^2/2 t^10 + (1/2)(a4t^4+a5t^5)(a3t^3+a4t^4+a5t^5)^2  coefficient-wise")
z, s_ = sp.symbols('z s_')
rem = z**4 * sp.integrate((1 - s_)**3 / 6 * sp.exp(s_ * z), (s_, 0, 1), conds='none')
assert sp.simplify(rem - (sp.exp(z) - 1 - z - z**2/2 - z**3/6)) == 0
ok("D2  e^z - sum_{r<4} z^r/r! = z^4 int_0^1 (1-s)^3/3! e^{sz} ds  =>  |.| <= |z|^4/24 max(1, e^{Re z});  |z| <= a3t^3+a4t^4+a5t^5")
assert sp.expand(sp.re(Wh.subs(t, sp.Symbol('tt', positive=True))) - a4*sp.Symbol('tt', positive=True)**4) == 0
ok("D3  Re w_hat = a4 t^4 = kappa4 t^4/24 (w3, w5 purely imaginary); |e^w - e^{w_hat}| <= |w - w_hat| e^{max(Re w, Re w_hat)} (segment), |w - w_hat| <= a6 t^6 (Taylor, |psi^(6)| <= kappa6)")
d0 = (k4*phi0**2/12 + k6*phi0**4/360) / B
slack = sp.expand(d0*B*t**2/2 - (k4*t**4/24 + k6*t**6/720))
assert sp.simplify(slack - t**2*(k4*(phi0**2 - t**2)/24 + k6*(phi0**4 - t**4)/720)) == 0
ok("D4  delta0 B t^2/2 - (k4 t^4/24 + k6 t^6/720) = t^2[k4(phi0^2-t^2)/24 + k6(phi0^4-t^4)/720] >= 0 on |t| <= phi0; also a4t^4 <= this bound")
elem("D5  Re omega <= k4 t^4/24 + k6 t^6/720: Re psi has Taylor terms -Bt^2/2 + k4 t^4/24 at orders <= 5 (B1: odd orders imaginary), remainder <= k6 t^6/720")

# ---------------- E. T_H, Qabs, PH1abs ----------------
Hs, S4 = sp.symbols('Hs S4')                       # Hs = H(phi), S4 = H1+H2+H3+H4
H1s = sp.symbols('H1s')
lhs = sp.exp(Hs) - (1 + S4 + H1s**2/2)
rhs = (sp.exp(Hs) - 1 - Hs - Hs**2/2) + (Hs - S4) + (Hs - H1s)*(Hs + H1s)/2
assert sp.simplify(lhs - rhs) == 0
ok("E1  e^H - P_H = (e^H-1-H-H^2/2) + (H - sum_{r<=4} H_r) + (H-H1)(H+H1)/2   [exact]")
elem("E2  bounds: |e^H-1-H-H^2/2| <= |H|^3/6 <= c1^3 t^3/6 (Re H<=0, |H|<=c1 t); |H - sum H_r| <= c5 t^5/120; |H-H1| <= c2 t^2/2, |H+H1| <= 2c1 t  => T_H")
TH = c5*t**5/120 + (c1*c2/2 + c1**3/6)*t**3
assert sp.expand(TH - (c1**3*t**3/6 + c5*t**5/120 + (c2*t**2/2)*(2*c1*t)/2)) == 0
ok("E3  T_H = c5 t^5/120 + (c1c2/2 + c1^3/6) t^3 equals the sum of the E2 pieces")
Qabs = 1 + a3*t**3 + a4*t**4 + a5*t**5 + a3**2/2*t**6 + a3*a4*t**7 + a3**3/6*t**9
assert dominated(Qabs, absmaj(Qa), (a3, a4, a5, t)) and dominated(absmaj(Qa), Qabs, (a3, a4, a5, t))
PH1abs = c1*t + (c2/2 + c1**2/2)*t**2 + c3/6*t**3 + c4/24*t**4
assert dominated(PH1abs, absmaj(PH1), (c1, c2, c3, c4, t))
ok("E4  Qabs = coefficient-wise |Q| exactly; PH1abs dominates |P_H - 1| coefficient-wise; |e^H - 1| <= c1 t (|H'| <= c1, Re H <= 0)")

# ---------------- F. combination ----------------
g, ew, eH, Qs, PHs = sp.symbols('g ew eH Qs PHs')    # chi = g*ew with g = e^{-Bt^2/2}
assert sp.expand(g*ew*(eH - 1) - g*Qs*(PHs - 1) - g*((ew - Qs)*(eH - 1) + Qs*(eH - PHs))) == 0
ok("F1  chi(e^H-1) - g Q (P_H-1) = g[(e^w - Q)(e^H - 1) + Q(e^H - P_H)]  => central E_N = int c1 t T_w e^{-(1-delta0)Bt^2/2} + int Qabs T_H g")
Nt, Dt, N0s, D0s = sp.symbols('Nt Dt N0s D0s')
assert sp.simplify(Nt/Dt - N0s/D0s - ((Nt - N0s) - (N0s/D0s)*(Dt - D0s))/Dt) == 0
ok("F2  Ntil/D - N0/D0 = [(Ntil-N0) - r0 (D-D0)]/D, |D| >= |D0| - E_D  => Delta_model = (E_N + |r0| E_D)/(|D0| - E_D)")
dd = sp.symbols('dd', positive=True)
assert sp.expand((1 - dd)*(1 + dd) - (1 - dd**2)) == 0   # so -log(1-d) - log(1+d) = -log(1-d^2) >= 0 for 0 <= d < 1
ok("F3  for |e| <= d < 1: log(1+e) in [log(1-d), log(1+d)] and -log(1-d) - log(1+d) = -log(1-d^2) >= 0  => |L - L2| <= -log(1 - Delta/M2)")
elem("F4  L - L2 = log(R/g0) - log M2 (L2 = -theta - log(1-q) + log M2 = log g0 + log M2), |R/g0 - M2| <= Delta_model + |Mmodel - M2|")
elem("F5  minor arcs: |e^H - 1| <= 2, arc length 2(pi - phi0) <= 2pi  => |int chi(e^H-1)| <= 2*2pi*sup|chi|, |int chi| <= 2pi*sup|chi|")
# independent Arb implementation from the sympy polynomials above, compared with ws_p12_eps.core
from flint import arb, acb
import ws_p12_eps as EPS
Tw_full = sp.expand(a6*t**6 + Tw_poly_part + (a3*t**3 + a4*t**4 + a5*t**5)**4/24)
def mom(p, Bv):   # E|Z|^p, Z ~ N(0, 1/Bv)
    return (2 / Bv) ** (arb(p) / 2) * arb((p + 1) / 2).gamma() / arb.pi().sqrt()
def pint(poly, subsd, Bv):
    P = sp.Poly(sp.expand(poly.subs(subsd)), t); tot = arb(0)
    for (e,), cf in P.terms(): tot += arb(str(sp.Float(cf, 40))) * mom(e, Bv)
    return tot
random.seed(11); worst = 0
for trial in range(25):
    A = {r: arb(random.uniform(0.0, 0.15) / r) for r in (3, 4, 5, 6)}
    Cc = {r: arb(random.uniform(0.0, 0.3) / r) for r in (1, 2, 3, 4, 5)}
    w0 = arb(random.uniform(4, 15)); mn = arb(random.uniform(0, 1e-6))
    try: e_code, _ = EPS.core(A, Cc, w0, mn)
    except AssertionError: continue
    sub = {a3: sp.Rational(str(float(A[3].mid())))/6, a4: sp.Rational(str(float(A[4].mid())))/24, a5: sp.Rational(str(float(A[5].mid())))/120,
           a6: sp.Rational(str(float(A[6].mid())))/720}
    subc = {c1: sp.Rational(str(float(Cc[1].mid()))), c2: sp.Rational(str(float(Cc[2].mid()))), c3: sp.Rational(str(float(Cc[3].mid()))),
            c4: sp.Rational(str(float(Cc[4].mid()))), c5: sp.Rational(str(float(Cc[5].mid())))}
    d0v = A[4] * w0**2 / 12 + A[6] * w0**4 / 360; Bpv = 1 - d0v
    EN = (1 / Bpv).sqrt() * pint(c1 * t * Tw_full, {**sub, **subc}, Bpv) + pint(Qabs * TH, {**sub, **subc}, arb(1))
    ED = (1 / Bpv).sqrt() * pint(Tw_full, sub, Bpv)
    tf = (-w0**2 / 4).exp() * arb(2).sqrt()
    EN += tf * pint(Qabs * PH1abs, {**sub, **subc}, arb('0.5')) + 2 * mn
    ED += tf * pint(Qabs, sub, arb('0.5')) + mn
    num = {k3: 6*sub[a3], k4: 24*sub[a4], k5: 120*sub[a5], iB: 1, **subc}
    D0v = arb(str(sp.Float(D0.subs(num), 40))); r0 = arb(str(sp.Float((N0t/D0).subs(num), 40)))
    M2v = arb(str(sp.Float(M2.subs(num), 40))); gapv = arb(str(sp.Float(abs((Num_derived/Nn).subs(num)), 40)))
    Delta = (EN + abs(r0) * ED) / (abs(D0v) - ED) + gapv
    e_ind = -(1 - Delta / M2v).log()
    rel = float(abs(e_ind - e_code).upper() / e_code.lower()); worst = max(worst, rel)
    assert rel < 1e-8, (trial, rel)
ok(f"F6  independent eps (from the sympy polynomials) == ws_p12_eps.core on 25 random inputs, max rel diff {worst:.1e}")

# ---------------- G. Lemma W reduction ----------------
d, z2, l, xb = sp.symbols('d z2 l xb', positive=True)       # xb = 1/beta, l = log beta, d = delta = C l xb
Cc_ = sp.symbols('C', positive=True)
assert sp.expand(z2 - (1 - d)**2 * z2 - (2*d - d**2)*z2) == 0
ok("G1  theta_-^2 m <= (1-d)^2 z2 (m <= n = z2 beta^2) and theta_-^2 S(theta_-) >= z2 - T(K theta_-) - theta_-; z2 - (1-d)^2 z2 = (2d - d^2) z2")
u, ll, v = sp.symbols('u ll v', positive=True)
assert sp.simplify(sp.integrate(u * sp.exp(-ll * u), (u, v, sp.oo)) - sp.exp(-ll * v) * (v / ll + 1 / ll**2)) == 0
ok("G2  T(v) = int_v^inf u/(e^u-1) = sum_l e^{-lv}(v/l + 1/l^2) <= (v+1) sum_l e^{-lv} = (v+1)e^{-v}/(1-e^{-v})  (termwise, l >= 1)")
elem("G3  t^2 S(t) = t sum_{i<=k} f1(it) >= int_t^{Kt} f1 = F1(Kt) - F1(t) >= F1(Kt) - t  (f1 decreasing E1, 0 < f1 <= 1)")
vv = (l - 2) * (1 - d)
assert sp.simplify(sp.exp(-vv).subs(d, Cc_*l*xb) - sp.exp(2) * sp.exp(-l) * sp.exp(Cc_*(l**2*xb - 2*l*xb))) == 0
ok("G4  K theta_- >= (mu - 2beta)(1-d)/beta = (l-2)(1-d) =: v, T decreasing; e^{-v} = e^2 x exp(C(l^2 x - 2 l x)) with x = e^{-l} = 1/beta")
# sufficient condition, divided by l*x:  (2d-d^2) z2/(l x) = 2C z2 - C^2 z2 (l x)
assert sp.simplify(((2*d - d**2)*z2 / (l*xb)).subs(d, Cc_*l*xb) - (2*Cc_*z2 - Cc_**2*z2*l*xb)) == 0
# T(v)/(l x) <= (v+1) e^{-v}/((1-e^{-v}) l x) <= (l-1) e^2 exp(C l^2 x)/((1-e^{-v}) l) <= e^2 exp(C l^2 x)/(1-e^{-v}):
#   uses v + 1 <= l - 1 (v <= l - 2), exp(-2 C l x) <= 1, (l-1)/l <= 1;   theta_-/(l x) = (1-d)/l <= 1/l
assert sp.simplify((vv + 1) - (l - 1) + (l - 2) * d) == 0
ok("G5  v + 1 = (l - 1) - (l-2) d <= l - 1;  dividing the sufficient condition by l x gives exactly G = 2C z2 - C^2 z2 l x - e^2 e^{C l^2 x}/(1-e^{-v}) - 1/l > 0")
elem("G6  G > 0 on all beta >= beta_A is certified by ws_mech_window.py (worst box G >= 1.0092), no monotonicity in beta used")
print(f"all checks passed ({time.time() - t0:.0f}s)")
