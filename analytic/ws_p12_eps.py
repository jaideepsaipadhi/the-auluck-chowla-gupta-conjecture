"""Uniform Lemma A: rigorous E(theta, nu') >= eps over boxes, in scaled (dimensionless) form.
Variables: w = phi*sqrt(B).  A_r := kappa_r B^{-r/2} = s_r b^{-r/2} theta^{r/2-1},  C_r := c_r B^{-r/2}
= g_r b^{-r/2} theta^{r/2} (r>=2), C_1 := c1 B^{-1/2},  w0 := phi0 sqrt(B) = 0.4 sqrt(b/theta).
Scaled ingredients (P1):  s_r := theta^{r+1} kappa_r in F_r(nu-theta) +- theta*TV_r  (b = s_2),
   g_r := theta^r c_r = nu^r Li_{1-r}(e^{-nu}) (exact),  c1 = 1 + nu e^{-nu'}/(1-q),  q = theta e^{-nu'}.
Minor arcs (P2): sup|chi| <= exp(-c_*/theta).
core() is lemma_A of lemA_bound.py with B := 1 (identical algebra; every term is dimensionless)."""
import sys, math
from flint import arb, acb
from ws_p12_common import F_fun, Li_neg
import lemA_bound as LB
TV = {1: arb('1.000001'), 2: arb('1.000001'), 3: arb('2.000001'), 4: arb('6.086765'),
      5: arb('25.691943'), 6: arb('138.425')}               # certified by ws_mech_tv.py (rounded up; enclosure upper ends
# 1+1.5e-12, 1+3.4e-12, 2+2.7e-11, 6.08676444, 25.6919422, 138.424985).  The former constants came from ws_p12_tv.py, whose
# charge g(b)^r h(a)-g(a)^r h(b) on sign-undetermined cells is NOT a valid variation bound (REFEREE A.1); ws_p12_tv.py is retired.

def core(A, C, w0, minor_norm):
    """A: {3..6}, C: {1..5} (arb), w0 (arb), minor_norm = 2pi*sup|chi|/sqrt(2pi/B) (arb).
    Returns (eps, pieces)."""
    one = arb(1); ii = acb(0, 1)
    k3, k4, k5, k6 = A[3], A[4], A[5], A[6]; c = C
    w3 = {3: -ii * k3 / 6}; w4 = {4: acb(k4 / 24)}; w5 = {5: ii * k5 / 120}
    Q = LB.padd({0: acb(1)}, w3, w4, w5, LB.pscale(LB.pmul(w3, w3), acb(0.5)), LB.pmul(w3, w4),
                LB.pscale(LB.pmul(LB.pmul(w3, w3), w3), acb(1) / 6))
    H1 = {1: ii * c[1]}; H2 = {2: acb(-c[2] / 2)}; H3 = {3: -ii * c[3] / 6}; H4 = {4: acb(c[4] / 24)}
    PH1 = LB.padd(H1, H2, H3, H4, LB.pscale(LB.pmul(H1, H1), acb(0.5)))
    D0 = LB.gauss_int(Q, one); N0 = LB.gauss_int(LB.pmul(Q, PH1), one)
    Mmodel = 1 + N0 / D0
    h1 = -c[1]; h2 = c[2]; h3 = -c[3]; h4 = c[4]; B = one
    Nn = 1 + k4 / 8 - 5 * k3 * k3 / 24
    E1 = (k3 / 2 - k5 / 8 + 105 * k3 * k4 / 144 - 945 * k3 ** 3 / 1296) / Nn
    E2 = (1 + k4 / 2 - 5 * k3 * k3 / 4)
    M2 = 1 - h1 * E1 - (h2 + h1 * h1) * E2 / 2 + 5 * h3 * k3 / 12 + h4 / 8
    delta0 = k4 * w0 ** 2 / 12 + k6 * w0 ** 4 / 360
    assert delta0 < 1
    Bp = 1 - delta0; norm_fac = (1 / Bp).sqrt()
    a3, a4, a5, a6 = k3 / 6, k4 / 24, k5 / 120, k6 / 720
    om = {3: a3, 4: a4, 5: a5}
    Tw = LB.nn_add({6: a6}, {8: a4 * a4 / 2 + a3 * a5, 10: a5 * a5 / 2, 9: a4 * a5},
                   LB.nn_mul(LB.nn_mul({4: a4, 5: a5}, LB.nn_pow(om, 2)), {0: arb(3) / 6}),
                   LB.nn_mul(LB.nn_pow(om, 4), {0: arb(1) / 24}))
    TH = {5: c[5] / 120, 3: c[1] * c[2] / 2 + c[1] ** 3 / 6}
    Qabs = {0: one, 3: a3, 4: a4, 5: a5, 6: a3 * a3 / 2, 7: a3 * a4, 9: a3 ** 3 / 6}
    eN_c = norm_fac * LB.poly_abs_integral(LB.nn_mul(Tw, {1: c[1]}), Bp) + LB.poly_abs_integral(LB.nn_mul(Qabs, TH), one)
    eD_c = norm_fac * LB.poly_abs_integral(Tw, Bp)
    PH1abs = {1: c[1], 2: c[2] / 2 + c[1] ** 2 / 2, 3: c[3] / 6, 4: c[4] / 24}
    tailfac = (-w0 ** 2 / 4).exp() * arb(2).sqrt()
    eN_t = tailfac * LB.poly_abs_integral(LB.nn_mul(Qabs, PH1abs), one / 2)
    eD_t = tailfac * LB.poly_abs_integral(Qabs, one / 2)
    EN = eN_c + eN_t + 2 * minor_norm
    ED = eD_c + eD_t + minor_norm
    r0 = N0 / D0
    Dm = (abs(D0) - ED)
    assert Dm > 0
    # exact identity (ws_p12_sym.py, B=1):  Mmodel - M2 = Num/Nn, Num = 11 monomials
    Num = (-arb(25)/48*c[4]*k3**2 + c[4]*k4/6 + (c[2] + c[1]**2)*(k4**2/32 + arb(25)/192*k3**4 - arb(25)/192*k3**2*k4)
           + arb(5)/4*c[3]*k3**3 + arb(7)/48*c[3]*k5 - arb(25)/24*c[3]*k3*k4)
    gap = abs(Num) / Nn
    Delta = (EN + abs(r0) * ED) / Dm + gap
    assert Delta < M2.real
    eps = -(1 - Delta / M2.real).log()
    return eps, dict(M2=M2.real, Delta=Delta, eN_c=eN_c, gap=gap)

def scaled_from_th_nup(th, nup, cstar):
    """th, nup: arb balls (box).  Rigorous enclosures of A, C, w0, minor_norm over the box."""
    nu = nup - th.log()
    kth = nu - th
    s = {r: F_fun(r, kth) + arb(0, 1) * th * TV[r] for r in range(2, 7)}
    b = s[2]
    q = th * (-nup).exp()
    g = {r: nu ** r * Li_neg(r - 1, (-nu).exp()) for r in range(2, 6)}
    A = {r: s[r] * b ** (-arb(r) / 2) * th ** (arb(r) / 2 - 1) for r in range(3, 7)}
    C = {r: g[r] * b ** (-arb(r) / 2) * th ** (arb(r) / 2) for r in range(2, 6)}
    c1 = 1 + nu * (-nup).exp() / (1 - q)
    C[1] = c1 * (th ** 3 / b).sqrt()
    w0 = arb('0.4') * (b / th).sqrt()
    minor_norm = 2 * arb.pi() * (-cstar / th).exp() * (b / (2 * arb.pi() * th ** 3)).sqrt()
    return A, C, w0, minor_norm, q

def scaled_from_nk(n, k, cstar=None):
    """Exact-ingredient version (pointwise sanity check against lemA_bound.lemma_A)."""
    I = LB.ingredients(n, k); B = I['B']
    A = {3: I['k3'] / B ** 1.5, 4: I['k4'] / B ** 2, 5: I['k5'] / B ** 2.5, 6: I['k6'] / B ** 3}
    C = {r: I['c'][r] / B ** (arb(r) / 2) for r in range(1, 6)}
    w0 = arb('0.4') * I['th'] * B.sqrt()
    return A, C, w0, I

if __name__ == "__main__":
    # sanity: same eps as lemA_bound (with its numeric minor sup) at a pointwise (n,k)
    n, k = int(sys.argv[1]), int(sys.argv[2])
    A, C, w0, I = scaled_from_nk(n, k)
    o = LB.lemma_A(n, k)
    mn = 2 * arb.pi() * arb(o['lsup']).exp() / (2 * arb.pi() / I['B']).sqrt()
    e, _ = core(A, C, w0, mn)
    print("scaled core eps:", e.str(6), "  lemA_bound eps:", o['eps'].str(6))
    th = I['th']; nup = I['K'] * th + th.log()
    A2, C2, w02, mn2, q = scaled_from_th_nup(arb(th.mid()), arb(nup.mid()), arb('0.159'))
    e2, d = core(A2, C2, w02, mn2)
    print("P1-ingredient eps (point box):", e2.str(6), " ratio eps/(theta q) =", (e2 / (th * q)).str(4))
