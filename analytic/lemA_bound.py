"""Rigorous evaluation of the Lemma A error bound at given (n, k).

Representation (exact):  R_k/g0 - 1 = Ntil/D,
   D    = int_{-pi}^{pi} chi(phi) dphi,
   Ntil = int chi(phi) (e^{H(phi)} - 1) dphi,
chi = characteristic function of X - m (Boltzmann sum at the saddle theta),
H(phi) = i phi - log(1 - q e^{iK phi}) + log(1-q),  g0 = e^{-theta}/(1-q).

Model (Gaussian with explicit polynomials):
   omega = psi + B phi^2/2 ~ w3 + w4 + w5,  w3=-i k3 phi^3/6, w4=k4 phi^4/24, w5=i k5 phi^5/120
   Q   = 1 + w3 + w4 + w5 + w3^2/2 + w3 w4 + w3^3/6
   P_H = 1 + H1+H2+H3+H4 + H1^2/2, H1=i c1 phi, H2=-c2 phi^2/2, H3=-i c3 phi^3/6, H4=c4 phi^4/24
   D0  = int_R e^{-B phi^2/2} Q,   N0 = int_R e^{-B phi^2/2} Q (P_H - 1)
Uniform derivative bounds (all real phi):  |psi^{(r)}| <= kappa_r (r>=2),
   |H'| <= c1 = 1 + Kq/(1-q),  |H^{(r)}| <= c_r = K^r Li_{1-r}(q) (r>=2),  Re H <= 0.
Error terms:
   central arc |phi|<=phi0:   |e^omega - Q| <= Tw(|phi|) e^{(Re omega)_+},  |e^H - P_H| <= TH(|phi|)
   model tails |phi|>phi0, minor arcs via sup |chi| on [phi0, pi].
Main term compared:  M2 (closed form, see lemA_M2.py).
All quantities in Arb balls (python-flint); minor-arc sup in float64 with a stated safety margin.
"""
import math, sys
import numpy as np
from flint import arb, acb

# ---------------- ingredients ----------------
def Li_neg(s, y):
    """Li_{-s}(y) for s=0..5 in closed form (y arb)."""
    if s == 0: return y / (1 - y)
    if s == 1: return y / (1 - y) ** 2
    if s == 2: return y * (1 + y) / (1 - y) ** 3
    if s == 3: return y * (1 + 4 * y + y * y) / (1 - y) ** 4
    if s == 4: return y * (1 + 11 * y + 11 * y ** 2 + y ** 3) / (1 - y) ** 5
    if s == 5: return y * (1 + 26 * y + 66 * y ** 2 + 26 * y ** 3 + y ** 4) / (1 - y) ** 6
    raise ValueError


def S_float(t, k):
    j = np.arange(1, k + 1, dtype=float); e = np.exp(-j * t)
    return float(np.sum(j * e / (1 - e)))


def saddle_float(m, k):
    j = np.arange(1, k + 1, dtype=float)
    t = math.pi / math.sqrt(6 * m)
    for _ in range(200):
        e = np.exp(-j * t)
        g = np.sum(j * e / (1 - e)) - m
        dg = -np.sum(j * j * e / (1 - e) ** 2)
        tn = t - g / dg
        if tn <= 0: tn = t / 2
        if abs(tn - t) < 1e-15 * t: t = tn; break
        t = tn
    return t


def S_arb(t, k):
    s = arb(0)
    for i in range(1, k + 1):
        y = (-i * t).exp()
        s += i * y / (1 - y)
    return s


def ingredients(n, k):
    m = n - k; K = k + 1
    tf = saddle_float(m, k)
    lo, hi = arb(tf * (1 - 1e-11)), arb(tf * (1 + 1e-11))
    # S decreasing in theta: need S(lo) > m > S(hi)
    assert S_arb(lo, k) > m and S_arb(hi, k) < m, "saddle enclosure failed"
    th = arb(tf, tf * 1e-11 * 1.0000001)   # ball containing [lo, hi]
    kap = {r: arb(0) for r in range(2, 7)}
    for i in range(1, k + 1):
        y = (-i * th).exp()
        for r in range(2, 7):
            kap[r] += arb(i) ** r * Li_neg(r - 1, y)
    q = (-K * th).exp()
    c = {1: 1 + K * q / (1 - q)}
    for r in range(2, 6):
        c[r] = arb(K) ** r * Li_neg(r - 1, q)
    return dict(th=th, q=q, K=K, B=kap[2], k3=kap[3], k4=kap[4], k5=kap[5], k6=kap[6], c=c, tf=tf)


# ---------------- polynomials (dict power -> acb) ----------------
def pmul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            out[i + j] = out.get(i + j, acb(0)) + x * y
    return out


def padd(*ps):
    out = {}
    for p in ps:
        for i, x in p.items():
            out[i] = out.get(i, acb(0)) + x
    return out


def pscale(p, s):
    return {i: x * s for i, x in p.items()}


def gauss_int(p, B):
    """int_R e^{-B t^2/2} p(t) dt / sqrt(2 pi / B)  (odd powers vanish)."""
    tot = acb(0)
    for i, x in p.items():
        if i % 2 == 0:
            r = i // 2
            dfact = arb(1)
            for j in range(1, 2 * r, 2): dfact *= j
            tot += x * dfact / B ** r
    return tot


def abs_moment(p, Bp):
    """int_R e^{-Bp t^2/2} |t|^p dt / sqrt(2 pi / Bp)  = E|Z|^p, Z~N(0,1/Bp)."""
    # E|Z|^p = (2/Bp)^{p/2} Gamma((p+1)/2)/sqrt(pi)
    return (2 / Bp) ** (arb(p) / 2) * arb((p + 1) / 2).gamma() / arb.pi().sqrt()


def poly_abs_integral(coeffs, Bp):
    """sum_p coeffs[p] * E|Z|^p with Z ~ N(0,1/Bp); coeffs arb >= 0."""
    return sum((x * abs_moment(p, Bp) for p, x in coeffs.items()), arb(0))


def nn_mul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            out[i + j] = out.get(i + j, arb(0)) + x * y
    return out


def nn_add(*ps):
    out = {}
    for p in ps:
        for i, x in p.items():
            out[i] = out.get(i, arb(0)) + x
    return out


def nn_pow(a, e):
    out = {0: arb(1)}
    for _ in range(e): out = nn_mul(out, a)
    return out


# ---------------- minor-arc sup of |chi| ----------------
def minor_sup_log(th, k, phi0, nint=None):
    """Upper bound on max_{phi0<=phi<=pi} log|chi(phi)|.
    |chi|^2 = prod_i (1 + 4 y_i sin^2(i phi/2)/(1-y_i)^2)^{-1}, y_i = e^{-i theta}.
    On each subinterval [a,b], use for every i the minimum of sin^2(i phi/2) over [a,b]
    (zero if i[a,b]/2 contains a multiple of pi, else min at an endpoint).
    float64 with a safety factor (each log term reduced by 1e-9 relative)."""
    i = np.arange(1, k + 1, dtype=float)
    y = np.exp(-i * th)
    A = 4 * y / (1 - y) ** 2
    if nint is None:
        nint = int(min(4e5, max(2e3, 8 * k)))
    # geometric-then-uniform grid: fine near phi0
    grid = np.unique(np.concatenate([np.geomspace(phi0, 1.0, nint // 2),
                                     np.linspace(1.0, math.pi, nint // 2)]))
    worst = -np.inf
    for a, b in zip(grid[:-1], grid[1:]):
        ua, ub = i * a / 2, i * b / 2
        # contains multiple of pi?
        contains = np.floor(ub / math.pi) > np.floor(ua / math.pi)
        smin = np.minimum(np.sin(ua) ** 2, np.sin(ub) ** 2)
        # if not containing a multiple of pi, sin^2 on [ua,ub] has min at endpoints
        # unless it contains a point where sin^2 attains a smaller interior value: sin^2 is
        # unimodal between consecutive multiples of pi (max at pi/2), so endpoints are the min.
        smin = np.where(contains, 0.0, smin)
        val = -0.5 * np.sum(np.log1p(A * smin)) * (1 - 1e-9)
        worst = max(worst, val)
    return worst


# ---------------- main bound ----------------
def lemma_A(n, k, s_phi0=0.4, verbose=True):
    I = ingredients(n, k)
    th, q, K, B, k3, k4, k5, k6, c = I['th'], I['q'], I['K'], I['B'], I['k3'], I['k4'], I['k5'], I['k6'], I['c']
    ii = acb(0, 1)
    # model polynomials
    w3 = {3: -ii * k3 / 6}; w4 = {4: acb(k4 / 24)}; w5 = {5: ii * k5 / 120}
    Q = padd({0: acb(1)}, w3, w4, w5, pscale(pmul(w3, w3), acb(0.5)), pmul(w3, w4),
             pscale(pmul(pmul(w3, w3), w3), acb(1) / 6))
    H1 = {1: ii * c[1]}; H2 = {2: acb(-c[2] / 2)}; H3 = {3: -ii * c[3] / 6}; H4 = {4: acb(c[4] / 24)}
    PH1 = padd(H1, H2, H3, H4, pscale(pmul(H1, H1), acb(0.5)))      # P_H - 1
    D0 = gauss_int(Q, B)               # normalised by sqrt(2pi/B)
    N0 = gauss_int(pmul(Q, PH1), B)
    Mmodel = 1 + N0 / D0
    # closed form M2
    h1 = -c[1]; h2 = c[2]; h3 = -c[3]; h4 = c[4]
    Nn = 1 + k4 / (8 * B * B) - 5 * k3 * k3 / (24 * B ** 3)
    E1 = (k3 / (2 * B * B) - k5 / (8 * B ** 3) + 105 * k3 * k4 / (144 * B ** 4)
          - 945 * k3 ** 3 / (1296 * B ** 5)) / Nn
    E2 = (1 + k4 / (2 * B * B) - 5 * k3 * k3 / (4 * B ** 3)) / B
    M2 = 1 - h1 * E1 - (h2 + h1 * h1) * E2 / 2 + 5 * h3 * k3 / (12 * B ** 3) + h4 / (8 * B * B)
    # ---- error terms (all normalised by sqrt(2 pi / B)) ----
    phi0 = arb(s_phi0) * th
    delta0 = (k4 * phi0 ** 2 / 12 + k6 * phi0 ** 4 / 360) / B
    assert delta0 < 1, "delta0 >= 1"
    Bp = (1 - delta0) * B
    norm_fac = (B / Bp).sqrt()          # converts N(0,1/Bp) moments to the sqrt(2pi/B) normalisation
    a3, a4, a5, a6 = k3 / 6, k4 / 24, k5 / 120, k6 / 720
    om_hat = {3: a3, 4: a4, 5: a5}                 # |omega_hat| <= sum a_r t^r
    # |e^w - Q| <= [ a6 t^6 + |w4^2/2 + w5^2/2 + w3 w5 + w4 w5| + |w^3 - w3^3|/6 + |w|^4/24 ] e^{(Re)+}
    Tw = nn_add({6: a6},
                {8: a4 * a4 / 2 + a3 * a5, 10: a5 * a5 / 2, 9: a4 * a5},
                nn_mul(nn_mul({4: a4, 5: a5}, nn_pow(om_hat, 2)), {0: arb(3) / 6}),
                nn_mul(nn_pow(om_hat, 4), {0: arb(1) / 24}))
    TH = {5: c[5] / 120, 3: c[1] * c[2] / 2 + c[1] ** 3 / 6}
    # |Q| <= Qabs(t)
    Qabs = {0: arb(1), 3: a3, 4: a4, 5: a5, 6: a3 * a3 / 2, 7: a3 * a4, 9: a3 ** 3 / 6}
    # central-arc errors
    eN_c = norm_fac * poly_abs_integral(nn_mul(Tw, {1: c[1]}), Bp) + poly_abs_integral(nn_mul(Qabs, TH), B)
    eD_c = norm_fac * poly_abs_integral(Tw, Bp)
    # model tails beyond phi0: int_{|t|>phi0} e^{-Bt^2/2} P(t) <= e^{-B phi0^2/4} int e^{-Bt^2/4} P(t)
    PH1abs = {1: c[1], 2: c[2] / 2, 3: c[3] / 6, 4: c[4] / 24}
    PH1abs[2] = PH1abs[2] + c[1] ** 2 / 2
    tailfac = (-B * phi0 ** 2 / 4).exp() * arb(2).sqrt()   # sqrt(2) converts B/2-normalisation
    eN_t = tailfac * poly_abs_integral(nn_mul(Qabs, PH1abs), B / 2)
    eD_t = tailfac * poly_abs_integral(Qabs, B / 2)
    # minor arcs: |int| <= 2pi * sup|chi| * {2 for N, 1 for D}; normalise by sqrt(2pi/B)
    lsup = minor_sup_log(float(th.mid()), k, float(phi0.mid()) * (1 - 1e-9))
    sup = arb(lsup).exp()
    nrm = (2 * arb.pi() / B).sqrt()
    eN_m = 2 * 2 * arb.pi() * sup / nrm
    eD_m = 2 * arb.pi() * sup / nrm
    EN = eN_c + eN_t + eN_m
    ED = eD_c + eD_t + eD_m
    # combine:  | Ntil/D - N0/D0 | <= (EN + |N0/D0| ED) / (|D0| - ED)
    r0 = N0 / D0
    Delta_model = (EN + abs(r0) * ED) / (abs(D0) - ED)
    Delta_M2 = abs(Mmodel - M2)
    Delta = Delta_model + Delta_M2
    M2r = M2
    eps = -(1 - Delta / M2r).log()          # |log(R/g0) - log M2| <= eps
    out = dict(n=n, k=k, eps=eps, Delta_model=Delta_model, Delta_M2=Delta_M2, eN_c=eN_c, eD_c=eD_c,
               eN_t=eN_t, eN_m=eN_m, lsup=lsup, delta0=delta0, M2=M2, Mmodel=Mmodel, th=th, q=q, K=K)
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]); ks = [int(a) for a in sys.argv[2:]]
    U = lambda x: float(x.upper())
    for k in ks:
        o = lemma_A(n, k)
        print(f"n={n} k={k}: eps <= {U(o['eps']):.3e}  [model {U(o['Delta_model']):.2e}, "
              f"M2-gap {U(o['Delta_M2']):.2e}; centralN {U(o['eN_c']):.2e} centralD {U(o['eD_c']):.2e}"
              f" tails {U(o['eN_t']):.1e} minor-log-sup {o['lsup']:.1f} delta0 {float(o['delta0'].mid()):.3f}]")
