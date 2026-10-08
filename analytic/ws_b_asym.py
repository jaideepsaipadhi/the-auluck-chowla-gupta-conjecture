"""Lemma B, Region II (theta <= TH1): graded-monomial majorants for
   (i)  |log M2_{k+1} - log M2_k| / (theta q)        -> Sigma >= T2low - that
   (ii) |log M2_k| / theta                           -> sign of L2 at the window edges (B3).
Usage: python3 ws_b_asym.py TH1 NUP_LO NUP_HI [h]
Monomials t^p q^bq nu^j (t = sqrt(theta_k), q = q_k, nu = nu_k; q = t^2 e^{-nu'}), coefficients are arb upper bounds
valid for theta <= TH1, nu' in the cell, at BOTH points k and k+1 (LemmaB.md sec. 4).  Every majorant has nonneg.
coefficients; after division by t^2 q (resp. t^2), each monomial t^P nu^j e^{..nu'} is nondecreasing in t on
(0, sqrt TH1] (asserted: P > 0 and nu_min > 2j/P), so its sup is its value at t1 = sqrt TH1, nu = v1 + log(1/TH1)."""
import sys
from flint import arb
from ws_p12_common import F_fun, P as Eul
from ws_b_core import smax, TAU
from ws_p12_eps import TV

class GP:
    def __init__(s, d=None): s.d = {k: v for k, v in (d or {}).items()}
    @staticmethod
    def m(c, p=0, bq=0, j=0): return GP({(p, bq, j): arb(c)})
    def __add__(s, o):
        if not isinstance(o, GP): o = GP.m(o)
        r = dict(s.d)
        for k, v in o.d.items(): r[k] = r.get(k, arb(0)) + v
        return GP(r)
    __radd__ = __add__
    def __mul__(s, o):
        if not isinstance(o, GP): return GP({k: v * arb(o) for k, v in s.d.items()})
        r = {}
        for k1, v1 in s.d.items():
            for k2, v2 in o.d.items():
                k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2]); r[k] = r.get(k, arb(0)) + v1 * v2
        return GP(r)
    __rmul__ = __mul__

def setup(TH1, v0, v1):
    t1 = TH1.sqrt(); L1 = (1 / TH1).log()
    nmin = arb(v0) + L1                  # nu = nu' + log(1/theta) >= v0 + log(1/TH1)
    assert nmin > 3
    tm = TH1 * (1 + TAU)
    S = {r: smax(r, tm) for r in range(2, 7)}
    bmin = F_fun(2, nmin - TH1) - tm * TV[2]; assert bmin > 1
    ymax = (-nmin).exp(); q1 = TH1 * (-arb(v0)).exp()
    # sup of c1*theta^2 (increasing in theta for theta <= TH1 < e^{-1/2}):
    c1t2 = (1 + (arb(v1) + L1) * (-arb(v0)).exp() / (1 - q1)) * TH1 ** 2
    assert TAU * bmin > (1 + TAU) ** 3 * c1t2          # Lemma B.1 (tau-step)
    rho = 1 + TH1 * (1 + TAU) / nmin + (1 + TAU) ** 3 * c1t2 * (1 + TH1 / nmin) / bmin   # nu_{k+1} <= rho nu_k
    return dict(t1=t1, L1=L1, nmin=nmin, S=S, bmin=bmin, ymax=ymax, q1=q1, rho=rho, v0=arb(v0), v1=arb(v1))

def majorants(P):
    S, bm, y, q1, rho, T1 = P['S'], P['bmin'], P['ymax'], P['q1'], P['rho'], 1 + TAU
    one = GP.m(1)
    th_p = lambda e: GP.m(T1 ** max(e, 0), 2 * e) if True else None   # sup of theta^e over both points (theta'<=(1+tau)theta, theta'>=theta)
    def thp(e):   # sup over [theta, (1+tau)theta] of theta^e, as GP in t (e may be half-integer)
        e = arb(e)
        return GP.m(T1 ** e if e > 0 else 1, int(round(float(2 * e.mid()))))
    c1k = one + GP.m(1 / (1 - q1), -2, 1, 1)                       # c1 at point k
    c1s = one + GP.m(rho / (1 - q1), -2, 1, 1)                     # sup over both points
    dth = c1k * GP.m(T1 ** 3 / bm, 6)                             # Lemma B.2
    dnu = GP.m(T1, 2) + (GP.m(1, -2, 0, 1) + one) * dth           # Lemma B.3
    Eu = lambda r: Eul(r - 1, y) / (1 - y) ** r
    fsup = lambda r: GP.m(rho ** r * Eu(r), 0, 1, r)               # sup of f_r on [nu, nu_{k+1}] (u^r e^{-u}... <= (rho nu)^r q Eu)
    ds = {r: dth * GP.m((1 + r) * S[r] + S[r + 1], -2) + GP.m(T1, 2) * fsup(r) for r in (2, 3, 4, 5)}  # Lemma B.4
    db = ds[2]
    # A_r = s_r b^{-r/2} theta^{r/2-1}
    sA, dA = {}, {}
    for r in (3, 4, 5):
        e = arb(r) / 2 - 1
        sA[r] = GP.m(S[r] * bm ** (-arb(r) / 2)) * thp(e)
        dA[r] = (ds[r] * GP.m(bm ** (-arb(r) / 2)) * thp(e)
                 + GP.m(S[r] * (arb(r) / 2) * bm ** (-arb(r) / 2 - 1)) * db * thp(e)
                 + GP.m(S[r] * bm ** (-arb(r) / 2) * e) * thp(e - 1) * dth)
    # g_r(nu) = nu^r Li_{1-r}(e^{-nu});  |g_r'| <= (r g_r + g_{r+1})/nu
    gsup = {r: fsup(r) for r in range(2, 6)}
    dg = {r: dnu * (GP.m(r) * gsup[r] + gsup[r + 1]) * GP.m(1, 0, 0, -1) for r in (2, 3, 4)}
    sC, dC = {}, {}
    for r in (2, 3, 4):
        hr = arb(r) / 2
        sC[r] = gsup[r] * GP.m(bm ** (-hr)) * thp(hr)
        dC[r] = (dg[r] * GP.m(bm ** (-hr)) * thp(hr) + gsup[r] * GP.m(hr * bm ** (-hr - 1)) * db * thp(hr)
                 + gsup[r] * GP.m(bm ** (-hr) * hr) * thp(hr - 1) * dth)
    # c1 = 1 + h(nu)/theta, h(u) = u e^{-u}/(1-e^{-u}),  |h'(u)| <= (u+1) e^{-u}/(1-e^{-u})^2 <= (rho nu+1) q/(1-y)^2 on [nu, nu_{k+1}]
    dc1 = dnu * (GP.m(rho, 0, 1, 1) + GP.m(1, 0, 1, 0)) * GP.m(1 / (1 - y) ** 2, -2) + GP.m(rho / (1 - q1), 0, 1, 1) * dth * GP.m(1, -4)
    sC[1] = c1s * GP.m(bm ** arb(-0.5)) * thp(arb(3) / 2)
    dC[1] = (dc1 * GP.m(bm ** arb(-0.5)) * thp(arb(3) / 2) + c1s * GP.m(bm ** arb(-1.5) / 2) * db * thp(arb(3) / 2)
             + c1s * GP.m(bm ** arb(-0.5) * arb(3) / 2) * thp(arb(1) / 2) * dth)
    return sA, dA, sC, dC

def maxval(g, P, shift_t, shift_q):
    """sup over theta<=TH1, nu' in cell of g / (t^{shift_t} q^{shift_q})."""
    tot = arb(0)
    for (p, bq, j), c in g.d.items():
        Pp = p + 2 * bq - shift_t - 2 * shift_q; be = bq - shift_q
        assert (Pp == 0 and j == 0) or (Pp > 0 and P['nmin'] > arb(2 * j) / Pp), ("non-monotone", p, bq, j)
        ex = (-P['v0'] * be).exp() if be >= 0 else (P['v1'] * (-be)).exp()
        nu_hi = P['v1'] + P['L1']
        tot += c * P['t1'] ** Pp * (nu_hi ** j if j >= 0 else P['nmin'] ** j) * ex
    return tot

def cell(TH1, v0, v1):
    P = setup(TH1, v0, v1)
    sA, dA, sC, dC = majorants(P)
    M = lambda g: maxval(g, P, 0, 0)
    # scalar sups (theta-monotone majorants evaluated at TH1)
    A3, A4, A5 = M(sA[3]), M(sA[4]), M(sA[5])
    Nnmin = 1 - arb(5) / 24 * A3 ** 2; assert Nnmin > 0
    num = sA[3] * GP.m(arb(1) / 2) + sA[5] * GP.m(arb(1) / 8) + sA[3] * sA[4] * GP.m(arb(105) / 144) + sA[3] * sA[3] * sA[3] * GP.m(arb(945) / 1296)
    dnum = dA[3] * GP.m(arb(1) / 2) + dA[5] * GP.m(arb(1) / 8) + (dA[3] * sA[4] + sA[3] * dA[4]) * GP.m(arb(105) / 144) \
        + dA[3] * sA[3] * sA[3] * GP.m(3 * arb(945) / 1296)
    dNn = dA[4] * GP.m(arb(1) / 8) + dA[3] * sA[3] * GP.m(arb(10) / 24)
    # E1 = num/Nn :  |dE1| <= dnum/Nnmin + num * dNn / Nnmin^2 ;  |E1| <= num/Nnmin
    sE1 = num * GP.m(1 / Nnmin)
    dE1 = dnum * GP.m(1 / Nnmin) + num * dNn * GP.m(1 / Nnmin ** 2)
    sE2 = GP.m(1) + sA[4] * GP.m(arb(1) / 2) + sA[3] * sA[3] * GP.m(arb(5) / 4)
    dE2 = dA[4] * GP.m(arb(1) / 2) + dA[3] * sA[3] * GP.m(arb(10) / 4)
    C1, C2, C3, C4 = sC[1], sC[2], sC[3], sC[4]
    # M2 - 1 = C1 E1 - (C2 + C1^2) E2/2 - 5 C3 A3/12 + C4/8
    sM = C1 * sE1 + (C2 + C1 * C1) * sE2 * GP.m(arb(1) / 2) + C3 * sA[3] * GP.m(arb(5) / 12) + C4 * GP.m(arb(1) / 8)
    dM = (dC[1] * sE1 + C1 * dE1 + (dC[2] + dC[1] * C1 * GP.m(2)) * sE2 * GP.m(arb(1) / 2) + (C2 + C1 * C1) * dE2 * GP.m(arb(1) / 2)
          + (dC[3] * sA[3] + C3 * dA[3]) * GP.m(arb(5) / 12) + dC[4] * GP.m(arb(1) / 8))
    Mdev = M(sM); assert Mdev < 1
    M2min = 1 - Mdev
    dlog_over_thq = maxval(dM, P, 2, 1) / M2min                    # |Delta log M2| / (theta q)
    logM_over_th = maxval(sM, P, 2, 0) / M2min                    # |log M2| / theta  (|log(1+x)| <= |x|/(1-|x|))
    T2low = (1 - TH1 / 2) * (1 - TH1 * P['q1'] / (2 * (1 - P['q1'])))   # x>=theta q(1-theta/2), x<=theta q/(1-q), log(1+x)>=x(1-x/2)
    sigma = T2low - dlog_over_thq
    return sigma, logM_over_th, P

if __name__ == "__main__":
    TH1 = arb(sys.argv[1]); VLO = float(sys.argv[2]); VHI = float(sys.argv[3])
    h = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
    v = VLO; worst = 1e9; wl = 0
    while v < VHI - 1e-12:
        w = min(v + h, VHI)
        sig, lm, P = cell(TH1, v, w)
        worst = min(worst, float(sig.lower())); wl = max(wl, float(lm.upper()))
        print(f"nu' in [{v:+.3f},{w:+.3f}]  Sigma >= {float(sig.lower()):.5f}   |log M2|/theta <= {float(lm.upper()):.5f}", flush=True)
        v = w
    print(f"Region II (theta <= {sys.argv[1]}): Sigma >= {worst:.5f};  sup |log M2|/theta <= {wl:.5f}")
