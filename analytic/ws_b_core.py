"""Lemma B, shared core.  M2 as an explicit function of the six base variables
    v = (theta, nu, b, s3, s4, s5),      q = e^{-nu},  nu = K theta,  b = s2,
    A_r = s_r b^{-r/2} theta^{r/2-1},  C_r = g_r(nu) b^{-r/2} theta^{r/2} (r=2..4),  g_r = nu^r Li_{1-r}(e^{-nu}),
    C_1 = c1 theta^{3/2} b^{-1/2},  c1 = 1 + (nu/theta) q/(1-q),
    Nn = 1 + A4/8 - 5A3^2/24,  E1 = (A3/2 - A5/8 + 105A3A4/144 - 945A3^3/1296)/Nn,  E2 = 1 + A4/2 - 5A3^2/4,
    M2 = 1 + C1 E1 - (C2 + C1^2) E2/2 - 5 C3 A3/12 + C4/8          (= LemmaA.md sec. 3, B-scaled; identical to ws_p12_eps.core)
Forward-mode dual numbers over arb give rigorous enclosures of grad M2 over a box (all six variables independent)."""
from flint import arb
from ws_p12_common import F_fun, P as Eul, f_fun
from ws_p12_eps import TV
VARS = ('th', 'nu', 'b', 's3', 's4', 's5')
NV = len(VARS)
z2 = arb(2).zeta()
TAU = arb('0.01')

class D:
    __slots__ = ('v', 'g')
    def __init__(s, v, g=None):
        s.v = arb(v); s.g = g if g is not None else [arb(0)] * NV
    @staticmethod
    def var(val, i):
        g = [arb(0)] * NV; g[i] = arb(1); return D(val, g)
    def _c(o): return o if isinstance(o, D) else D(o)
    def __add__(s, o): o = D._c(o); return D(s.v + o.v, [a + b for a, b in zip(s.g, o.g)])
    __radd__ = __add__
    def __neg__(s): return D(-s.v, [-a for a in s.g])
    def __sub__(s, o): return s + (-D._c(o))
    def __rsub__(s, o): return D._c(o) - s
    def __mul__(s, o):
        o = D._c(o); return D(s.v * o.v, [a * o.v + s.v * b for a, b in zip(s.g, o.g)])
    __rmul__ = __mul__
    def recip(s):
        iv = 1 / s.v; return D(iv, [-a * iv * iv for a in s.g])
    def __truediv__(s, o): return s * D._c(o).recip()
    def __rtruediv__(s, o): return D._c(o) * s.recip()
    def powr(s, a):            # real exponent a (arb), s.v > 0
        a = arb(a)
        if a == 0: return D(1)
        pv = s.v ** a; d = a * s.v ** (a - 1); return D(pv, [x * d for x in s.g])
    def exp(s):
        e = s.v.exp(); return D(e, [x * e for x in s.g])

def Eul_gen(r, y):
    from ws_p12_common import EUL
    out = D(0) if isinstance(y, D) else arb(0)
    for c in reversed(EUL[r]): out = out * y + c
    return out

def M2_formula(th, nu, b, s3, s4, s5):
    """Works for D or arb inputs.  Returns (M2, dict of intermediates)."""
    pw = (lambda x, a: x.powr(a)) if isinstance(th, D) else (lambda x, a: x ** arb(a) if arb(a) != 0 else arb(1))
    ex = (lambda x: x.exp())
    q = ex(-nu)
    s = {3: s3, 4: s4, 5: s5}
    A = {r: s[r] * pw(b, -arb(r) / 2) * pw(th, arb(r) / 2 - 1) for r in (3, 4, 5)}
    C = {}
    for r in (2, 3, 4):
        Li = q * Eul_gen(r - 1, q) / pw(1 - q, r)          # Li_{1-r}(q)
        C[r] = pw(nu, r) * Li * pw(b, -arb(r) / 2) * pw(th, arb(r) / 2)
    c1 = 1 + nu / th * q / (1 - q)
    C[1] = c1 * pw(th, arb(3) / 2) * pw(b, -arb(1) / 2)
    A3, A4, A5 = A[3], A[4], A[5]
    Nn = 1 + A4 / 8 - 5 * A3 * A3 / 24
    E1 = (A3 / 2 - A5 / 8 + 105 * A3 * A4 / 144 - 945 * A3 * A3 * A3 / 1296) / Nn
    E2 = 1 + A4 / 2 - 5 * A3 * A3 / 4
    M2 = 1 + C[1] * E1 - (C[2] + C[1] * C[1]) * E2 / 2 - 5 * C[3] * A3 / 12 + C[4] / 8
    return M2, dict(A=A, C=C, q=q, c1=c1)

def smax(r, tmax):
    """sup of s_r(k,t) = t sum_{i<=k} f_r(it) for t <= tmax (f_r >= 0, int_0^inf f_r = r! zeta(2))."""
    return arb.fac_ui(r) * z2 + tmax * TV[r]

def deltas_box(th, nup):
    """th, nup: arb balls (a box).  Returns rigorous upper bounds of the k -> k+1 increments, valid for every
    (n,k) whose (theta_k, nu'_k) lies in the box.  See LemmaB.md sec. 2."""
    nu = nup - th.log()
    q = th * (-nup).exp()
    c1 = 1 + nu * (-nup).exp() / (1 - q)
    bmin = F_fun(2, nu - th) - th * (1 + TAU) * TV[2]           # lower bound of b(k, t), t in [theta, theta(1+tau)]
    assert bmin > 0
    # Lemma B.1: theta' <= theta(1+tau)  <=  tau*bmin >= (1+tau)^3 c1 theta^2
    assert TAU * bmin > (1 + TAU) ** 3 * c1 * th ** 2, "tau-step check failed"
    dth = (1 + TAU) ** 3 * c1 * th ** 3 / bmin                   # Lemma B.2
    dnu = th * (1 + TAU) + (nu / th + 1) * dth                    # Lemma B.3
    dnu_hi = arb(dnu.upper())
    tm = th * (1 + TAU)
    ds = {}
    for r in (2, 3, 4, 5):
        nseg = arb(nu.lower()).union(arb((nu + dnu_hi).upper()))
        ds[r] = dth * ((1 + r) * smax(r, tm) + smax(r + 1, tm)) / th + tm * f_fun(r, nseg)   # Lemma B.4
    return dict(nu=nu, q=q, c1=c1, bmin=bmin, dth=dth, dnu=dnu, ds=ds)

def hull_point(th, nup, dl):
    """Box containing both base points (theta_k, nu_k, b_k, s_r,k) and (theta_{k+1}, ...)."""
    nu = dl['nu']
    th_h = arb(th.lower()).union(arb((th * (1 + TAU)).upper()))
    nu_h = arb(nu.lower()).union(arb((nu + dl['dnu']).upper()))
    lo = nu - th; hi = nu + dl['dnu']
    tm = th * (1 + TAU)
    s = {}
    for r in (2, 3, 4, 5):
        Fl = F_fun(r, arb(lo.lower())); Fh = F_fun(r, arb(hi.upper()))
        w = arb(0, 1) * tm * TV[r]
        s[r] = (Fl + w).union(Fh + w)
    return th_h, nu_h, s

def dlogM2_box(th, nup):
    """Rigorous upper bound of |log M2_{k+1} - log M2_k| over the box (mean value theorem on the segment)."""
    dl = deltas_box(th, nup)
    th_h, nu_h, s = hull_point(th, nup, dl)
    xs = [D.var(th_h, 0), D.var(nu_h, 1), D.var(s[2], 2), D.var(s[3], 3), D.var(s[4], 4), D.var(s[5], 5)]
    M2, _ = M2_formula(*xs)
    assert M2.v > 0
    grad = [abs(gi / M2.v) for gi in M2.g]
    dv = [dl['dth'], dl['dnu'], dl['ds'][2], dl['ds'][3], dl['ds'][4], dl['ds'][5]]
    tot = sum((arb(g.upper()) * arb(abs(d).upper()) for g, d in zip(grad, dv)), arb(0))
    return tot, dl, [arb(g.upper()) * arb(abs(d).upper()) for g, d in zip(grad, dv)]

def L2_over_theta_box(th, nup):
    """Enclosure of L2/theta over the box, using P1.3 for s_r (k fixed)."""
    nu = nup - th.log(); kth = nu - th
    s = {r: F_fun(r, kth) + arb(0, 1) * th * TV[r] for r in (2, 3, 4, 5)}
    M2, _ = M2_formula(th, nu, s[2], s[3], s[4], s[5])
    assert M2 > 0
    q = th * (-nup).exp()
    return -1 + (-(1 - q).log()) / th + M2.log() / th, q

def T2_over_thq(th, nup):
    """Lower bound of log(1 + q(1-e^{-theta})/(1-q)) / (theta q)."""
    q = th * (-nup).exp()
    x = q * (-(-th).expm1()) / (1 - q)
    return (1 + x).log() / (th * q)

def load_gamma(paths=('ws_p12_regionI.out', 'ws_p12_regionII.out')):
    import re
    cells = []
    for p in paths:
        for line in open(p):
            m = re.match(r"nu' in \[([-+0-9.]+),([-+0-9.]+)\]\s+Gamma <= ([0-9.infa]+)", line)
            if m: cells.append((float(m.group(1)), float(m.group(2)), float(m.group(3))))
    return cells

def gamma_on(cells, a, b):
    """max Gamma over all cells meeting [a, b] (both regions)."""
    vals = [g for (x, y, g) in cells if y >= a - 1e-12 and x <= b + 1e-12]
    assert vals
    return max(vals)
