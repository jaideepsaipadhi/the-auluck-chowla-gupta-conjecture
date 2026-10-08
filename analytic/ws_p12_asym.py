"""Region II (theta <= TH1): rigorous majorant of Gamma := eps/(theta q) via a graded-monomial argument.
Usage: python3 ws_p12_asym.py TH1 NUP_LO NUP_HI CSTAR [h]
Scaled quantities (B=1 normalisation, t := sqrt(theta), q = t^2 e^{-nu'}, nu = nu' + 2 log(1/t)):
   A_r = a_r t^{r-2}          (r=3..6),  a_r = s_r b^{-r/2} <= smax_r bmin^{-r/2}
   C_r = gam_r q (t nu)^r     (r=2..5),  gam_r = P_{r-1}(y) (1-y)^{-r} b^{-r/2},  y = e^{-nu}
   C_1 = beta1 t^3 + gam1 q (t nu),  beta1 = b^{-1/2},  gam1 = b^{-1/2}/(1-q)
Every error polynomial of Lemma A (nonnegative coefficients) becomes a sum of monomials
   c * t^p q^bq nu^j,   c independent of theta (upper-bounded on the cell).
Divided by theta q = t^2 q:  c t^{P} e^{-(bq-1) nu'} nu^j with P = p + 2bq - 4.
Lemma II-mono: t^P nu^j is nondecreasing in t on (0, sqrt(TH1)] iff P=j=0, or P>0 and nu_min > 2j/P.
The script ASSERTS this for every monomial, then bounds it by its value at t1 = sqrt(TH1), nu = nu'+log(1/TH1)."""
import sys
from flint import arb, ctx
import lemA_bound as LB
from ws_p12_common import F_fun, P as Eul
from ws_p12_eps import TV
TH1 = arb(sys.argv[1]); VLO = float(sys.argv[2]); VHI = float(sys.argv[3]); CS = arb(sys.argv[4])
h = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05
z2 = arb(2).zeta(); L1 = (1 / TH1).log(); t1 = TH1.sqrt()

class GP:
    """graded polynomial: {(p, bq, j): arb >= 0}"""
    def __init__(s, d=None): s.d = dict(d or {})
    def __add__(s, o):
        if not isinstance(o, GP): o = GP({(0, 0, 0): arb(o)})
        r = dict(s.d)
        for k, v in o.d.items(): r[k] = r.get(k, arb(0)) + v
        return GP(r)
    __radd__ = __add__
    def __mul__(s, o):
        if not isinstance(o, GP): return GP({k: v * o for k, v in s.d.items()})
        r = {}
        for k1, v1 in s.d.items():
            for k2, v2 in o.d.items():
                k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2]); r[k] = r.get(k, arb(0)) + v1 * v2
        return GP(r)
    __rmul__ = __mul__
    def __truediv__(s, o): return s * (1 / arb(o))
    def __pow__(s, e):
        r = GP({(0, 0, 0): arb(1)})
        for _ in range(e): r = r * s
        return r

def maxval(g, nmin, v0, v1, shift):
    """sup over theta<=TH1, nu' in [v0,v1] of g / (t^2 q)^shift   (shift in {0,1})."""
    tot = arb(0)
    for (p, bq, j), c in g.d.items():
        Pp = p + 2 * bq - 4 * shift; be = bq - shift
        assert (Pp == 0 and j == 0) or (Pp > 0 and nmin > arb(2 * j) / Pp), ("non-monotone monomial", p, bq, j)
        ex = (-arb(v0) * be).exp() if be >= 0 else (arb(v1) * (-be)).exp()
        tot += c * t1 ** Pp * (arb(v1) + L1) ** j * ex
    return tot

def cell(v0, v1):
    nmin = arb(v0) + L1                          # min of nu over the cell (attained at theta = TH1)
    assert nmin > 3
    smax = {r: arb.fac_ui(r) * z2 + TH1 * TV[r] for r in range(2, 7)}
    bmin = F_fun(2, nmin - TH1) - TH1 * TV[2]; assert bmin > 1
    y = (-nmin).exp(); q1 = TH1 * (-arb(v0)).exp()
    a = {r: smax[r] * bmin ** (-arb(r) / 2) for r in range(3, 7)}
    A = {r: GP({(r - 2, 0, 0): a[r]}) for r in range(3, 7)}
    C = {r: GP({(r, 1, r): Eul(r - 1, y) / (1 - y) ** r * bmin ** (-arb(r) / 2)}) for r in range(2, 6)}
    C[1] = GP({(3, 0, 0): bmin ** arb(-0.5), (1, 1, 1): bmin ** arb(-0.5) / (1 - q1)})
    delta0 = arb('0.16') * smax[4] / bmin / 12 + arb('0.0256') * smax[6] / bmin / 360   # theta-free
    assert delta0 < 1
    Bp = 1 - delta0; nf = (1 / Bp).sqrt(); one = arb(1)
    a3, a4, a5, a6 = A[3] / 6, A[4] / 24, A[5] / 120, A[6] / 720
    om = {3: a3, 4: a4, 5: a5}; c = C
    Z = GP()
    def nadd(*ps):
        out = {}
        for p in ps:
            for i, x in p.items(): out[i] = out.get(i, Z) + x
        return out
    def nmul(p1, p2):
        out = {}
        for i, x in p1.items():
            for j, yv in p2.items(): out[i + j] = out.get(i + j, Z) + x * yv
        return out
    def npow(p, e):
        out = {0: GP({(0, 0, 0): one})}
        for _ in range(e): out = nmul(out, p)
        return out
    def pint(co, B_):
        return sum((x * LB.abs_moment(p, B_) for p, x in co.items()), Z)
    def evenint(co):
        tot = Z
        for e, x in co.items():
            if e % 2 == 0:
                d = arb(1)
                for jj in range(1, e, 2): d *= jj
                tot = tot + x * d
        return tot
    Tw = nadd({6: a6}, {8: a4 * a4 / 2 + a3 * a5, 10: a5 * a5 / 2, 9: a4 * a5},
              nmul(nmul({4: a4, 5: a5}, npow(om, 2)), {0: GP({(0, 0, 0): arb(3) / 6})}),
              nmul(npow(om, 4), {0: GP({(0, 0, 0): arb(1) / 24})}))
    TH = {5: c[5] / 120, 3: c[1] * c[2] / 2 + c[1] ** 3 / 6}
    Qabs = {0: GP({(0, 0, 0): one}), 3: a3, 4: a4, 5: a5, 6: a3 * a3 / 2, 7: a3 * a4, 9: a3 ** 3 / 6}
    PH1abs = {1: c[1], 2: c[2] / 2 + c[1] ** 2 / 2, 3: c[3] / 6, 4: c[4] / 24}
    eNc = pint(nmul(Tw, {1: c[1]}), Bp) * nf + pint(nmul(Qabs, TH), one)
    eDc = pint(Tw, Bp) * nf
    polyNt = pint(nmul(Qabs, PH1abs), one / 2) * arb(2).sqrt()
    polyDt = pint(Qabs, one / 2) * arb(2).sqrt()
    k3, k4, k5 = A[3], A[4], A[5]
    Num = (c[4] * k3 ** 2 * (arb(25) / 48) + c[4] * k4 / 6 + (c[2] + c[1] ** 2) * (k4 ** 2 / 32 + k3 ** 4 * (arb(25) / 192)
           + k3 ** 2 * k4 * (arb(25) / 192)) + c[3] * k3 ** 3 * (arb(5) / 4) + c[3] * k5 * (arb(7) / 48) + c[3] * k3 * k4 * (arb(25) / 24))
    N0abs = evenint(nmul(Qabs, PH1abs))                  # |N0| <= sum_even |coef| (e-1)!!
    # theta-monotone scalar factors (Lemma II-tail / II-minor in md): value at TH1
    tailG = (-arb('0.04') * bmin / TH1).exp() / TH1 ** 2 * arb(v1).exp()          # e^{-w0^2/4}/(theta q)
    tail0 = (-arb('0.04') * bmin / TH1).exp()
    mnG = 2 * arb.pi() * (-CS / TH1).exp() * (smax[2] / (2 * arb.pi() * TH1 ** 3)).sqrt() / TH1 ** 2 * arb(v1).exp()
    mn0 = 2 * arb.pi() * (-CS / TH1).exp() * (smax[2] / (2 * arb.pi() * TH1 ** 3)).sqrt()
    assert TH1 < arb('0.02') * bmin and TH1 < CS / arb('3.5')
    # sup of plain quantities (shift 0) and of quantities/(theta q) (shift 1)
    M = lambda g, s: maxval(g, nmin, v0, v1, s)
    ED0 = M(eDc, 0) + tail0 * M(polyDt, 0) + mn0
    EDG = M(eDc, 1) + tailG * M(polyDt, 0) + mnG
    ENG = M(eNc, 1) + tailG * M(polyNt, 0) + 2 * mnG
    EN0 = M(eNc, 0) + tail0 * M(polyNt, 0) + 2 * mn0
    D0low = 1 - arb(5) / 24 * a[3] ** 2 * TH1          # D0 = 1 + A4/8 - 5A3^2/24 (scaled), A4 >= 0
    r0 = M(N0abs, 0) / D0low
    Dm = D0low - ED0; assert Dm > 0
    Nnlow = D0low
    DeltaG = (ENG + r0 * EDG) / Dm + M(Num, 1) / Nnlow
    Delta0 = (EN0 + r0 * ED0) / Dm + M(Num, 0) / Nnlow
    # M2 lower bound: M2 = 1 + C1 E1 - (C2 + C1^2) E2/2 - 5 C3 A3/12 + C4/8 >= 1 - |C1||E1| - (C2+C1^2)|E2|/2 - 5 C3 A3/12
    A3m, A4m, A5m = [M(A[r], 0) for r in (3, 4, 5)]
    C1m, C2m, C3m = M(C[1], 0), M(C[2], 0), M(C[3], 0)
    E1m = (A3m / 2 + A5m / 8 + 105 * A3m * A4m / 144 + 945 * A3m ** 3 / 1296) / D0low
    E2m = 1 + A4m / 2
    M2low = 1 - C1m * E1m - (C2m + C1m ** 2) * E2m / 2 - 5 * C3m * A3m / 12
    assert M2low > Delta0
    return DeltaG / (M2low - Delta0)

if __name__ == "__main__":
    v = VLO; worst = 0
    while v < VHI - 1e-12:
        w = min(v + h, VHI)
        g = float(cell(v, w).upper()); worst = max(worst, g)
        print(f"nu' in [{v:+.3f},{w:+.3f}]  Gamma <= {g:.5f}", flush=True)
        v = w
    print(f"Region II (theta <= {sys.argv[1]}): Gamma <= {worst:.5f} over nu' in [{VLO},{VHI}]")
