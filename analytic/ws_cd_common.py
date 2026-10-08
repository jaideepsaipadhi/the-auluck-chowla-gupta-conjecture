"""Shared rigorous pieces for Lemmas C/D (python-flint arb).  Proofs: LemmaCD.md sections 2-3.

Scaled notation: v = k*theta, nu = K*theta = v + theta, s_r = theta^{r+1} kappa_r = theta*sum_{i<=k} f_r(i theta),
b = s_2 = theta^3 B, c1*theta = theta + f1(nu), c2*theta^2 = f2(nu), f1(u) = u/(e^u-1), f2(u) = (u/2)^2/sinh^2(u/2).

Lemma C0 (crude saddle bound), phi0 = c0*theta, Re psi <= -(1-delta) B phi^2/2 on |phi| <= phi0:
   |R_k/g0 - 1| <= rho = (T1/(1-d)^{5/2} + T2/(1-d)^{3/2} + 2M) / (1 - erfc(w0/sqrt2) - sqrt(2/pi) A3/(3(1-d)^2) - M)
   T1 = c1 kappa3/(2B^2) = theta * (c1 theta) s3/(2 b^2),  T2 = (c2 + c1^2)/(2B) = theta * (f2 + (c1 theta)^2)/(2b),
   A3 = kappa3 B^{-3/2} = s3 theta^{1/2} b^{-3/2},  w0 = phi0 B^{1/2} = c0 (b/theta)^{1/2},
   M  = sqrt(2 pi B) sup_{phi0<=|phi|<=pi}|chi| <= sqrt(2 pi b) theta^{-3/2} exp(-Psi/(2 theta)),  theta*S >= Psi.
"""
from flint import arb
from ws_p12_common import g_fun, h_fun, F_fun
from ws_p12_eps import TV            # TV_r of f_r on [0,inf) (certified by ws_mech_tv.py, rounded up)
PI = arb.pi()
SQ2PI = (2 / PI).sqrt()
Z2 = arb(2).zeta()
PB = arb('0.6')

def f1(u): return u / u.expm1()
def f2(u): h = u / 2; return (h / h.sinh()) ** 2
def h3(u): return (1 + 1 / (u / 2).sinh() ** 2).log()

def sup_f(r, U=60, N=6000):
    """sup_{u>=0} f_r(u): on [a,b], f_r = g^r h_r <= g(b)^r h_r(a); f_r decreasing on [U,inf) (P1.2)."""
    best = arb(0)
    for j in range(N):
        a = arb(U) * j / N; b = arb(U) * (j + 1) / N
        v = g_fun(b) ** r * h_fun(r, a)
        if v.upper() > best.upper(): best = arb(v.upper())
    return best
MS = {3: arb('2.0302009'), 4: arb('6.1253939')}     # rounded-up output of sup_f(3), sup_f(4) (checked in __main__)

def lower_int(fun, a, b, N):
    """right-endpoint Riemann sum of a DECREASING nonnegative fun on [a,b]: rigorous lower bound (0 if empty)."""
    a = arb(arb(a).upper()); b = arb(arb(b).lower())
    if not (b > a): return arb(0)
    h = (b - a) / N
    return h * sum((fun(a + h * j) for j in range(1, N + 1)), arb(0))

def psi_minor(thmax, v1, v2, c0, NA=50, NR=50, verbose=False):
    """Lemma C-minor.  Psi with theta*S(phi) >= Psi for every phi in [c0 theta, pi], valid for all
    theta <= thmax and all k with v1 <= k theta <= v2 (v2 = None: no upper bound)."""
    thmax = arb(thmax); v1 = arb(v1); c0 = arb(c0)
    best = arb(0)
    for cMf in ([30] + ([4 * float(PI / v1.lower()), 8 * float(PI / v1.lower())] if v1 < 0.8 else [])):
        cM = arb(cMf)
        if not cM > c0: continue
        # (A) c in [c0, cM]
        worstA = None
        r = (cM / c0) ** (arb(1) / NA)
        for j in range(NA):
            ca = c0 * r ** j; cb = c0 * r ** (j + 1)
            bA = arb(0)
            for T in [PI * jj / 32 for jj in range(1, 17)]:
                s2 = (T.sin() / T) ** 2
                U = v1.min(2 * T / cb)
                val = lower_int(lambda u: (1 + s2 * ca * ca * f2(u)).log(), thmax, U, NR)
                if val.lower() > bA.lower(): bA = val
            if worstA is None or bA.lower() < worstA.lower(): worstA = bA
        # tail phi in [cM theta, pi]: min(B1, B2)
        bB1 = arb(0)
        for jj in range(2, 16):
            al = arb(jj) / 20
            rho = 1 - (2 * al + PB / 2) / PI
            if not rho > 0: continue
            C0th = (1 + al / PI) * (4 * al / cM + thmax)     # >= theta*C0, C0 = (1+alpha/pi)(2alpha/d+1), d >= cM theta/2
            sa = al.sin() ** 2
            u1 = (thmax + C0th) / rho + thmax
            val = rho * lower_int(lambda u: (1 + sa / (u / 2).sinh() ** 2).log(), u1, v1.min(arb(40)), NR * 3)
            if val.lower() > bB1.lower(): bB1 = val
        sb = (PB / 4).sin() ** 2
        bB2 = lower_int(lambda w: (1 + sb / w.sinh() ** 2).log(), thmax, (v1 / 2).min(arb(20)), NR * 6)
        tail = min(bB1.lower(), bB2.lower())
        psi = arb(min(worstA.lower(), tail))
        if verbose: print("  cM=%.1f A %.5f B1 %.5f B2 %.5f" % (cMf, float(worstA.lower()), float(bB1.lower()), float(bB2.lower())))
        if psi > best: best = psi
    return best

def X_cell(th, v1, v2, c0, psi=None):
    """Lemma C1.  Upper bound X >= rho/theta valid for every saddle theta <= th with k theta in [v1, v2]
    (v2 = None: [v1, inf)).  All ingredient bounds are nondecreasing in theta (proved in LemmaCD.md C1)."""
    th = arb(th); v1 = arb(v1); c0 = arb(c0)
    if v2 is None:
        s3 = 6 * Z2 + th * TV[3]; s4 = 24 * Z2 + th * TV[4]; bup = 2 * Z2
    else:
        v2 = arb(v2)
        s3 = (F_fun(3, v2) + th * TV[3]).min(v2 * MS[3])
        s4 = (F_fun(4, v2) + th * TV[4]).min(v2 * MS[4])
        bup = F_fun(2, v2)
    blo = F_fun(2, v1) - th
    if v2 is not None: blo = blo.max(v1 * f2(v2))
    assert blo > 0
    c1th = th + f1(v1)
    g2 = f2(v1)
    dT = c0 * c0 * s4 / (12 * blo)
    delta = dT
    if v2 is not None and c0 * v2 < 2:
        dC = 1 - (1 + c0 * c0).log() / (c0 * c0) * (1 - (c0 * v2) ** 2 / 12)
        delta = dT.min(dC)
    if not delta < 1: return arb('inf'), {}
    T1 = c1th * s3 / (2 * blo ** 2)
    T2 = (g2 + c1th ** 2) / (2 * blo)
    A3 = s3 * th.sqrt() / blo ** arb(1.5)
    w0 = c0 * (blo / th).sqrt()
    if psi is None: psi = psi_minor(th, v1, v2, c0)
    if not psi > 5 * th: return arb('inf'), {}
    M = (2 * PI * bup).sqrt() * th ** arb(-1.5) * (-psi / (2 * th)).exp()
    num = T1 / (1 - delta) ** arb(2.5) + T2 / (1 - delta) ** arb(1.5) + 2 * M / th
    den = 1 - (w0 / arb(2).sqrt()).erfc() - SQ2PI * A3 / (3 * (1 - delta) ** 2) - M
    if not den > 0: return arb('inf'), {}
    return num / den, dict(T1=T1, T2=T2, delta=delta, A3=A3, w0=w0, M=M, psi=psi)

C0LIST = [arb(x) for x in ('0.4', '0.6', '0.8', '1.0', '1.25', '1.5', '2.0')]
def X_best(th, v1, v2):
    best = (arb('inf'), None, None)
    for c0 in C0LIST:
        X, d = X_cell(th, v1, v2, c0)
        if X.is_finite() and (not best[0].is_finite() or X.upper() < best[0].upper()): best = (X, c0, d)
    return best

if __name__ == "__main__":
    for r in (3, 4):
        s = sup_f(r); print("sup f_%d <=" % r, s.str(9)); assert s.upper() <= MS[r].lower()
    for v1, v2, th in [(0.05, 0.055, 0.0004), (0.3, 0.33, 0.0041), (1, 1.1, 0.0041), (3, 3.3, 0.0041), (6, None, 0.0058)]:
        X, c0, d = X_best(th, v1, v2)
        print(v1, v2, th, "X =", X.str(5), "c0 =", c0, {kk: float(vv.mid()) for kk, vv in d.items()} if d else None)
