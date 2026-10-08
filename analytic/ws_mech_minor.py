"""MECH-1: fully mechanical minor-arc bound (replaces block counting / pair arguments / hand Abel step).
Usage: python3 ws_mech_minor.py THMAX UTOP C0      (P2 of LemmaA_uniform: 0.0042 3.0 0.4)
CLAIM.  For every theta in (0, THMAX], every k with (k+1)theta >= UTOP, and every phi in [C0 theta, pi]:
        theta * S(phi) >= PSI,   S(phi) = sum_{i<=k} log(1 + sin^2(i phi/2)/sinh^2(i theta/2))  ( = -2 log|chi(phi)| ).
PROOF (all ingredients verified in ws_mech_elem.py; every number below is an arb ball):
 (a) c = phi/theta in [C0, c_a] (geometric cells [c_lo, c_hi]).  For T in (0, pi/2] and i <= 2T/phi, i phi/2 <= T, so
     sin(i phi/2) >= (sin T/T)(i phi/2) (E3).  Summand >= l(i theta),  l(u) = log(1 + s^2 c_lo^2 f2(u)),  f2(u) = ((u/2)/sinh(u/2))^2,
     l decreasing (E2).  With I = min(k, floor(2T/phi)):  S >= sum_{i<=I} l(i theta) >= theta^{-1} int_theta^{(I+1)theta} l
     >= theta^{-1} int_{THMAX}^{min(UTOP, 2T/c_hi)} l     ((I+1)theta >= min((k+1)theta, 2T/c), l >= 0).
     The integral is bounded below by a right Riemann sum (l decreasing), evaluated in arb.  No dependence on theta remains.
 (D) c >= c_a (i.e. phi in [c_a theta, pi], d := phi/2 in [c_a theta/2, pi/2]).
     (D1) concavity (E4): log(1 + sin^2(id) A_i) >= sin^2(id) w_i,  w_i := log(1 + 1/sinh^2(i theta/2)) = 2 log coth(i theta/2), decreasing (E8).
     (D2) Abel (E16): sum_{i<=k} sin^2(id) w_i = sum_{I<=k} G_I (w_I - w_{I+1}),  w_{k+1} := 0,  G_I = sum_{i<=I} sin^2(id).
     (D3) Dirichlet (E9): G_I = I/2 - sin(Id)cos((I+1)d)/(2 sin d) >= (I - N0)/2 with N0 := ceil(1/sin d); also G_I >= 0.
          So S >= (1/2) sum_I (w_I - w_{I+1})(I - N0)_+ = (1/2) sum_{i=N0+1}^{k} w_i >= (1/2) theta^{-1} int_{(N0+1)theta}^{(k+1)theta} w.
     (D4) (N0+1) theta <= theta/sin d + 2 theta and sin d >= sin(c_a theta/2) >= (c_a theta/2) sinc(c_a THMAX/2) (sin increasing on [0,pi/2],
          sinc decreasing E3; asserted c_a THMAX/2 <= pi/2), so (N0+1)theta <= L := 2 THMAX + 2/(c_a sinc(c_a THMAX/2)).
     Hence theta S >= int_L^{UTOP} log coth(u/2) du  (right Riemann sum, integrand decreasing).  Valid for ALL phi in [c_a theta, pi] at once.
 PSI = max over c_a in a list of min( min_{cells in [C0, c_a]} (a), (D)(c_a) ).
No block counting, no pair argument, no density claim, no theta box cover, no theta -> 0 asymptotics: both (a) and (D) are theta-uniform
because theta enters only through the bounds theta <= THMAX and (k+1)theta >= UTOP."""
import sys, time
from flint import arb
PI = arb.pi()
def f2(u): h = u / 2; return (h / h.sinh()) ** 2
def lower_int(fun, a, b, N):
    a = arb(arb(a).upper()); b = arb(arb(b).lower())
    if not (b > a): return arb(0)
    h = (b - a) / N
    return h * sum((fun(a + h * j) for j in range(1, N + 1)), arb(0))
def logcoth(u): return (1 / (u / 2).tanh()).log()
def regime_a_cells(thmax, utop, c0, cmax, ncell=120, NR=600):
    thmax, utop, c0, cmax = map(arb, (thmax, utop, c0, cmax))
    r = (cmax / c0) ** (arb(1) / ncell); out = []
    Ts = [PI / 2 * jj / 48 for jj in range(8, 49)]
    for j in range(ncell):
        clo = c0 * r ** j; chi = c0 * r ** (j + 1)
        best = arb(0)
        for T in Ts:
            s2 = (T.sin() / T) ** 2
            U = utop.min(2 * T / chi)
            v = lower_int(lambda u: (1 + s2 * clo * clo * f2(u)).log(), thmax, U, NR)
            if v.lower() > best.lower(): best = arb(v.lower())
        out.append((clo, chi, best))
    return out
def regime_D(thmax, utop, ca, NR=4000):
    thmax, utop, ca = arb(thmax), arb(utop), arb(ca)
    x = ca * thmax / 2
    assert x <= PI / 2
    L = 2 * thmax + 2 / (ca * x.sin() / x)
    return lower_int(logcoth, L, utop, NR), L
def psi_mech(thmax, utop, c0, cmax=60, ncell=400, verbose=False):
    cells = regime_a_cells(thmax, utop, c0, cmax, ncell)
    best = (arb(0), None)
    run = None
    for idx, (clo, chi, va) in enumerate(cells):
        run = va if run is None else run.min(va)        # min of (a) over [C0, chi]
        if idx % 4 == 3 or idx == len(cells) - 1:
            if chi * arb(thmax) / 2 > PI / 2: break
            vd, L = regime_D(thmax, utop, chi, 3000)
            cand = run.min(vd)
            if cand.lower() > best[0].lower(): best = (arb(cand.lower()), (float(chi.mid()), float(run.lower()), float(vd.lower()), float(L.mid())))
    return best
if __name__ == "__main__":
    t0 = time.time()
    TH, UT, C0 = sys.argv[1], sys.argv[2], sys.argv[3]
    psi, info = psi_mech(TH, UT, C0)
    ca, va, vd, L = info
    print(f"theta <= {TH}, (k+1)theta >= {UT}, phi in [{C0} theta, pi]")
    print(f"  regime (a), c in [{C0}, {ca:.3f}]:  theta*S >= {va:.5f}")
    print(f"  regime (D), phi in [{ca:.3f} theta, pi] (L = {L:.4f}):  theta*S >= {vd:.5f}")
    print(f"  PSI = {float(psi.lower()):.5f}   =>  c_* = PSI/2 >= {float(psi.lower())/2:.5f}")
    print(f"time {time.time()-t0:.1f}s")
