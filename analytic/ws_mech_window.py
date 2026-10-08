"""MECH-2/3b: Lemma W, Lemma CW and the Lemma-B edge bounds, recomputed from scratch and certified for ALL n >= N_A
WITHOUT any monotonicity-in-beta argument.  Imports nothing from the project.  Usage: python3 ws_mech_window.py N_A C
Derivation re-done (MECH.md item 3b).  With x = 1/beta, l = log beta (so x = e^{-l}), every quantity below is written as a
function of the four "atoms"  x, lx = l x, l2x = l^2 x, il = 1/l   (plus l itself only in nu_min, which is bounded below by l).
  delta = C lx,   eps = (lx + 2x + x^2)/zeta2,   r = (1-eps)^{-1/2},  eps*l = (l2x + 2 lx + x lx)/zeta2
  G      = 2C z2 - C^2 z2 lx - e^2 exp(C l2x)/(1 - e^{-v}) - il ,   e^{-v} = e^2 x exp(C (l2x - 2 lx))   [v = (l-2)(1-delta)]
  thA    = r x                                  (theta <= r/beta)
  nuhi   = (r-1) l + (2+x) r + log r            (nu' upper bound on W)
  nulo   = -2 - C(l2x - 2 lx) + log(1 - C lx)   (nu' lower bound on W)
  numin  = l - 2 - C l2x + 2C lx                (nu lower bound on W)
  VL     = -2 + (r-1)(l-2) + 2 r x + log r      (nu' at the left end of W; LemmaB sec 5, sharper than the printed majorant)
  VLprinted = -2 + eps(l-1)/(1-eps) + 2x/(1-eps)
  VR = nuCR = 2 - C l2x - 2C lx + log(1 - C lx) (nu' at the right end of W / in region CR)
  nuCL   = (r-1) l - 2 r + x r + log r          (nu' on region CL, K <= mu - 2beta + 1)
  (r-1) l is computed as ((r-1)/eps) * (eps l).
Sweep: l in [l_A, 60] by boxes of width h (each atom enclosed by direct ball evaluation of its formula on the l-box),
tail l in [60, inf): x in [0, e^{-60}], lx in [0, 60 e^{-60}], l2x in [0, 3600 e^{-60}] (l^j e^{-l} decreasing for l > j, E6), il in [0, 1/60],
and l >= 60 where l enters alone (numin).
Output: sup / inf of each quantity over ALL beta >= beta_A, to compare with the project constants."""
import sys, time
from flint import arb
NA = arb(sys.argv[1]) if len(sys.argv) > 1 else arb(100000)
C = arb(sys.argv[2]) if len(sys.argv) > 2 else arb(5)
z2 = arb(2).zeta(); e = arb(1).exp()
bA = (6 * NA).sqrt() / arb.pi(); lA = bA.log()
def quantities(x, lx, l2x, il, l_low):
    delta = C * lx
    eps = (lx + 2 * x + x * x) / z2
    epsl = (l2x + 2 * lx + x * lx) / z2
    assert eps < 1 and delta < 1
    r = 1 / (1 - eps).sqrt()
    rm1 = (r - 1) / eps if eps.lower() > 0 else (r - 1)  # see below
    # (r-1)/eps = 1/(sqrt(1-eps)(1+sqrt(1-eps))) : stable form, no division by a ball containing 0
    s = (1 - eps).sqrt(); rm1_over_eps = 1 / (s * (1 + s))
    rm1l = rm1_over_eps * epsl
    emv = e * e * x * (C * (l2x - 2 * lx)).exp()
    assert emv < 1
    G = 2 * C * z2 - C * C * z2 * lx - e * e * (C * l2x).exp() / (1 - emv) - il
    out = dict(
        G=G, thA=r * x,
        nuhi=rm1l + (2 + x) * r + r.log(),
        nulo=-2 - C * (l2x - 2 * lx) + (1 - delta).log(),
        numin=l_low - 2 - C * l2x + 2 * C * lx,
        VL=-2 + rm1l - 2 * rm1_over_eps * eps + 2 * r * x + r.log(),
        VLprinted=-2 + (epsl - eps) / (1 - eps) + 2 * x / (1 - eps),
        nuCR=2 - C * l2x - 2 * C * lx + (1 - delta).log(),
        nuCL=rm1l - 2 * r + x * r + r.log())
    return out
SUP = ['thA', 'nuhi', 'VL', 'VLprinted', 'nuCL']; INF = ['G', 'nulo', 'numin', 'nuCR']
PROJ = dict(G=('>=', '0'), thA=('<=', '0.0040939'), nuhi=('<=', '2.0840'), VL=('<=', '-1.9067'), VLprinted=('<=', '-1.906674'),
            nuCL=('<=', '-1.95361'), nulo=('>=', '-2.5102'), numin=('>=', '3.1158'), nuCR=('>=', '1.04305'))
def box(a, b):
    l = a.union(b); x = (-l).exp()
    return quantities(x, l * x, l * l * x, 1 / l, a)
def holds(q, k):
    op, val = PROJ[k]; return q[k] <= arb(val) if op == '<=' else q[k] >= arb(val)
def certify(a, b, keys, depth, stats):
    """prove every claim in keys on l in [a,b] (adaptive bisection, depth <= 60)"""
    q = box(a, b); bad = [k for k in keys if not holds(q, k)]
    stats['boxes'] += 1
    for k in keys:
        if k not in bad:
            ext = q[k].upper() if PROJ[k][0] == '<=' else q[k].lower()
            stats['ext'][k] = ext if k not in stats['ext'] else (max(stats['ext'][k], ext) if PROJ[k][0] == '<=' else min(stats['ext'][k], ext))
    if not bad: return
    assert depth < 60, ("cannot certify", bad, a.str(12), b.str(12))
    m = (a + b) / 2
    certify(a, m, bad, depth + 1, stats); certify(m, b, bad, depth + 1, stats)
if __name__ == "__main__":
    t0 = time.time()
    h = arb(sys.argv[3]) if len(sys.argv) > 3 else arb('0.01')
    stats = dict(boxes=0, ext={})
    a = arb(lA.lower()); LT = arb(60); keys = list(PROJ)
    while a < LT:
        b = (a + h).min(LT); certify(a, b, keys, 0, stats); a = b
    xt = arb(0).union((-LT).exp()); lxt = arb(0).union(LT * (-LT).exp()); l2xt = arb(0).union(LT * LT * (-LT).exp())
    qt = quantities(xt, lxt, l2xt, arb(0).union(1 / LT), LT)
    for k in keys: assert holds(qt, k), ("tail fails", k)
    x = 1 / bA; qA = quantities(x, lA * x, lA * lA * x, 1 / lA, lA)
    print(f"beta_A = {bA.str(10)}, l_A = {lA.str(8)}, C = {C.str(3)}; sweep l in [l_A, 60] ({stats['boxes']} adaptive boxes) + tail [60, inf)")
    print("claim                          value at beta_A          worst certified box bound")
    for k, (op, val) in PROJ.items():
        print(f"  {k:10s} {op} {val:11s}   {qA[k].str(10):26s} {float(stats['ext'][k]):.10f}   CERTIFIED for all n >= N_A")
    print(f"time {time.time() - t0:.1f}s")
    # independent sanity (not part of the proof): actual saddles at the window ends for a few n
    def S(t, k): return sum((arb(i) / (arb(i) * t).expm1() for i in range(1, k + 1)), arb(0))
    for n in (100000, 1000000, 10000000):
        be = (6 * arb(n)).sqrt() / arb.pi(); l = be.log(); mu = be * l
        for k in (int((mu - 2 * be).floor().unique_fmpz()) + 1, int((mu + 2 * be).ceil().unique_fmpz()) - 1):
            m = n - k; lo, hi = arb('1e-6'), arb('0.1')
            for _ in range(60):
                mid = (lo + hi) / 2
                if S(mid, k) > m: lo = mid
                else: hi = mid
            th = lo; nup = (k + 1) * th + th.log()
            dl = C * l / be; eps = (l + 2 + 1 / be) / (z2 * be)
            print(f"  sanity n={n} k={k}: theta*beta={float((th*be).mid()):.6f} in [(1-delta)={float((1-dl).mid()):.5f}, r={float((1/(1-eps)).sqrt().mid()):.5f}],  nu'={float(nup.mid()):+.5f}")
    # Lemma W reduction chain, re-derived (MECH.md 3b); arb spot checks of each inequality at random admissible (beta, k, theta) -- consistency only
    import random; random.seed(3)
    def T_exact(v): y = (-v).exp(); return -v * (-y).log1p() + y.polylog(2)
    for _ in range(300):
        v = arb(random.uniform(0.5, 15)); y = (-v).exp()
        assert T_exact(v) <= (v + 1) * y / (1 - y)
    print("  T(v) <= (v+1)e^{-v}/(1-e^{-v}) spot-checked (proof: Li1(y), Li2(y) <= y/(1-y) termwise)")
