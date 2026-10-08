"""Shared rigorous primitives for the uniform Lemma A (P1/P2).  python-flint arb throughout.

f_r(u) := u^r Li_{1-r}(e^{-u}) = g(u)^r * h_r(u),  g(u) = u/(1-e^{-u}) (increasing, g(0)=1),
          h_r(u) = y P_{r-1}(y), y = e^{-u}, P_s the Eulerian polynomial (Li_{-s}(y) = y P_s(y)/(1-y)^{s+1}),
          h_r decreasing in u (positive coefficients, y decreasing).
F_r(v) := int_0^v f_r = r! ( zeta(2) - sum_{j=0}^r v^j/j! Li_{2-j}(e^{-v}) ).
theta^{r+1} kappa_r = theta * sum_{i<=k} f_r(i theta);  c_r theta^r = nu^r Li_{1-r}(q)  (exact).
"""
from flint import arb, fmpz
EUL = {0: [1], 1: [1], 2: [1, 1], 3: [1, 4, 1], 4: [1, 11, 11, 1], 5: [1, 26, 66, 26, 1], 6: [1, 57, 302, 302, 57, 1]}

def P(s, y):
    out = arb(0)
    for c in reversed(EUL[s]): out = out * y + c
    return out

def Li_neg(s, y):            # Li_{-s}(y), s >= 0
    return y * P(s, y) / (1 - y) ** (s + 1)

def Li(sidx, y):             # Li_{sidx}(y) for sidx <= 2
    if sidx == 2: return y.polylog(2) if hasattr(y, 'polylog') else arb(y).polylog(2)
    if sidx == 1: return -(1 - y).log()
    return Li_neg(-sidx, y)

def g_fun(u):
    u = arb(u)
    if u == 0: return arb(1)
    return u / (-(-u).expm1())

def h_fun(r, u):
    y = (-arb(u)).exp(); return y * P(r - 1, y)

def f_fun(r, u):
    return g_fun(u) ** r * h_fun(r, u)

def F_fun(r, v):
    """int_0^v f_r, v an arb ball with v > 0."""
    y = (-v).exp()
    z2 = arb(2).zeta()
    s = arb(0); vj = arb(1); fj = arb(1)
    for j in range(0, r + 1):
        if j > 0: vj *= v; fj *= j
        s += vj / fj * Li(2 - j, y)
    return arb.fac_ui(r) * (z2 - s)

def fact(r): return arb.fac_ui(r)
