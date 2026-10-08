"""Pointwise consistency: M2_formula(exact ingredients) vs ws_p12_eps.core's M2, and exact step/(theta q) vs bounds.
Usage: python3 ws_b_check.py n k1 k2 ..."""
import sys
from flint import arb
import lemA_bound as LB
from ws_b_core import M2_formula
from ws_p12_eps import core, scaled_from_nk
n = int(sys.argv[1])
def L2(k):
    I = LB.ingredients(n, k); th = I['th']; B = I['B']
    s = {r: I[{2:'B',3:'k3',4:'k4',5:'k5'}[r]] * th ** (r + 1) for r in (2, 3, 4, 5)}
    M2, _ = M2_formula(th, I['K'] * th, s[2], s[3], s[4], s[5])
    A, C, w0, _ = scaled_from_nk(n, k)
    _, d = core(A, C, w0, arb(0))
    return -th - (1 - I['q']).log() + M2.log(), th, I['q'], M2, d['M2']
for k in map(int, sys.argv[2:]):
    a, th, q, M, Mc = L2(k); b = L2(k + 1)[0]
    print(k, "M2 agree:", (M - Mc).str(3), " L2=", a.str(6), " step/(theta q) =", ((a - b) / (th * q)).str(5))
