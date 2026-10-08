"""Lemma B3: sign of L2 at the window edges, all n >= N_A.  Usage: python3 ws_b_edges.py N_A TH1 THA LOGM2_II
Left edge k_L = min W (K <= mu-2beta+2):  nu'_L <= -2+(r-1)(lb-2)+2r/beta+log r <= -2 + eps(lb-1)/(1-eps) + 2/(beta(1-eps)) =: VL
Right edge k_R = max W (K >= mu+2beta):   nu'_R >= 2 - delta(lb+2) + log(1-delta) =: VR
(both monotone in beta in the favourable direction, LemmaB.md sec. 5).  Then:
 L2/theta > 0 for nu' <= VL and < 0 for nu' >= VR, over all theta <= theta_A (Region I boxes + Region II bound)."""
import sys
from flint import arb
from ws_b_core import L2_over_theta_box
NA, TH1, THA, LM = [arb(a) for a in sys.argv[1:5]]
z2 = arb(2).zeta()
b = (6 * NA).sqrt() / arb.pi(); lb = b.log()
eps = (lb + 2 + 1 / b) / (z2 * b); r = (1 / (1 - eps)).sqrt(); delta = 5 * lb / b
VL = -2 + eps * (lb - 1) / (1 - eps) + 2 / (b * (1 - eps))     # >= -2+(r-1)(lb-2)+2r/beta+log r; decreasing in beta
VR = 2 - delta * (lb + 2) + (1 - delta).log()                   # increasing in beta
print("nu'_L <=", float(VL.upper()), "  nu'_R >=", float(VR.lower()))
VL = arb(VL.upper()); VR = arb(VR.lower()); VLO = arb('-2.52'); VHI = arb('2.10')
# Region II (theta <= TH1): L2/theta = -1 + (-log(1-q))/theta + log M2/theta, |log M2|/theta <= LM;
#  q/theta <= -log(1-q)/theta <= q/(theta(1-q)), q/theta = e^{-nu'}
q1 = TH1 * (-VLO).exp()
left2 = -1 + (-VL).exp() - LM
right2 = -1 + (-VR).exp() / (1 - TH1 * (-VR).exp()) + LM
print("Region II: left L2/theta >=", left2.str(5), "  right L2/theta <=", right2.str(5)); assert left2 > 0 and right2 < 0
# Region I boxes
ths = []; t = float(TH1.lower())
while t < float(THA.upper()): ths.append(t); t *= 1.02
ths.append(float(THA.upper()))
def scan(v0, v1, sign):
    worst = None; h = 0.02; v = v0
    while v < v1:
        w = min(v + arb(h), v1)
        for a, c in zip(ths[:-1], ths[1:]):
            L, _ = L2_over_theta_box(arb(a).union(arb(c)), v.union(w))
            val = float(L.lower()) if sign > 0 else float(L.upper())
            assert (sign > 0 and L > 0) or (sign < 0 and L < 0), (a, float(v.mid()), L.str(4))
            worst = val if worst is None else (min(worst, val) if sign > 0 else max(worst, val))
        v = w
    return worst
print("Region I: min L2/theta on nu' in [-2.52, VL] =", scan(VLO, VL, +1))
print("Region I: max L2/theta on nu' in [VR, 2.10] =", scan(VR, VHI, -1))
print("B3 edges: OK")
