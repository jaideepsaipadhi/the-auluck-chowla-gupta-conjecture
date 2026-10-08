"""MECH-2: Lemma B (B3) window edges re-run with the edge values CERTIFIED FOR ALL n >= 10^5 by ws_mech_window.py
(adaptive arb sweep over beta in [beta_A, inf), no monotonicity-in-beta argument):
   nu'_L <= -1.906674  (printed majorant; the sharp expression is <= -1.94915),   nu'_R >= 1.043052.
The scan code is ws_b_edges.py's (Region I boxes via ws_b_core.L2_over_theta_box; Region II via |log M2|/theta <= 0.14236)."""
import sys
from flint import arb
from ws_b_core import L2_over_theta_box
TH1, THA, LM = arb('0.002'), arb('0.0040939'), arb('0.14236')
for VLs, tag in (('-1.906674', 'printed majorant'), ('-1.94915', 'sharp expression')):
    VL = arb(VLs); VR = arb('1.043052'); VLO = arb('-2.52'); VHI = arb('2.10')
    left2 = -1 + (-VL).exp() - LM
    right2 = -1 + (-VR).exp() / (1 - TH1 * (-VR).exp()) + LM
    assert left2 > 0 and right2 < 0
    ths = []; t = float(TH1.lower())
    while t < float(THA.upper()): ths.append(t); t *= 1.02
    ths.append(float(THA.upper()))
    def scan(v0, v1, sign):
        worst = None; h = 0.02; v = v0
        while v < v1:
            w = min(v + arb(h), v1)
            for a, c in zip(ths[:-1], ths[1:]):
                L, _ = L2_over_theta_box(arb(a).union(arb(c)), v.union(w))
                assert (sign > 0 and L > 0) or (sign < 0 and L < 0)
                val = float(L.lower()) if sign > 0 else float(L.upper())
                worst = val if worst is None else (min(worst, val) if sign > 0 else max(worst, val))
            v = w
        return worst
    print(f"[{tag}] nu'_L <= {VLs}, nu'_R >= 1.043052")
    print(f"  Region II: left L2/theta >= {left2.str(5)}, right L2/theta <= {right2.str(5)}")
    print(f"  Region I : left min L2/theta = {scan(VLO, VL, +1):.5f}, right max L2/theta = {scan(VR, VHI, -1):.5f}")
print("B3 edges: OK for all n >= 10^5 (edge constants certified in ws_mech_window.out)")
