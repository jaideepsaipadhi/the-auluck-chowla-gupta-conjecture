"""Lemma B, Region I (theta in [TH1, THA]): rigorous lower bound of
   Sigma := (L2_k - L2_{k+1}) / (theta_k q_k)  >=  T2/(theta q) - |Delta log M2|/(theta q)
over boxes (theta, nu') covering the window (Lemma W), clipped to nu >= NU_MIN.
Usage: python3 ws_b_step.py TH1 THA NUP_LO NUP_HI NU_MIN [rt] [h]"""
import sys, time
from flint import arb
from ws_b_core import dlogM2_box, T2_over_thq, TAU
TH1, THA, VLO, VHI, NUMIN = [arb(a) for a in sys.argv[1:6]]
rt = float(sys.argv[6]) if len(sys.argv) > 6 else 1.02
h = float(sys.argv[7]) if len(sys.argv) > 7 else 0.02
ths = []; t = float(TH1.lower())
while t < float(THA.upper()): ths.append(t); t *= rt
ths.append(float(THA.upper()))
nv = int(round((float(VHI.upper()) - float(VLO.lower())) / h + 0.4999)) + 0
t0 = time.time(); worst = 1e9; nbox = 0
for j in range(nv):
    v0 = arb(VLO.lower()) + arb(h) * j; v1 = v0 + arb(h)
    cw = 1e9
    for a, b in zip(ths[:-1], ths[1:]):
        ta, tb = arb(a), arb(b)
        vlo_eff = v0.max(NUMIN + ta.log())
        if not (v1 > vlo_eff): continue
        th = ta.union(tb); nup = vlo_eff.union(v1)
        tot, dl, _ = dlogM2_box(th, nup)
        T2 = T2_over_thq(th, nup)
        sig = T2 - tot / (th * dl['q'])
        cw = min(cw, float(sig.lower())); nbox += 1
    worst = min(worst, cw)
    print(f"nu' in [{float(v0.mid()):+.3f},{float(v1.mid()):+.3f}]  Sigma >= {cw:.5f}  ({time.time()-t0:.0f}s)", flush=True)
need = 2 * 0.4282 * (1 + float(TAU.mid()))
print(f"Region I: min Sigma >= {worst:.5f} over {nbox} boxes;  need 2*Gamma_max*(1+tau) = {need:.4f}:", "OK" if worst > need else "FAIL")
