"""Region I: rigorous Gamma(nu'-cell) := max over the box of eps/(theta q), theta in [TH1, THA].
Usage: python3 ws_p12_boxes.py TH1 THA NUP_LO NUP_HI NU_MIN CSTAR [rθ=1.02] [h=0.02]
Box: theta in [t, t*rθ], nu' in [v, v+h]; boxes with nu = nu'+log(1/theta) < NU_MIN everywhere are
unreachable (Lemma W) and skipped; in partially-reachable boxes nu' is clipped to nu >= NU_MIN.
Output: per nu'-cell the max over theta of the rigorous upper bound eps_box/(theta_lo * q_lo)."""
import sys, time
from flint import arb
from ws_p12_eps import core, scaled_from_th_nup
TH1, THA, VLO, VHI, NUMIN, CS = [arb(a) for a in sys.argv[1:7]]
rt = float(sys.argv[7]) if len(sys.argv) > 7 else 1.02
h = float(sys.argv[8]) if len(sys.argv) > 8 else 0.02
ths = []; t = float(TH1.lower())
while t < float(THA.upper()): ths.append(t); t *= rt
ths.append(float(THA.upper()))
nv = int((float(VHI.upper()) - float(VLO.lower())) / h) + 1
t0 = time.time(); out = {}
for j in range(nv):
    v0 = arb(VLO.lower()) + arb(h) * j; v1 = v0 + arb(h)
    worst = 0.0
    for a, b in zip(ths[:-1], ths[1:]):
        ta, tb = arb(a), arb(b)
        vlo_eff = v0.max(NUMIN + ta.log()) if True else v0      # nu >= NUMIN  <=>  nu' >= NUMIN - log(1/theta)
        if not (v1 > NUMIN + tb.log() or v1 > NUMIN + ta.log()): continue
        vlo_eff = v0.max(NUMIN + ta.log())                       # valid lower end for all theta in [a,b]
        if not (v1 > vlo_eff): continue
        th = ta.union(tb); nup = vlo_eff.union(v1)
        try:
            A, C, w0, mn, q = scaled_from_th_nup(th, nup, CS)
            e, d = core(A, C, w0, mn)
            g = float((e / (ta * ta * (-v1).exp())).upper())
        except AssertionError:
            g = float('inf')
        worst = max(worst, g)
    out[j] = worst
    print(f"nu' in [{float(v0.mid()):+.3f},{float(v1.mid()):+.3f}]  Gamma <= {worst:.5f}   ({time.time()-t0:.0f}s)", flush=True)
