"""MECH-1 (LemmaCD part): re-run Lemmas CR and CL-b with the C-minor constant Psi replaced by the MECHANICAL bound of
ws_mech_minor.py (no block counting (B1), no pair argument (B2)).  Usage: python3 ws_mech_cd.py N_A
Setting of Lemma C-minor: theta <= thmax, k theta >= v1 (so (k+1)theta >= v1), phi in [c0 theta, pi]  ->  ws_mech_minor with UTOP = v1.
Everything else (Lemma C0/C1 algebra, X_cell, window constants) is the project code; window constants are re-certified in ws_mech_window.py,
TV constants in ws_mech_tv.py (project values are upper bounds of the certified ones, so remain valid)."""
import sys, time
from flint import arb
import ws_cd_common as CC
import ws_mech_minor as MM
PI = arb.pi()
_cache = {}
def psi_mech_cd(thmax, v1, v2, c0, **kw):
    thmax = arb(thmax); v1 = arb(v1); c0 = arb(c0)
    key = (float(thmax.upper()), float(v1.lower()))
    if key not in _cache:
        cmax = max(60.0, 30.0 / float(v1.lower()))
        cells = MM.regime_a_cells(thmax, v1, arb('0.4'), cmax, ncell=160, NR=200)
        Ds = {}
        for idx, (clo, chi, va) in enumerate(cells):
            if idx % 4 == 3 and chi * thmax / 2 <= PI / 2:
                Ds[idx] = MM.regime_D(thmax, v1, chi, 1500)[0]
        _cache[key] = (cells, Ds)
    cells, Ds = _cache[key]
    best = arb(0); run = None
    for idx, (clo, chi, va) in enumerate(cells):
        if chi <= c0: continue
        run = va if run is None else run.min(va)          # (a) over cells meeting [c0, chi]  (c_lo <= c0 is conservative: l increases in c)
        if idx in Ds:
            cand = run.min(Ds[idx])
            if cand.lower() > best.lower(): best = arb(cand.lower())
    return best
CC.psi_minor = psi_mech_cd
if __name__ == "__main__":
    t0 = time.time()
    NA = sys.argv[1] if len(sys.argv) > 1 else '100000'
    print("=== Lemma CR with mechanical Psi ==="); sys.argv = ['ws_cd_CR.py', NA]
    exec(open('ws_cd_CR.py').read().split('"""', 2)[2], {'__name__': 'x'})
    print("=== Lemma CL-b with mechanical Psi ==="); sys.argv = ['ws_cd_CLb.py', NA, '0.1', '78']
    exec(open('ws_cd_CLb.py').read().split('"""', 2)[2], {'__name__': 'x'})
    print(f"time {time.time()-t0:.1f}s")
