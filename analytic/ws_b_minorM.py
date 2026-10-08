"""Audit of Lemma P2 regime (M): 30 theta <= phi <= 0.6.  Usage: python3 ws_b_minorM.py THETA_A NU_A
Independent re-derivation (LemmaB.md sec. 7).  d = phi/2 in [15 theta, 0.3], alpha in (0, pi/2), d < pi - 2 alpha.
 (i)   good := {i : dist(i d, pi Z) >= alpha}.  Maximal bad runs have <= N_b := floor(2alpha/d)+1 elements; two runs are
       separated by >= N_g := ceil((pi-2alpha)/d) - 1 good indices.  Hence G(I) := #good in [1,I] >= rho (I - N_b),
       rho = N_g/(N_g+N_b) >= (pi - 2alpha - d)/pi >= (pi - 2alpha - 0.3)/pi.
 (ii)  Abel: w_i := log(1 + sin^2 alpha / sinh^2(i theta/2)) decreasing, w_{k+1} := 0:
       S >= sum_good w_i = sum_I (w_I - w_{I+1}) G(I) >= rho sum_{i > N_b} w_i >= (rho/theta) int_{(N_b+1)theta}^{k theta} w(u/theta) du.
 (iii) (N_b+1) theta <= 2 theta + 2 alpha theta/d <= 2 theta_A + 4 alpha/30;  k theta = nu - theta >= NU_A - THETA_A.
Prints both the audited bound (lower limit u1 = 2thA + 4alpha/30) and the original script's (more conservative) one."""
import sys
from flint import arb
thA = arb(sys.argv[1]); nuA = arb(sys.argv[2]); pi = arb.pi()
def lower_int(fun, a, b, N=4000):
    a = arb(a.upper()); b = arb(b.lower())
    if not (b > a): return arb(0)
    h = (b - a) / N
    return h * sum((fun(a + h * j) for j in range(1, N + 1)), arb(0))
best = {}
for j in range(1, 60):
    al = arb(j) / 40
    rho = (pi - 2 * al - arb('0.3')) / pi   # density, LemmaB.md sec. 7 (i)
    if not rho > 0: continue
    sa2 = al.sin() ** 2
    w = lambda u: (1 + sa2 / (u / 2).sinh() ** 2).log()
    for tag, lo in (('audited', 2 * thA + 4 * al / 30), ('original', (2 * thA + 4 * al / 30) / rho + thA), ('CD-count', (thA + (1 + al / pi) * (4 * al / 30 + thA)) / rho + thA)):
        v = rho * lower_int(w, lo, nuA - thA)
        if tag not in best or float(v.lower()) > float(best[tag][0].lower()): best[tag] = (v, float(al.mid()))
for tag, (v, a) in best.items():
    print(f"(M) {tag}: theta*S >= {float(v.lower()):.5f}  (alpha = {a:.3f})")
