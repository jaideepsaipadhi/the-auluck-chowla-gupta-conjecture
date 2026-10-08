"""(P2) Rigorous minor-arc constant.  Usage: python3 ws_p12_minor.py THETA_A NU_A
Proves: for all theta <= THETA_A, all k with nu=(k+1)theta >= NU_A, all phi in [0.4 theta, pi]:
   S(phi) := sum_{i<=k} log(1 + sin^2(i phi/2)/sinh^2(i theta/2)) >= 2 c_* / theta,
i.e.  log|chi(phi)| = -S/2 <= -c_*/theta.  Prints c_* (rigorous lower bound).
Three regimes (c = phi/theta):
 (a) 0.4 <= c <= CA :  S >= (1/theta) int_{THETA_A}^{min(2T/c, NU_A)} log(1 + s^2 c^2 (u/2)^2/sinh^2(u/2)) du,
     s = sin T / T, any T in (0, pi/2].
 (M) CA*theta <= phi <= PB :  S >= (1/theta) rho int_{(2THETA_A+4alpha/CA)/rho+THETA_A}^{NU_A-THETA_A} log(1+sin^2 alpha/sinh^2(u/2)) du,
     rho(alpha,phi) = (pi - 2alpha - PB/2)/pi  (block counting, see LemmaA_uniform.md).
 (b) PB <= phi <= pi :  S >= (1/theta) int_{THETA_A}^{(NU_A-THETA_A)/2} log(1 + sin^2(PB/4)/sinh^2 v) dv.
All integrands are decreasing in the integration variable, so right-endpoint Riemann sums
(evaluated in arb) are rigorous lower bounds.  Monotonicity in the parameter (c or phi) is used to
reduce each regime to finitely many evaluations (see LemmaA_uniform.md, Lemma P2)."""
import sys, math
from flint import arb
thA = arb(sys.argv[1]); nuA = arb(sys.argv[2])
CA = 30; PB = arb('0.6'); NR = 3000
pi = arb.pi()

def lower_int(fun, a, b, N=NR):
    """Right-endpoint Riemann sum of a DECREASING nonneg fun on [a,b] (a<b): rigorous lower bound."""
    a = arb(a.upper()); b = arb(b.lower())        # shrink to exact endpoints inside the true interval
    if not (b > a): return arb(0)
    h = (b - a) / N
    return h * sum((fun(a + h * j) for j in range(1, N + 1)), arb(0))

def ell(u, s, c):
    return (1 + (s * c * u / 2) ** 2 / (u / 2).sinh() ** 2).log()

# regime (a): c in [0.4, CA], geometric boxes
Ts = [arb(j) / 40 * pi / 2 for j in range(12, 41)]
best_a = None
cs = [0.4 * (CA / 0.4) ** (j / 300) for j in range(301)]
worst_a = arb(1e9)
for c1f, c2f in zip(cs[:-1], cs[1:]):
    c1 = arb(c1f).lower() if False else arb(c1f); c2 = arb(c2f)
    vals = []
    for T in Ts:
        s = T.sin() / T
        U = (2 * T / c2).min(nuA)
        vals.append(lower_int(lambda u: ell(u, s, c1), thA, U, 400))
    v = max(vals, key=lambda z: float(z.lower()))
    if float(v.lower()) < float(worst_a.lower()): worst_a = v; arg_a = c1f
print(f"(a) min over c in [0.4,{CA}] of theta*S >= {float(worst_a.lower()):.5f} (at c ~ {arg_a:.3f})", flush=True)

# regime (M)
best_M = arb(0)
for j in range(1, 60):
    al = arb(j) / 40
    rho = (pi - 2 * al - PB / 2) / pi           # d = phi/2 <= PB/2
    if not rho > 0: continue
    sa2 = al.sin() ** 2
    lo = (2 * thA + 4 * al / CA) / rho + thA      # >= theta((1 + c0)/rho + 1),  c0 = 2alpha/d + 1 <= 4alpha/(CA theta) + 1
    v = rho * lower_int(lambda u: (1 + sa2 / (u / 2).sinh() ** 2).log(), lo, nuA - thA)
    if float(v.lower()) > float(best_M.lower()): best_M = v; arg_M = float(al.mid())
print(f"(M) theta*S >= {float(best_M.lower()):.5f} (alpha={arg_M:.3f})", flush=True)

# regime (b)
a = (PB / 4).sin() ** 2
vb = lower_int(lambda v: (1 + a / v.sinh() ** 2).log(), thA, (nuA - thA) / 2, 20000)
print(f"(b) theta*S >= {float(vb.lower()):.5f}", flush=True)
m = min(float(worst_a.lower()), float(best_M.lower()), float(vb.lower()))
print(f"c_* = min/2 >= {m/2:.5f}   [valid for theta <= {sys.argv[1]}, nu >= {sys.argv[2]}]")
