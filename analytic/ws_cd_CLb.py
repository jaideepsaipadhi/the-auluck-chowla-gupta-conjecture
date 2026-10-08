"""(C-left, v >= V0) Lemma CL-b.  Usage: python3 ws_cd_CLb.py N_A V0
For n >= N_A and k <= mu - 2beta with v = k theta >= V0:   R_k > 1.
 log R_k >= -theta - log(1-q) + log(1 - theta X) >= -theta + q - theta X/(1 - theta X),
 so R_k > 1 if  q/theta > 1 + X/(1 - th X).  Two lower bounds for q/theta:
   (i)  q/theta = e^{-nu'} >= e^{-nuCL}  (ws_cd_window.py),
   (ii) q/theta = e^{-nu}/theta >= e^{-(v2+th)}/th  for theta <= th, v <= v2.
 v-cells [v1, v2] geometric from V0, last cell [v1, inf).  theta <= th := min(thA, v2/K0) (k >= K0);
 X = X_best(th, v1, v2) valid for all theta <= th (monotone in theta)."""
import sys
from flint import arb
from ws_cd_common import X_best
NA = arb(sys.argv[1]) if len(sys.argv) > 1 else arb(100000)
V0 = arb(sys.argv[2]) if len(sys.argv) > 2 else arb('0.1')
K0 = int(sys.argv[3]) if len(sys.argv) > 3 else 78     # every k in CL satisfies k > k_D(n) >= K0
b = (6 * NA).sqrt() / arb.pi(); lb = b.log(); z2 = arb(2).zeta()
eps = (lb + 2 + 1 / b) / (z2 * b); r = (1 / (1 - eps)).sqrt()
thA = arb((r / b).upper())
nuCL = (lb - 2 + 1 / b) * r - lb + r.log()
QL = (-nuCL).exp()
print("thA =", thA.str(6), " nuCL =", nuCL.str(6), " e^{-nuCL} =", QL.str(6))
v1 = V0; worst = arb(1e9); ratio = arb('1.15')
while True:
    last = v1 > 8
    v2 = None if last else v1 * ratio
    th = thA if last else thA.min(v2 / K0)              # theta = v/k <= v2/K0
    th = arb(th.upper())
    X, c0, d = X_best(th, v1, v2)
    assert th * X < 1
    need = 1 + X / (1 - th * X)
    qt = QL if last else QL.max((-(v2 + th)).exp() / th)
    marg = qt / need
    ok = marg > 1
    print(f"v in [{float(v1.mid()):.4f}, {'inf' if last else '%.4f' % float(v2.mid())}]  X <= {float(X.upper()):9.4f}  c0={float(c0.mid()):.2f}  q/theta >= {float(qt.lower()):9.3f}  margin {float(marg.lower()):.3f}", flush=True)
    assert ok
    worst = worst.min(marg)
    if last: break
    v1 = v2
print("CL-b PROVED for n >=", sys.argv[1] if len(sys.argv) > 1 else 100000, "; min margin (q/theta)/(1+X/(1-th X)) >=", float(worst.lower()))
