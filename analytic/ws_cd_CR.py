"""(C-right) Usage: python3 ws_cd_CR.py N_A
For n >= N_A and mu+2beta <= k <= ceil(n/2):  log R_k <= -theta - log(1-q) + log(1+rho) <= -theta + q/(1-q) + theta X,
so R_k < 1 if  e^{-nu'}/(1-q) + X < 1   (q/theta = e^{-nu'}).
Inputs (ws_cd_window.py): theta <= thCR, nu' >= nuCR, hence nu = nu' - log theta >= nuCR + log(1/thCR), v = k theta = nu - theta.
X from X_cell(thCR, vlo, None)."""
import sys
from flint import arb
from ws_cd_common import X_best, Z2
NA = arb(sys.argv[1]) if len(sys.argv) > 1 else arb(100000)
C = arb(5)
b = (6 * NA).sqrt() / arb.pi(); lb = b.log(); delta = C * lb / b
nuCR = (lb + 2) * (1 - delta) - lb + (1 - delta).log()
thCR = arb((2 * Z2 / (NA - 1)).sqrt().upper())
vlo = arb((nuCR - thCR.log() - thCR).lower())
X, c0, d = X_best(thCR, vlo, None)
q = (-vlo).exp()                       # q = e^{-nu} <= e^{-v}
lhs = (-nuCR).exp() / (1 - q) + X
print("thCR =", thCR.str(6), " nuCR =", nuCR.str(6), " v >=", vlo.str(6))
print("X <=", X.str(6), " (c0 =", c0, ")  psi =", d['psi'].str(5), " M =", d['M'].str(3))
print("e^{-nu'}/(1-q) + X <=", lhs.str(6))
assert lhs < 1
print("C-right PROVED for n >=", sys.argv[1] if len(sys.argv) > 1 else 100000, " margin: log R_k <= -theta*(1 - %.4f)" % float(lhs.upper()))
