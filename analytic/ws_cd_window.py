"""(C-W) Rigorous parameter ranges for regions CL and CR, all n >= N_A.  Usage: python3 ws_cd_window.py N_A C
beta = sqrt(6n)/pi, lb = log beta, delta = C lb/beta.  Facts used (LemmaCD.md, Lemma CW):
 (i)  theta_k is increasing in k (S(t;k) increasing in k, m = n-k decreasing, S decreasing in t).
 (ii) theta^2 m <= zeta2 (t^2 S(t) <= F1(kt) <= zeta2), hence theta <= sqrt(zeta2/m).
 (iii) for k >= mu - 2 beta: theta > (1-delta)/beta, valid when G > 0 (proof of Lemma W, LemmaA_uniform.md;
      the proof only uses m <= n and K theta_- >= (lb-2)(1-delta)).
CL (k <= mu-2beta): theta <= r/beta (by (i) and Lemma W), nu' = K theta + log theta <= nuCL := (lb-2+1/beta) r - lb + log r.
CR (mu+2beta <= k <= ceil(n/2)): m >= (n-1)/2 so theta <= thCR := sqrt(2 zeta2/(n-1));
     nu' >= nuCR := (lb+2)(1-delta) - lb + log(1-delta)  (nu' increasing in theta and K; K >= mu+2beta).
Monotonicity in beta (beta >= beta_A): nuCL = -2 + 1/beta + (lb-2+1/beta)(r-1) + log r, with r-1 = O(lb/beta) decreasing => nuCL decreasing;
 nuCR = 2 - delta(lb+2) + log(1-delta): delta(lb+2) = C lb(lb+2)/beta and delta decreasing for lb >= 2 => nuCR increasing.
 thA = r/beta and thCR decreasing in n.  So the values at N_A are valid for all n >= N_A."""
import sys
from flint import arb
NA = arb(sys.argv[1]); C = arb(sys.argv[2]) if len(sys.argv) > 2 else arb(5)
z2 = arb(2).zeta(); e = arb(1).exp()
b = (6 * NA).sqrt() / arb.pi(); lb = b.log()
delta = C * lb / b
v = (lb - 2) * (1 - delta)
G = 2 * C * z2 - C * C * z2 * lb / b - e ** 2 * (delta * lb).exp() / (1 - (-v).exp()) - 1 / lb
assert G > 0, "Lemma W lower bound fails for this C"
eps = (lb + 2 + 1 / b) / (z2 * b); r = (1 / (1 - eps)).sqrt()
thA = r / b
nuCL = (lb - 2 + 1 / b) * r - lb + r.log()
nuCR = (lb + 2) * (1 - delta) - lb + (1 - delta).log()          # K >= mu + 2 beta
thCR = (2 * z2 / (NA - 1)).sqrt()
print("G =", G.str(5))
print("thA  <=", float(thA.upper()))
print("nuCL <=", float(nuCL.upper()))
print("thCR <=", float(thCR.upper()))
print("nuCR >=", float(nuCR.lower()))
