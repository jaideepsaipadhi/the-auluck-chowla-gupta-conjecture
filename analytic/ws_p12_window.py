"""(P1-W) Rigorous (theta, nu') range of the window for all n >= N_A.  Usage: python3 ws_p12_window.py N_A [C]
beta = sqrt(6n)/pi (n = zeta2 beta^2), lb = log beta, window mu-2beta < k < mu+2beta, m = n-k, K = k+1.
nu := K theta, nu' := nu + log theta (increasing in theta at fixed K).
 S(t) := sum_{i<=k} i/(e^{it}-1) is decreasing; theta is the root of S = m.
 t^2 S(t) = t sum f1(i t), f1(u) = u/(e^u-1) decreasing, so F1(Kt) - t <= t^2 S(t) <= F1(kt) <= zeta2,
 F1(v) = int_0^v f1 = zeta2 - T(v),  T(v) <= (v+1) e^{-v}/(1-e^{-v}).
UPPER: theta <= theta_+ := sqrt(zeta2/m) <= r/beta, r := (1 - (lb+2+1/beta)/(zeta2 beta))^{-1/2}.
LOWER: theta > theta_- := (1-delta)/beta, delta = C lb/beta, provided S(theta_-) > m, which follows from
   (2delta - delta^2) zeta2 > T(K theta_-) + theta_-   (since theta_-^2 m <= (1-delta)^2 zeta2),
 equivalent (x beta/lb) to  G := 2C zeta2 - C^2 zeta2 lb/beta - e^2 beta^delta/(1-e^{-v}) - 1/lb > 0,
 v := (lb-2)(1-delta) <= K theta_-   [uses T(v) <= (lb-1) e^2 beta^{delta-1}/(1-e^{-v})].
 G is increasing in beta for beta >= e^2 (each term monotone), so checking G(beta_A) > 0 suffices.
All other outputs are monotone in beta in the favourable direction (see LemmaA_uniform.md, Lemma W)."""
import sys
from flint import arb
NA = arb(sys.argv[1]); C = arb(sys.argv[2]) if len(sys.argv) > 2 else arb(5)
z2 = arb(2).zeta(); e = arb(1).exp()
b = (6 * NA).sqrt() / arb.pi(); lb = b.log()
assert b > e ** 2
delta = C * lb / b
v = (lb - 2) * (1 - delta)
G = 2 * C * z2 - C * C * z2 * lb / b - e ** 2 * (delta * lb).exp() / (1 - (-v).exp()) - 1 / lb
print("beta_A =", b.str(8), " delta =", delta.str(5), " G =", G.str(5))
assert G > 0
eps = (lb + 2 + 1 / b) / (z2 * b)
r = (1 / (1 - eps)).sqrt()
thA = r / b
nuhi = (lb + 2 + 1 / b) * r - lb + r.log()
nulo = -2 - delta * (lb - 2) + (1 - delta).log()
print("theta in [(1-delta)/beta, r/beta], r =", r.str(6))
print("theta_A <=", float(thA.upper()), "  nu'_hi <=", float(nuhi.upper()), "  nu'_lo >=", float(nulo.lower()))
print("nu_min = (lb-2)(1-delta) >=", float(v.lower()))
