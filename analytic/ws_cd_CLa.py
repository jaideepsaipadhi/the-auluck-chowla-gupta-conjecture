"""(C-left, small v) Lemma CL-a.  Usage: python3 ws_cd_CLa.py K0 V0
Claim: for every n and every k with k >= K0 and v := k*theta <= V0 (theta = saddle of (m=n-k, k)),  R_k > 1.
No dependence on n.  Proof: Lemma C0 with phi0 = c0*theta; all ingredients bounded uniformly in (v <= V0, k >= K0):
  theta = v/k <= V0/K0 =: th0,  t := 1/k <= 1/K0
  b = s_2 in [v f2(v), v] (f2 decreasing, f2 <= 1),  s3 <= v*sup f_3,  s4 <= v*sup f_4,  c1*theta <= 1 + th0,  c2 theta^2 <= 1
  T1 = theta*(c1 theta) s3/(2b^2) <= t (1+th0) MS3 / (2 f2(V0)^2)
  T2 = theta*(f2 + (c1 theta)^2)/(2b) <= t (1 + (1+th0)^2)/(2 f2(V0))
  A3 = s3 theta^{1/2} b^{-3/2} <= MS3 t^{1/2}/f2(V0)^{3/2},  w0 = c0 (b/theta)^{1/2} >= c0 (K0 f2(V0))^{1/2}
  delta = min(c0^2 MS4/(12 f2(V0)),  1 - log(1+c0^2)/c0^2 * (1 - (c0 V0)^2/12))
 minor arcs (normalised N := sqrt(B/2pi) int_{phi0<=|phi|<=pi}|chi|), N <= N1+N2+N3+N4:
  N1: phi in [c_j theta, c_{j+1} theta] (c0 = c_0 < ... < c_J = CA):  |chi| <= exp(-lam_j B phi^2/2),
      lam_j = log(1+s_j^2 c_{j+1}^2)/c_{j+1}^2, s_j = sin(T_j)/T_j, T_j = c_{j+1} V0/2 (< pi/2);
      N1 <= sum_j erfc(c_j sqrt(lam_j K0 f2(V0)/2))/sqrt(lam_j)                       (theta^2 B = b/theta >= k f2(v))
  N2: phi in [CA theta, pi/k] (s = 2/pi):  S >= k log(1 + s^2 c^2 f2(v)),
      N2 <= sqrt(2/pi) sqrt(k) (s^2 f2(V0))^{-k/2} CA^{1-k}/(k-1)                          (theta sqrt B = sqrt(b/theta) <= sqrt k)
  N3: phi in [pi/k, PB]: block counting, #good >= k*gam - 2, gam = 1 - (2 al + PB/2)/pi - (1+al/pi) 4 al/pi;
      N3 <= sqrt(2 pi) k^{3/2} v^{-1} (v sig0/(2 sin al))^{(k gam - 2)/2}                   (sqrt B <= k^{3/2}/v)
  N4: phi in [PB, pi]: pairs:  N4 <= sqrt(2 pi) k^{3/2} v^{-1} (v sig0/(2 sin(PB/4)))^{floor(k/2)/2}
      sig0 = sinh(V0/2)/(V0/2) >= sinh(v/2)/(v/2).   N3,N4 increasing in v (exponents >= 2), decreasing in k (checked).
 rho = (T1/(1-d)^{5/2} + T2/(1-d)^{3/2} + 2N)/(1 - erfc(w0/sqrt2) - sqrt(2/pi)A3/(3(1-d)^2) - N),
 R_k/g0 >= 1 - rho and log R_k >= -theta - log(1-q) + log(1-rho) >= -log(nu e^theta) + log(1-rho),  nu = v+theta <= V0(1+1/K0).
 Need 1 - rho > e^{th0} V0 (1 + 1/K0)."""
import sys
from flint import arb
from ws_cd_common import f2, MS, PI, SQ2PI
K0 = int(sys.argv[1]) if len(sys.argv) > 1 else 79
V0 = arb(sys.argv[2]) if len(sys.argv) > 2 else arb('0.1')
PB = arb('0.6'); CA = arb(10)
th0 = V0 / K0; t = arb(1) / K0
F2 = f2(V0); sig0 = (V0 / 2).sinh() / (V0 / 2)
best = None
for c0 in [arb(x) / 10 for x in range(3, 13)]:
    T1 = t * (1 + th0) * MS[3] / (2 * F2 ** 2)
    T2 = t * (1 + (1 + th0) ** 2) / (2 * F2)
    A3 = MS[3] * t.sqrt() / F2 ** arb(1.5)
    w0 = c0 * (K0 * F2).sqrt()
    delta = (c0 ** 2 * MS[4] / (12 * F2)).min(1 - (1 + c0 ** 2).log() / c0 ** 2 * (1 - (c0 * V0) ** 2 / 12))
    # N1
    J = 60; r = (CA / c0) ** (arb(1) / J); N1 = arb(0)
    for j in range(J):
        ca = c0 * r ** j; cb = c0 * r ** (j + 1)
        T = cb * V0 / 2; assert T < PI / 2
        s = T.sin() / T
        lam = (1 + s * s * cb * cb).log() / (cb * cb)
        N1 += (ca * (lam * K0 * F2 / 2).sqrt()).erfc() / lam.sqrt()
    s2 = (2 / PI) ** 2
    def N2f(k): return SQ2PI * arb(k).sqrt() * (s2 * F2) ** (-arb(k) / 2) * CA ** (1 - k) / (k - 1)
    y3 = V0 * sig0 / (2 * arb('0.25').sin()); y4 = V0 * sig0 / (2 * (PB / 4).sin())
    al = arb('0.25'); gam = 1 - (2 * al + PB / 2) / PI - (1 + al / PI) * 4 * al / PI
    def N3f(k, v): return (2 * PI).sqrt() * arb(k) ** arb(1.5) / v * (v * sig0 / (2 * al.sin())) ** ((k * gam - 2) / 2)
    def N4f(k, v): return (2 * PI).sqrt() * arb(k) ** arb(1.5) / v * (v * sig0 / (2 * (PB / 4).sin())) ** (arb(k // 2) / 2)
    assert (K0 * gam - 2) / 2 >= 1 and K0 // 2 >= 2 and y3 < 1 and y4 < 1   # => N3, N4 increasing in v on (0, V0]
    # decreasing in k for k >= K0: log N2 slope: 1/(2k) - log(...)/... check via ratios N(k+1)/N(k) < 1 bound:
    x2 = 1 / (CA * (s2 * F2).sqrt()); assert x2 < 1
    #   N2(k+1)/N2(k) = sqrt((k+1)/k) (k-1)/k * x2 <= sqrt(2) x2 < 1
    assert arb(2).sqrt() * x2 < 1
    #   N3(k+1)/N3(k) = ((k+1)/k)^{1.5} y3^{gam/2};  N4(k+2)/N4(k) = ((k+2)/k)^{1.5} y4^{1/2}; N4(k+1)/N4(k) <= (1+1/K0)^1.5 (taken care of by max over K0, K0+1)
    assert (1 + arb(1) / K0) ** arb(1.5) * y3 ** (gam / 2) < 1 and (1 + arb(2) / K0) ** arb(1.5) * y4.sqrt() < 1
    N4max = N4f(K0, V0).max(N4f(K0 + 1, V0))
    N = N1 + N2f(K0) + N3f(K0, V0) + N4max
    if not delta < 1: continue
    den = 1 - (w0 / arb(2).sqrt()).erfc() - SQ2PI * A3 / (3 * (1 - delta) ** 2) - N
    if not den > 0: continue
    rho = (T1 / (1 - delta) ** arb(2.5) + T2 / (1 - delta) ** arb(1.5) + 2 * N) / den
    if best is None or rho.upper() < best[0].upper(): best = (rho, c0, N1, N2f(K0), N3f(K0, V0), N4max, delta)
rho, c0, N1, N2, N3, N4, delta = best
need = 1 - th0.exp() * V0 * (1 + arb(1) / K0)
print(f"K0={K0} V0={V0.str(4)}  best c0={c0.str(3)}  delta={delta.str(4)}")
print("  N1..N4 =", [z.str(3) for z in (N1, N2, N3, N4)])
print("  rho <=", rho.str(6), "   need rho < 1 - e^th0 V0 (1+1/K0) =", need.str(6))
assert rho < need
print("CL-a PROVED: R_k > 1 whenever k >=", K0, "and k*theta <=", V0.str(4))
