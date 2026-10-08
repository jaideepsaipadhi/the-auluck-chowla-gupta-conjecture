"""Lemma D (elementary small-k region).  Usage: python3 ws_cd_D.py [N_A=100000] [gamma=1.7]
Claim: for n >= N_A and 1 <= k <= gamma n^{1/3}:  p(n,k+1) > p(n,k).
Bounds (LemmaCD.md, Lemma D1):  C(n-1,k-1)/k! <= p(n,k) <= C(n-1+k(k-1)/2, k-1)/k!.
Sufficient (Lemma D2):  Phi(k,n) := log(n-k) - log(k(k+1)) - k(k-1)^2/(2(n-k+1)) > 0.
Phi decreasing in real k >= 1 (each of the three terms is);  Phi(gamma s, s^3) >= Psi(s) :=
   log((s^2-gamma)/(gamma(gamma s+1))) - gamma^3/(2(1-gamma/s^2)),  increasing in s for s^2 > gamma+... (checked below).
Also prints, for comparison, the exact-integer check of the sufficient condition at N_A (largest k)."""
import sys
from flint import arb
from math import lgamma, log
NA = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
g = arb(sys.argv[2]) if len(sys.argv) > 2 else arb('1.7')
s = arb(NA) ** (arb(1) / 3)
Psi = ((s * s - g) / (g * (g * s + 1))).log() - g ** 3 / (2 * (1 - g / (s * s)))
print("s_A =", s.str(8), "  Psi(s_A) =", Psi.str(6)); assert Psi > 0
# monotonicity of Psi in s on [s_A, inf): d/ds log(s^2-g) = 2s/(s^2-g) > d/ds log(g s+1) = g/(g s+1) iff 2s(gs+1) > g(s^2-g), true for s>0;
# gamma^3/(2(1-g/s^2)) decreasing in s.  Hence Psi(s) >= Psi(s_A) > 0 for all n >= N_A.
# Phi(gamma s, s^3) >= Psi(s): log(n-k) = log(s^3 - g s) ; k(k+1) = g s(g s+1);
#   k(k-1)^2/(2(n-k+1)) <= k^3/(2(n-k)) = g^3 s^3/(2(s^3 - g s)) = g^3/(2(1-g/s^2)).
# float illustration: largest k with sufficient condition at a few n
def lC(a, b): return lgamma(a + 1) - lgamma(b + 1) - lgamma(a - b + 1)
for n in [NA, 10 ** 6, 10 ** 9]:
    k = 1
    while lC(n - 1, k) - log(k + 1) - lC(n - 1 + k * (k - 1) // 2, k - 1) > 0: k += 1
    print(f"n={n}: binomial criterion holds for k <= {k-1}  (gamma n^(1/3) = {float(g.mid())*n**(1/3):.1f})")
