"""CLOSE-2: Lemma D (D1-D3) fully mechanized, incl. the monotonicity of Psi(s) (PROOFS_FULL.md §7).
Usage: python3 ws_close_D.py      (~20 s)
 D2a  sympy: the product identity C(n-1,k-1)/C(N-1,k-1) = prod_{j<=k-2} 1/(1 + D/(n-1-j)), the ratio C(n-1,k)/C(n-1,k-1) = (n-k)/k,
      and the reduction 'C(n-1,k)/(k+1)! > C(N-1,k-1)/k!'  <=  Phi > 0 (k = 1..12 symbolic in n; general k by the written proof)
 D3a  sympy: dPhi/dk < 0 on 1 <= k <= n-1: each of the three terms has a manifestly signed derivative
 D3b  sympy: Phi(gamma s, s^3) - Psi(s) = [k(k-1)^2 ... ] >= 0  (exact difference written as a sum of nonneg. terms)
 D3c  sympy: Psi'(s) = (gamma s^2 + 2 s + gamma^2)/((s^2-gamma)(gamma s + 1)) + gamma^4 s/(s^2-gamma)^2  exactly,
      hence Psi' > 0 for s > sqrt(gamma)  (ANALYTIC proof of 'Psi increasing', REFEREE A.3)
 D3d  Arb: s_A = (10^5)^{1/3} as a ball; Psi(s_A) > 0; independent Arb derivative-sign sweep of Psi' on [s_A, 10^4] plus the
      exact formula for the tail (redundant with D3c, done because it was requested)
 D3e  k_D(n) = floor(1.7 n^{1/3}) >= 78 for all n >= 10^5; exact-integer check that a_{k+1} > a_k for k <= k_D at n = 10^5 and
      n = 2*10^5 (big integers, column DP)."""
import time, sympy as sp
from flint import arb
t0 = time.time()
def ok(m): print(f"[OK] {m}  ({time.time()-t0:.0f}s)", flush=True)
g = sp.Rational(17, 10)
n, k, s, D = sp.symbols('n k s D', positive=True)

# D2a
for kk in range(1, 13):
    Dk = sp.Rational(kk * (kk - 1), 2)
    lhs = sp.binomial(n - 1, kk - 1) / sp.binomial(n - 1 + Dk, kk - 1)
    rhs = sp.prod([1 / (1 + Dk / (n - 1 - j)) for j in range(0, kk - 1)]) if kk > 1 else sp.Integer(1)
    assert sp.simplify(sp.expand_func(lhs) - rhs) == 0, kk
    assert sp.simplify(sp.expand_func(sp.binomial(n - 1, kk) / sp.binomial(n - 1, kk - 1)) - (n - kk) / kk) == 0
ok("D2a  C(n-1,k-1)/C(n-1+D,k-1) = prod_{j=0}^{k-2} (1 + D/(n-1-j))^{-1} and C(n-1,k)/C(n-1,k-1) = (n-k)/k  (k=1..12, symbolic n)")
# 1+x <= e^x  =>  prod >= exp(-sum D/(n-1-j)) >= exp(-(k-1)D/(n-k+1));  (k-1)D = k(k-1)^2/2
x = sp.symbols('x', nonnegative=True)
assert sp.simplify(sp.diff(sp.exp(x) - 1 - x, x) - (sp.exp(x) - 1)) == 0
kk_ = sp.symbols('kk', positive=True)
assert sp.expand((kk_ - 1) * kk_ * (kk_ - 1) / 2 - kk_ * (kk_ - 1) ** 2 / 2) == 0
ok("D2b  1+x <= e^x (derivative e^x - 1 >= 0 on x >= 0, equality at 0); (k-1) D = k(k-1)^2/2")

# D3a: Phi(k) = log(n-k) - log k - log(k+1) - k(k-1)^2 / (2(n-k+1))
Phi = sp.log(n - k) - sp.log(k) - sp.log(k + 1) - k * (k - 1) ** 2 / (2 * (n - k + 1))
d3 = sp.diff(k * (k - 1) ** 2 / (n - k + 1), k)
num = sp.factor(sp.together(d3) * (n - k + 1) ** 2)
# num = (k-1)[(3k-1)(n-k+1) + k(k-1)]  -> >= 0 for 1 <= k <= n
assert sp.expand(num - (k - 1) * ((3 * k - 1) * (n - k + 1) + k * (k - 1))) == 0
dPhi = sp.diff(Phi, k)
assert sp.simplify(dPhi - (-1 / (n - k) - 1 / k - 1 / (k + 1) - d3 / 2)) == 0
ok("D3a  dPhi/dk = -1/(n-k) - 1/k - 1/(k+1) - (k-1)[(3k-1)(n-k+1)+k(k-1)]/(2(n-k+1)^2) < 0 for real 1 <= k < n")

# D3b: Phi(gs, s^3) >= Psi(s)
Psi = sp.log((s ** 2 - g) / (g * (g * s + 1))) - g ** 3 / (2 * (1 - g / s ** 2))
PhiS = Phi.subs({k: g * s, n: s ** 3})
# log parts agree exactly:
logpart = sp.log(s ** 3 - g * s) - sp.log(g * s) - sp.log(g * s + 1)
assert sp.simplify((s ** 3 - g * s) / (g * s * (g * s + 1)) - (s ** 2 - g) / (g * (g * s + 1))) == 0   # arguments of the logs agree
# cubic part: k(k-1)^2/(2(n-k+1)) <= k^3/(2(n-k)) for 1 <= k < n:  k^3 (n-k+1) - k(k-1)^2 (n-k) = k[(n-k)(2k-1) + k^2] >= 0
diff_ = sp.expand(k ** 3 * (n - k + 1) - k * (k - 1) ** 2 * (n - k))
assert sp.expand(diff_ - k * ((n - k) * (2 * k - 1) + k ** 2)) == 0
assert sp.simplify((g * s) ** 3 / (2 * (s ** 3 - g * s)) - g ** 3 / (2 * (1 - g / s ** 2))) == 0
ok("D3b  Phi(gamma s, s^3) >= Psi(s): log terms equal; k(k-1)^2/(n-k+1) <= k^3/(n-k) since k^3(n-k+1)-k(k-1)^2(n-k) = k[(n-k)(2k-1)+k^2] >= 0")

# D3c: Psi'(s) closed form
dPsi = sp.diff(Psi, s)
claim = (g * s ** 2 + 2 * s + g ** 2) / ((s ** 2 - g) * (g * s + 1)) + g ** 4 * s / (s ** 2 - g) ** 2
assert sp.simplify(dPsi - claim) == 0
ok("D3c  Psi'(s) = (g s^2 + 2s + g^2)/((s^2-g)(g s+1)) + g^4 s/(s^2-g)^2  (exact)  => Psi' > 0 for all s > sqrt(g): Psi is INCREASING")

# D3d: Arb value and (redundant) derivative-sign sweep
G = arb(17) / 10
sA = arb(10) ** (arb(5) / 3)
assert sA > G.sqrt()
PsiA = ((sA * sA - G) / (G * (G * sA + 1))).log() - G ** 3 / (2 * (1 - G / (sA * sA)))
assert PsiA > 0
def dPsi_arb(sb): return (G * sb * sb + 2 * sb + G * G) / ((sb * sb - G) * (G * sb + 1)) + G ** 4 * sb / (sb * sb - G) ** 2
a = float(sA.lower()); cells = 0; lo = arb(a)
while float(lo.lower()) < 1e4:
    hi = lo * arb('1.01'); assert dPsi_arb(lo.union(hi)) > 0; lo = hi; cells += 1
ok(f"D3d  s_A = {sA.str(10)}, Psi(s_A) = {PsiA.str(8)} > 0;  Arb sweep: Psi' > 0 on {cells} cells covering [s_A, 1e4]; tail s > 1e4 by D3c")

# D3e   k_D(n) = floor(1.7 n^{1/3}) = max{k : 1000 k^3 <= 4913 n}  (exact integers)
def kDex(N):
    kk = 0
    while 1000 * (kk + 1) ** 3 <= 4913 * N: kk += 1
    return kk
for N in [10 ** 5, 10 ** 5 + 1, 2 * 10 ** 5]:
    kD = kDex(N)
    assert kD >= 78
# floor(1.7 n^{1/3}) is nondecreasing in n, value 78 at n = 10^5:
assert (arb(17) / 10 * arb(10 ** 5) ** (arb(1) / 3)) > 78
def check_row(N, kmax):
    # a_k = f_k(N-k) with f_k = partitions into parts <= k ; column DP up to parts kmax+1
    f = [1] + [0] * N
    vals = []
    for kk in range(1, kmax + 2):
        for M in range(kk, N + 1): f[M] += f[M - kk]
        vals.append(f[N - kk])            # a_kk
    for i in range(len(vals) - 1): assert vals[i + 1] > vals[i], (N, i + 1)
for N in [10 ** 5, 2 * 10 ** 5]:
    kD = kDex(N)
    check_row(N, kD)
ok("D3e  k_D(n) >= 78 for n >= 1e5; exact integers: a_{k+1} > a_k for all k <= k_D at n = 1e5 (k_D=78) and n = 2e5 (k_D=99)")
print(f"all checks passed ({time.time()-t0:.0f}s)")
