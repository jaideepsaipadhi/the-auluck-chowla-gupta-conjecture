"""CLOSE-1: sanity tests for the classical items of PROOFS_FULL.md §§1, 2, 8 (representation, derivative bounds, Region E)
and Lemma D1 (exact integers).  The PROOFS are in PROOFS_FULL.md; this script checks the statements on many instances.
Usage: python3 ws_close_rep.py          (~1-2 min; exact big integers + mpmath 40 digits)
 R1  p(n,k) = f_k(n-k) (partitions of n-k into parts <= k), f_{k+1}(N) = sum_s f_k(N - sK)        n <= 300, exact
 R2  representation (1) with EXACT RATIONAL x:  sum_s x^{1+sK} P(X=m-1-sK)/P(X=m) == p(n,k+1)/p(n,k),
     P(X=N) := f_k(N) x^N prod_{i<=k}(1-x^i)  (exact rationals; 3 values of x per (n,k))               n <= 60, all 2<=k<=n-2
 R3  P(X=.) is a probability law with mean m at the saddle (mpmath, truncated sum + explicit geometric tail bound)
 R4  representation (2) by quadrature: R/g0 - 1 = Ntil/D, D = 2pi P(X=m)                                (mpmath quad)
 R5  |psi^(r)(phi)| <= kappa_r (r=2..6), |H'| <= c1, |H^(r)| <= c_r (r=2..5), Re H <= 0 at random real phi (mpmath diff)
 E1  p(n,k) = p(n-k) for n/2 <= k <= n,  p(j) <= p(j+1),  a_{k+1} <= a_k for ceil(n/2) <= k <= n-1          n <= 400, exact
 E2  every row n <= 400 is weakly unimodal (exact)
 D1  C(n-1,k-1) <= k! p(n,k) <= C(n-1+k(k-1)/2, k-1)                                                n <= 250, exact
 D2' whenever the exact D1-criterion C(n-1,k)/(k+1)! > C(N-1,k-1)/k! (N = n+k(k-1)/2) holds, a_{k+1} > a_k   n <= 250, exact"""
import time, random
from fractions import Fraction as Fr
from math import comb, factorial
import mpmath as mp
t0 = time.time()
def ok(m): print(f"[OK] {m}  ({time.time()-t0:.0f}s)", flush=True)

NMAX = 400
# p(n,k) table by p(n,k) = p(n-1,k-1) + p(n-k,k)
p =[[0] * (NMAX + 1) for _ in range(NMAX + 1)]; p[0][0] = 1
for n in range(1, NMAX + 1):
    for k in range(1, n + 1):
        p[n][k] = p[n - 1][k - 1] + p[n - k][k]
pt = [sum(p[j][1:j + 1]) if j else 1 for j in range(NMAX + 1)]       # p(j)
# f_k(N): partitions of N into parts <= k
def fcol(k, Nmax):
    f = [1] + [0] * Nmax
    for i in range(1, k + 1):
        for N in range(i, Nmax + 1): f[N] += f[N - i]
    return f
F = {k: fcol(k, 300) for k in range(1, 302)}
for n in range(1, 301):
    for k in range(1, n + 1): assert p[n][k] == F[k][n - k]
for k in range(1, 300):
    K = k + 1
    for N in range(0, 301): assert F[k + 1][N] == sum(F[k][N - s * K] for s in range(0, N // K + 1))
ok("R1  p(n,k) = f_k(n-k) (n<=300) and f_{k+1}(N) = sum_s f_k(N-sK) (k<300, N<=300), exact")

cnt = 0
for n in range(4, 61):
    for k in range(2, n - 1):
        m = n - k; K = k + 1
        for x in (Fr(1, 3), Fr(7, 10), Fr(19, 20)):
            Z = 1
            for i in range(1, k + 1): Z *= (1 - x ** i)
            PX = lambda N: F[k][N] * x ** N * Z if N >= 0 else 0
            rhs = sum(x ** (1 + s * K) * PX(m - 1 - s * K) for s in range(0, (m - 1) // K + 1)) / PX(m)
            assert rhs == Fr(p[n][k + 1], p[n][k]); cnt += 1
ok(f"R2  representation (1) == p(n,k+1)/p(n,k) exactly, {cnt} (n,k,x) triples, x in {{1/3, 7/10, 19/20}}, n <= 60")

mp.mp.dps = 30
def saddle(m, k):
    S = lambda t: mp.fsum(i / mp.expm1(i * t) for i in range(1, k + 1)) - m
    lo, hi = mp.mpf(10) ** -6, mp.mpf(20)          # S decreasing: S(lo) > 0 > S(hi)
    for _ in range(120):
        mid = (lo + hi) / 2
        if S(mid) > 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2
for (n, k) in [(40, 6), (60, 9), (90, 12), (120, 7)]:
    m = n - k; th = saddle(m, k); x = mp.e ** (-th); K = k + 1
    Z = mp.fprod(1 - x ** i for i in range(1, k + 1))
    Nmax = 40 * m
    f = fcol(k, Nmax)
    probs = [f[N] * x ** N * Z for N in range(Nmax + 1)]
    # tail bound: f_k(N) <= C(N+k-1,k-1) and x^N decays; check the computed mass is within 1e-12 of 1 and mean ~ m
    tot = mp.fsum(probs); mean = mp.fsum(N * probs[N] for N in range(Nmax + 1))
    assert abs(tot - 1) < mp.mpf(10) ** -12 and abs(mean - m) < mp.mpf(10) ** -9, (n, k, tot, mean)
    # (2) by quadrature
    q = x ** K; g0 = x / (1 - q)
    chi = lambda ph: mp.e ** (-1j * m * ph) * mp.fprod((1 - x ** i) / (1 - x ** i * mp.e ** (1j * i * ph)) for i in range(1, k + 1))
    H = lambda ph: 1j * ph - mp.log(1 - q * mp.e ** (1j * K * ph)) + mp.log(1 - q)
    D = mp.quad(chi, mp.linspace(-mp.pi, mp.pi, 4 * k + 1))
    Nt = mp.quad(lambda ph: chi(ph) * (mp.e ** H(ph) - 1), mp.linspace(-mp.pi, mp.pi, 4 * k + 1))
    R = mp.mpf(p[n][k + 1]) / p[n][k]
    assert abs(D - 2 * mp.pi * probs[m]) < mp.mpf(10) ** -20
    assert abs((R / g0 - 1) - Nt / D) < mp.mpf(10) ** -20, (n, k)
    # (1) also numerically with the irrational saddle x (sanity of P(X=N) normalisation)
    r1 = mp.fsum(x ** (1 + s * K) * probs[m - 1 - s * K] for s in range(0, (m - 1) // K + 1)) / probs[m]
    assert abs(r1 - R) < mp.mpf(10) ** -20
    # R5 derivative bounds at random phi
    kap = {r: mp.fsum(mp.mpf(i) ** r * mp.polylog(1 - r, x ** i) for i in range(1, k + 1)) for r in range(2, 7)}
    cc = {1: 1 + K * q / (1 - q)}; cc.update({r: mp.mpf(K) ** r * mp.polylog(1 - r, q) for r in range(2, 6)})
    psi = lambda ph: -1j * m * ph + mp.fsum(mp.log(1 - x ** i) - mp.log(1 - x ** i * mp.e ** (1j * i * ph)) for i in range(1, k + 1))
    assert abs(mp.diff(psi, 0, 1)) < mp.mpf(10) ** -15                 # E X = m  <=>  psi'(0) = 0
    for r in range(2, 7): assert abs(mp.diff(psi, 0, r) - (1j) ** r * kap[r]) < mp.mpf(10) ** -10 * kap[r]
    random.seed(n * 1000 + k)
    for _ in range(6):
        ph = mp.mpf(random.uniform(-mp.pi, mp.pi))
        for r in range(2, 7): assert abs(mp.diff(psi, ph, r)) <= kap[r] * (1 + mp.mpf(10) ** -10)
        assert abs(mp.diff(H, ph, 1)) <= cc[1] * (1 + mp.mpf(10) ** -10)
        for r in range(2, 6): assert abs(mp.diff(H, ph, r)) <= cc[r] * (1 + mp.mpf(10) ** -10)
        assert mp.re(H(ph)) <= 0
        assert abs(mp.e ** psi(ph) - chi(ph)) < mp.mpf(10) ** -20        # chi = e^psi with factorwise principal logs
ok("R3/R4/R5  law of X (mass 1, mean m at saddle), (2) by quadrature, psi^(r)(0)=i^r kappa_r, |psi^(r)|<=kappa_r, |H^(r)|<=c_r, Re H<=0, chi=e^psi")

for n in range(1, NMAX + 1):
    for k in range(1, n + 1):
        if 2 * k >= n: assert p[n][k] == pt[n - k]
    for k in range((n + 1) // 2, n): assert p[n][k + 1] <= p[n][k] and (p[n][k + 1] < p[n][k] or k == n - 1)
for j in range(NMAX): assert pt[j] <= pt[j + 1]
ok("E1  p(n,k) = p(n-k) for 2k >= n; p(j) <= p(j+1); a_{k+1} <= a_k on ceil(n/2) <= k <= n-1, strict except k = n-1 (n <= 400)")
for n in range(1, NMAX + 1):
    row = p[n][1:n + 1]; i = 0
    while i + 1 < n and row[i + 1] >= row[i]: i += 1
    while i + 1 < n and row[i + 1] <= row[i]: i += 1
    assert i == n - 1, n
ok("E2  rows n <= 400 weakly unimodal (exact)")

for n in range(2, 251):
    for k in range(1, n + 1):
        D_ = k * (k - 1) // 2
        assert comb(n - 1, k - 1) <= factorial(k) * p[n][k] <= comb(n - 1 + D_, k - 1)
ok("D1  C(n-1,k-1) <= k! p(n,k) <= C(n-1+k(k-1)/2, k-1), all 1<=k<=n<=250 (exact)")
c = 0
for n in range(3, 251):
    for k in range(1, n - 1):
        N_ = n + k * (k - 1) // 2
        if comb(n - 1, k) * factorial(k) > factorial(k + 1) * comb(N_ - 1, k - 1):    # lower(a_{k+1}) > upper(a_k)
            assert p[n][k + 1] > p[n][k]; c += 1
ok(f"D2' exact binomial criterion => a_{{k+1}} > a_k in all {c} instances (n <= 250)")
print(f"all checks passed ({time.time()-t0:.0f}s)")
