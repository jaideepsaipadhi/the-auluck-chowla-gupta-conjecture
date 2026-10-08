"""MECH-3a: independent rigorous recomputation of TV_r = total variation of f_r(u) = u^r Li_{1-r}(e^{-u}) on [0, inf), r = 1..6.
Written from scratch (imports nothing from the project).  Arb (python-flint) throughout.
Evaluation of f_r and f_r':
  * u > US = 2.5 :  f_r = u^r Li_{1-r}(e^{-u}) via arb's own polylog;  f_r' = (r f_r - f_{r+1})/u   (identity E14, ws_mech_elem.py).
  * u in [0, US] (ball upper end <= US)  :  Bernoulli expansion.  Li_{1-r}(e^{-u}) = (-d/du)^{r-1} [1/(e^u-1)],  1/(e^u-1) = sum_n B_n u^{n-1}/n!  (B_1 = -1/2),
       f_r(u)  = sum_{n>=0} a_n u^n,  a_n = (-1)^{r-1} B_n/n! * (n-1)(n-2)...(n-r+1),
       truncated at n = NT with tail bound |B_n|/n! <= 4/(2 pi)^n (n >= 2 even; = 2 zeta(n)/(2pi)^n), |falling| <= n^{r-1}:
       |tail f|  <= sum_{n>NT} 4 n^{r-1} x^n,  |tail f'| <= sum_{n>NT} 4 n^r x^{n-1}/(2pi),  x = u/(2pi) <= 2.5/(2pi),
       each majorised by 2x its first term (term ratio <= ((NT+2)/(NT+1))^r x < 1/2, asserted).
Total variation:  adaptive cells on [0, 60].  If the arb enclosure of f' on the cell has constant sign, Var = |f(b) - f(a)|;
  otherwise bisect (to depth 14); at the bottom Var <= width * sup|f'|.
  [60, inf):  f_r = sum_l l^{-1} (lu)^r e^{-lu} is decreasing for u > r (E6), and f_r -> 0, so Var = f_r(60).
Also verifies at 200 random points that both evaluation branches agree with arb polylog (consistency only)."""
from flint import arb, fmpq
import sympy as sp, random, time
US = arb('2.5'); NT = 160; PI = arb.pi()
B = [fmpq(int(sp.bernoulli(n).p), int(sp.bernoulli(n).q)) for n in range(NT + 1)]
B[1] = fmpq(-1, 2)
def coeffs(r):
    out = []
    for n in range(NT + 1):
        fall = 1
        for j in range(1, r): fall *= (n - j)
        fac = 1
        for j in range(2, n + 1): fac *= j
        out.append(arb(B[n] * fall * (-1) ** (r - 1) / fac))
    return out
C = {r: coeffs(r) for r in range(1, 8)}
def ser(r, u, deriv):
    x = u / (2 * PI)
    xu = arb(x.abs_upper())
    rat = (arb(NT + 2) / (NT + 1)) ** r * xu
    assert rat < arb('0.5')
    if not deriv:
        val = arb(0)
        for c in reversed(C[r]): val = val * u + c
        tail = 2 * 4 * arb(NT + 1) ** (r - 1) * xu ** (NT + 1)
    else:
        val = arb(0)
        for n in range(NT, 0, -1): val = val * u + n * C[r][n]
        tail = 2 * 4 * arb(NT + 1) ** r * xu ** NT / (2 * PI)
    return val + arb(0, tail.upper())
def fr(r, u):
    if u.upper() <= US.upper() + 0: return ser(r, u, False)
    return u ** r * (-u).exp().polylog(1 - r)
def dfr(r, u):
    if u.upper() <= US.upper(): return ser(r, u, True)
    return (r * fr(r, u) - fr(r + 1, u)) / u
def var_cell(r, a, b, depth=0):
    ab = a.union(b)
    d = dfr(r, ab)
    if d > 0 or d < 0: return abs(fr(r, b) - fr(r, a))
    if depth >= 14: return (b - a) * arb(d.abs_upper())
    m = (a + b) / 2
    return var_cell(r, a, m, depth + 1) + var_cell(r, m, b, depth + 1)
if __name__ == "__main__":
    random.seed(7)
    for _ in range(200):          # consistency of the two branches near the switch and against arb polylog
        r = random.randint(1, 6); u = arb(random.uniform(0.05, 2.5))
        assert abs(ser(r, u, False) - u ** r * (-u).exp().polylog(1 - r)) < 1e-10
        assert abs(ser(r, u, True) - (r * ser(r, u, False) - ser(r + 1, u, False)) / u) < 1e-8
    OLD = {1: '1.062650', 2: '1.705253', 3: '7.708993', 4: '51.518986', 5: '236.048', 6: '1626.90'}
    t0 = time.time(); TV = {}
    for r in range(1, 7):
        N = 6000; tot = arb(0)
        for j in range(N):
            a = arb(60) * j / N; b = arb(60) * (j + 1) / N
            tot += var_cell(r, a, b)
        tot += fr(r, arb(60))
        TV[r] = tot; print(f"  [r={r} done at {time.time()-t0:.0f}s]", flush=True)
        print(f"TV_{r} <= {float(tot.upper()):.6f}   (true value enclosed in {tot.str(8)};  project constant {OLD[r]}:"
              f" {'OK, project >= new' if arb(OLD[r]) >= arb(tot.upper()) else 'PROJECT CONSTANT TOO SMALL'})", flush=True)
    print(f"time {time.time() - t0:.1f}s")
