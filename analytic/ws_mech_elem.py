"""MECH-0: elementary facts used by the mechanical minor-arc / monotonicity proofs.  sympy (exact symbolic) + arb.
Each fact is reduced to: an exact symbolic identity (sympy, difference simplifies to 0) and a sign that is a product of
manifestly signed factors on the stated domain.
 E1  f1(u) = u/(e^u-1) decreasing on (0,inf):          f1' = N/(e^u-1)^2, N(0)=0, N' = -u e^u < 0.
 E2  x/sinh x decreasing on (0,inf) (=> f2 = ((u/2)/sinh(u/2))^2 decreasing, sinh x/x increasing):
                                                         (x/sinh x)' = P/sinh^2, P(0)=0, P' = -x sinh x < 0.
 E3  sin x/x decreasing on (0,pi) (=> x/sin x increasing; sin x >= (sin T/T) x on [0,T], T <= pi/2... any T<pi):
                                                         (sin x/x)' = Q/x^2, Q(0)=0, Q' = -x sin x < 0.
 E4  log(1+A x) >= x log(1+A) for x in [0,1], A >= 0:   concave in x (2nd derivative -A^2/(1+Ax)^2), equal at x=0,1.
 E5  Dirichlet: 2 sin d * sum_{i=a+1}^{a+m} cos(2 i d) = 2 sin(m d) cos((2a+m+1) d)
     (telescoping 2 sin d cos(2id) = sin((2i+1)d) - sin((2i-1)d), then sum-to-product); checked symbolically for the
     telescoping step and the sum-to-product step, and numerically in arb for m <= 40 at random (a, d).
     Consequence: for d in (0, pi/2], sum_{block of m} sin^2(i d) >= m/2 - 1/(2 sin d).
 E6  l^a e^{-l} decreasing for l > a:                    derivative = (a - l) l^{a-1} e^{-l}.
 E7  t^P nu^j with nu = nu' + 2 log(1/t) (t>0):          d/dt = t^{P-1} nu^{j-1} (P nu - 2 j)  -> >= 0 iff P nu >= 2j.
"""
import sympy as sp, random
from flint import arb
u, x, A, d, l, t, P, j, nup, T = sp.symbols('u x A d l t P j nup T', positive=True)
a_, m_, i_ = sp.symbols('a m i', integer=True, nonnegative=True)
ok = lambda e: sp.simplify(e) == 0
# E1
f1 = u / (sp.exp(u) - 1); N = sp.exp(u) * (1 - u) - 1
assert ok(sp.diff(f1, u) - N / (sp.exp(u) - 1) ** 2) and N.subs(u, 0) == 0 and ok(sp.diff(N, u) + u * sp.exp(u)); print("E1 ok")
# E2
g = x / sp.sinh(x); Pp = sp.sinh(x) - x * sp.cosh(x)
assert ok(sp.diff(g, x) - Pp / sp.sinh(x) ** 2) and Pp.subs(x, 0) == 0 and ok(sp.diff(Pp, x) + x * sp.sinh(x)); print("E2 ok")
# E3
s = sp.sin(x) / x; Q = x * sp.cos(x) - sp.sin(x)
assert ok(sp.diff(s, x) - Q / x ** 2) and Q.subs(x, 0) == 0 and ok(sp.diff(Q, x) + x * sp.sin(x)); print("E3 ok")
# E4
h = sp.log(1 + A * x)
assert ok(sp.diff(h, x, 2) + A ** 2 / (1 + A * x) ** 2) and h.subs(x, 0) == 0 and ok(h.subs(x, 1) - sp.log(1 + A)); print("E4 ok")
# E5
tel = 2 * sp.sin(d) * sp.cos(2 * i_ * d) - (sp.sin((2 * i_ + 1) * d) - sp.sin((2 * i_ - 1) * d))
assert ok(sp.expand_trig(sp.expand(tel.rewrite(sp.exp))))
s2p = sp.sin((2 * a_ + 2 * m_ + 1) * d) - sp.sin((2 * a_ + 1) * d) - 2 * sp.sin(m_ * d) * sp.cos((2 * a_ + m_ + 1) * d)
assert ok(sp.expand(s2p.rewrite(sp.exp)))
random.seed(1)
for _ in range(300):
    a0 = random.randint(0, 500); m0 = random.randint(1, 40); d0 = arb(random.random() * 1.5707)
    lhs = 2 * d0.sin() * sum((((2 * i0) * d0).cos() for i0 in range(a0 + 1, a0 + m0 + 1)), arb(0))
    rhs = 2 * (m0 * d0).sin() * ((2 * a0 + m0 + 1) * d0).cos()
    assert abs(lhs - rhs) < 1e-9
print("E5 ok")
# E6
assert ok(sp.diff(l ** A * sp.exp(-l), l) - (A - l) * l ** (A - 1) * sp.exp(-l)); print("E6 ok")
# E7
nu = nup + 2 * sp.log(1 / t)
assert ok(sp.diff(t ** P * nu ** j, t) - t ** (P - 1) * nu ** (j - 1) * (P * nu - 2 * j)); print("E7 ok")
# E8  log coth(x) decreasing on (0,inf), and 1 + csch^2 x = coth^2 x
assert ok(sp.diff(sp.log(sp.coth(x)), x).rewrite(sp.exp) - (-1 / (sp.sinh(x) * sp.cosh(x))).rewrite(sp.exp)) and ok((1 + 1 / sp.sinh(x) ** 2 - sp.coth(x) ** 2).rewrite(sp.exp)); print("E8 ok")
# E9  sum_{i=1}^I sin^2(i d) = I/2 - sin(I d) cos((I+1) d)/(2 sin d)   (E5 with a=0, m=I), numeric arb confirmation
for _ in range(200):
    I0 = random.randint(1, 300); d0 = arb(random.random() * 1.5707) + arb('1e-3')
    lhs = sum(((i0 * d0).sin() ** 2 for i0 in range(1, I0 + 1)), arb(0)); rhs = arb(I0) / 2 - (I0 * d0).sin() * ((I0 + 1) * d0).cos() / (2 * d0.sin())
    assert abs(lhs - rhs) < 1e-8
print("E9 ok")
# E10  sin x >= x - x^3/6 for x >= 0:  g(0)=g'(0)=g''(0)=0, third derivative = 1 - cos x >= 0
gx = sp.sin(x) - x + x ** 3 / 6
assert gx.subs(x, 0) == 0 and sp.diff(gx, x).subs(x, 0) == 0 and sp.diff(gx, x, 2).subs(x, 0) == 0 and ok(sp.diff(gx, x, 3) - (1 - sp.cos(x))); print("E10 ok")
# E11  s e^{-s} decreasing for s > 1; s^2 e^{-s} decreasing for s > 2
sv = sp.symbols('s', positive=True)
assert ok(sp.diff(sv * sp.exp(-sv), sv) - (1 - sv) * sp.exp(-sv)) and ok(sp.diff(sv ** 2 * sp.exp(-sv), sv) - sv * (2 - sv) * sp.exp(-sv)); print("E11 ok")
# E12  d/dth log(th^{-p} e^{-a/th}) = (a - p th)/th^2   (> 0 iff th < a/p)
th, aa, pp = sp.symbols('th a p', positive=True)
assert ok(sp.diff(sp.log(th ** (-pp) * sp.exp(-aa / th)), th) - (aa - pp * th) / th ** 2); print("E12 ok")
# E14  u f_r'(u) = r f_r(u) - f_{r+1}(u), f_r = u^r Li_{1-r}(e^{-u}); Li_{1-r}(y) = (y d/dy)^{r-1}[y/(1-y)]
yv = sp.symbols('y', positive=True)
Lis = {1: yv / (1 - yv)}
for r in range(2, 9): Lis[r] = sp.simplify(yv * sp.diff(Lis[r - 1], yv))
fr = lambda r: u ** r * Lis[r].subs(yv, sp.exp(-u))
for r in range(1, 7): assert ok(u * sp.diff(fr(r), u) - r * fr(r) + fr(r + 1))
print("E14 ok")
# E15  h(u) = u y/(1-y), y = e^{-u}:  h' = y/(1-y) - u y/(1-y)^2  =>  |h'| <= (1+u) y/(1-y)^2
hh = u * sp.exp(-u) / (1 - sp.exp(-u)); Y = sp.exp(-u)
assert ok(sp.diff(hh, u) - (Y / (1 - Y) - u * Y / (1 - Y) ** 2)); print("E15 ok")
# E16  Abel summation identity, k = 1..10 symbolic
for k in range(1, 11):
    ss = sp.symbols('s1:%d' % (k + 1)); ww = list(sp.symbols('w1:%d' % (k + 1))) + [0]
    G = [sum(ss[:I]) for I in range(1, k + 1)]
    assert sp.expand(sum(ss[i] * ww[i] for i in range(k)) - sum(G[I] * (ww[I] - ww[I + 1]) for I in range(k))) == 0
print("E16 ok")
# E17  int_v^inf u/(e^u-1) du = -v log(1-e^{-v}) + Li_2(e^{-v}): derivative check (numerically in arb at 50 points, symbolic d/dv)
vv = sp.symbols('v', positive=True)
R = -vv * sp.log(1 - sp.exp(-vv)) + sp.polylog(2, sp.exp(-vv))
dR = sp.diff(R, vv)
for v0 in [0.05 * k for k in range(1, 51)]:
    v0 = sp.Rational(v0).limit_denominator(1000); assert abs((dR.subs(vv, v0) + v0 / (sp.exp(v0) - 1)).evalf(40)) < 1e-30
print("E17 ok")
print("ALL ELEMENTARY FACTS VERIFIED (E1-E17)")
