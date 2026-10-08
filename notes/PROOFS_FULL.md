# PROOFS_FULL — written proofs of every step of the ACG proof that is not machine-derived

Date: 2026-10-08. This file closes items 1–8 of the trust list in STATUS.md (as of the morning of 2026-10-08).
For each item, it gives either a complete written proof or a pointer to an exact or rigorous machine check. In most
cases it gives both. The style rule is: every inequality is justified, and nothing is called "clear" or "easy".

Machine checks introduced here (all run in < 1 min; outputs are next to the scripts):

| script | what it checks | kind |
|---|---|---|
| ws_close_rep.py → .out | (1) exactly with rational x; (2) by quadrature; the law of X; derivative bounds at random φ; Region E; Lemma D1; rows n ≤ 400 | exact integers / mpmath sanity |
| ws_close_D.py → .out | Lemma D2–D3 symbolically, including **Ψ′(s) > 0 (exact closed form)**; Arb Ψ(s_A); exact a_{k+1} > a_k for k ≤ k_D at n = 10⁵, 2·10⁵ | sympy + Arb + exact |
| ws_close_elem.py → .out | the [ELEM] Taylor identities; the F_r telescoping and its limit at 0; Arb quadrature of F_r; the B.1/B.3/B.4/T₂ identities; Region II majorant facts; the C0 Gaussian integrals; the C1 scaled forms; the CL-a normalisations; B.1–B.4 and the hull at real saddles | sympy (S) / Arb (A) / numeric sanity (N) |
| ws_close_thresh.py → .out | every lemma's hypotheses hold for all n ≥ 10⁵, and the regions cover 1 ≤ k ≤ n−1 (reads the certified numbers from the .out files) | bookkeeping + Arb |

Notation follows LemmaA.md, LemmaA_uniform.md, LemmaB.md and LemmaCD.md.
- a_k = p(n,k), m = n−k, K = k+1.
- θ is the saddle, x = e^{−θ}, q = x^K, ν = Kθ, ν′ = ν + log θ, v = kθ.
- f_r(u) = u^r Li_{1−r}(e^{−u}). F_r(v) = ∫₀^v f_r.
- s_r = θ^{r+1}κ_r and b = s₂.

All logarithms of complex numbers are principal logarithms of numbers with positive real part. Where the variable is
φ ∈ ℝ, "Taylor's theorem with integral remainder" means: for a C^{N+1} function g: ℝ → ℂ,

  g(φ) − Σ_{r≤N} g^{(r)}(0)φ^r/r! = ∫₀^φ g^{(N+1)}(u)(φ−u)^N/N! du,

so |g(φ) − T_N g(φ)| ≤ sup|g^{(N+1)}|·|φ|^{N+1}/(N+1)!.
This is the real Taylor theorem applied to Re g and Im g separately, then recombined; the integral form holds for
complex-valued g verbatim. The integral ∫₀^p (p−u)^N/N! du = p^{N+1}/(N+1)! is checked in ws_close_elem X2.

---

## 1. Representation (1)/(2) (LemmaA.md §1) — trust item 1

**Setting.** Fix n and 2 ≤ k ≤ n−2. Put m = n−k ≥ 2 and K = k+1.

**1.1 Combinatorial identity.** Let f_k(N) be the number of partitions of N into parts of size ≤ k, with f_k(N) = 0 for N < 0.
- *Claim A:* p(n,k) = f_k(n−k).
  *Proof.* Take a partition λ of n into exactly k parts. Its conjugate λ′ has largest part exactly k. Remove one part equal to k
  from λ′; what remains is a partition of n−k into parts ≤ k. Conversely, adjoin a part k. The two maps are inverse to each other. □
- *Claim B:* f_{k+1}(N) = Σ_{s≥0} f_k(N − sK).
  *Proof.* Classify the partitions of N into parts ≤ K by the number s of parts equal to K. Removing those parts is a bijection onto
  partitions of N − sK into parts ≤ k. □

Hence a_{k+1} = f_{k+1}(n−k−1) = f_{k+1}(m−1) = Σ_{s≥0} f_k(m−1−sK). This is a finite sum, because terms with m−1−sK < 0 vanish.

**1.2 The random variable X.** Fix x ∈ (0,1). Let Y₁, …, Y_k be independent, with P(Y_i = y) = (1−x^i)x^{iy} for y ≥ 0. These are
geometric laws, and each sums to 1. Put X = Σ_i iY_i.
- *Claim C:* P(X = N) = f_k(N)x^N Z with Z := Π_{i≤k}(1−x^i).
  *Proof.* P(X = N) = Σ over (y_1, …, y_k) with Σ i y_i = N of Π(1−x^i)x^{iy_i} = Z·x^N·#{(y_i) : Σ i y_i = N}. The count is f_k(N),
  where y_i is the multiplicity of the part i. □
- *Claim D:* E X = Σ_i i·x^i/(1−x^i) = Σ_i i/(e^{iθ}−1) for x = e^{−θ}.
  *Proof.* E Y_i = x^i/(1−x^i) for a geometric law with ratio x^i. □

So "θ is the saddle" (Σ_i i/(e^{iθ}−1) = m) is exactly E X = m. The left side is continuous and strictly decreasing in θ, from +∞ to 0,
so the saddle exists and is unique.

**1.3 Proof of (1).** f_k(m) ≥ 1 (the all-ones partition), so P(X = m) > 0. By Claim C,
f_k(N) = P(X = N)x^{−N}/Z. Then

  R_k = a_{k+1}/a_k = Σ_s f_k(m−1−sK)/f_k(m) = Σ_s x^{m−(m−1−sK)} P(X = m−1−sK)/P(X = m) = Σ_s x^{1+sK} P(X = m−1−sK)/P(X = m). □

**1.4 Fourier form.** Put χ(φ) := E e^{iφ(X−m)} = Σ_N P(X = N)e^{iφ(N−m)}. The series converges absolutely and uniformly, since
Σ_N P(X = N) = 1. Integrating term by term (allowed by uniform convergence):

  (1/2π)∫_{−π}^{π} χ(φ)e^{ijφ} dφ = Σ_N P(X = N)·(1/2π)∫e^{iφ(N−m+j)} dφ = P(X = m−j).

Insert this into (1). Since Σ_s x^{1+sK}|χ(φ)| ≤ Σ_s x^{1+sK} < ∞ uniformly in φ, sum and integral can be exchanged:

  R_k = ∫χ(φ)G(φ)dφ / ∫χ(φ)dφ,  G(φ) := Σ_s x^{1+sK}e^{i(1+sK)φ} = xe^{iφ}/(1 − qe^{iKφ}).

The last step is the geometric series with ratio |qe^{iKφ}| = q < 1. The denominator is D := ∫χ = 2πP(X = m) > 0.

**1.5 Proof of (2).** Put g₀ = G(0) = x/(1−q) and H(φ) := iφ − log(1 − qe^{iKφ}) + log(1−q). The logarithm is defined because
Re(1 − qe^{iKφ}) ≥ 1−q > 0, and exp(H) = e^{iφ}(1−q)/(1−qe^{iKφ}) = G/g₀. Therefore

  R_k/g₀ − 1 = (∫χG/g₀ − ∫χ)/∫χ = ∫χ(e^H − 1)/D = Ñ/D. □

**1.6 Factorisation of χ.** By independence and Claim D,

  χ(φ) = e^{−imφ}Π_i E e^{iφiY_i} = e^{−imφ}Π_i (1−x^i)/(1 − x^i e^{iiφ}).

Define ψ(φ) := −imφ + Σ_i[log(1−x^i) − log(1−x^i e^{iiφ})]. Each log is principal, with argument of real part ≥ 1−x^i > 0.
Then e^ψ = χ, and ψ is C^∞ (real-analytic).

**1.7 Machine sanity (ws_close_rep.out).**
- R1: Claims A and B exactly, for n ≤ 300.
- R2: (1) checked in exact rational arithmetic for all 4 ≤ n ≤ 60, all 2 ≤ k ≤ n−2 and x ∈ {1/3, 7/10, 19/20}, against p(n,k+1)/p(n,k). That is 4959 triples, all identical.
  Identity (1) holds for every x, because the x-dependence cancels as in 1.3. The x-free form is what makes it exact.
- R3: the law of X has mass 1 and mean m at the saddle (mpmath, 30 digits).
- R4: (2) by quadrature at (n,k) = (40,6), (60,9), (90,12), (120,7), agreeing to 10⁻²⁰. This includes D = 2πP(X = m).

---

## 2. Uniform derivative bounds (LemmaA.md §2) — trust item 2

**Lemma 2.1.** Let κ_r := Σ_{i≤k} i^r Li_{1−r}(x^i) for r ≥ 2. Then ψ(0) = ψ′(0) = 0, ψ^{(r)}(0) = i^rκ_r, and |ψ^{(r)}(φ)| ≤ κ_r for every φ ∈ ℝ and r ≥ 2.

*Proof.*
1. For 0 < y < 1 and a > 0, set ℓ_{y,a}(φ) := −log(1 − ye^{iaφ}). Because |ye^{iaφ}| = y < 1, the principal branch satisfies
   ℓ_{y,a}(φ) = Σ_{l≥1} y^l e^{ilaφ}/l, with absolute and uniform convergence.
2. The r-times differentiated series Σ_l (ila)^r y^{l}e^{ilaφ}/l = (ia)^r Σ_l l^{r−1}y^l e^{ilaφ} is dominated by a^r Σ_l l^{r−1}y^l < ∞
   (ratio test), uniformly in φ. Hence term-by-term differentiation is valid (Weierstrass M-test, applied inductively), and

   ℓ_{y,a}^{(r)}(φ) = (ia)^r Σ_l l^{r−1}y^l e^{ilaφ},  |ℓ_{y,a}^{(r)}(φ)| ≤ a^r Σ_l l^{r−1}y^l = a^r Li_{1−r}(y).

   For r ≥ 1, Li_{1−r}(y) := Σ_l l^{r−1}y^l.
3. ψ = −imφ + Σ_i[log(1−x^i) + ℓ_{x^i,i}]. For r ≥ 2 the linear and constant parts vanish after differentiation, so
   |ψ^{(r)}(φ)| ≤ Σ_i i^r Li_{1−r}(x^i) = κ_r. At φ = 0, ψ^{(r)}(0) = Σ_i (ii)^r Li_{1−r}(x^i) = i^rκ_r.
4. ψ(0) = 0 is immediate. ψ′(0) = −im + iΣ_i i·Li₀(x^i) = i(−m + Σ_i i x^i/(1−x^i)) = 0 by the saddle equation (Claim D). □

**Lemma 2.2.** Put c₁ := 1 + Kq/(1−q) and c_r := K^r Li_{1−r}(q) for r ≥ 2. Then H(0) = 0, H′(0) = ic₁, H″(0) = −c₂, H‴(0) = −ic₃,
H⁗(0) = c₄ and H^{(5)}(0) = ic₅. Also |H′| ≤ c₁ and |H^{(r)}| ≤ c_r (r ≥ 2) on ℝ, and Re H ≤ 0.

*Proof.*
1. H = iφ + ℓ_{q,K}(φ) + log(1−q). By step 2 of Lemma 2.1, H′ = i + iKΣ_l q^l e^{ilKφ}, so |H′| ≤ 1 + KΣ_l q^l = 1 + Kq/(1−q) = c₁.
2. For r ≥ 2, H^{(r)} = ℓ_{q,K}^{(r)}, so |H^{(r)}| ≤ K^r Li_{1−r}(q) = c_r.
3. At 0, H^{(r)}(0) = (iK)^r Li_{1−r}(q) = i^r c_r for r ≥ 2: i² = −1, i³ = −i, i⁴ = 1, i⁵ = i. Also H′(0) = i(1 + Kq/(1−q)) = ic₁.
4. Re H = log|G/g₀|. Since |1 − qe^{iKφ}|² − (1−q)² = 2q(1 − cos Kφ) ≥ 0 (sympy, ws_mech_lemA B3), |G| ≤ g₀, so Re H ≤ 0. □

**Corollary 2.3 (used in LemmaCD C0).** Let ω := ψ + Bφ²/2 with B = κ₂.
- ω(0) = ω′(0) = 0, ω″(0) = ψ″(0) + B = −κ₂ + B = 0, and |ω‴| = |ψ‴| ≤ κ₃.
- Re ψ is even and Im ψ is odd. This holds because χ(−φ) = conj χ(φ): X is real-valued, and ψ(−φ) = conj ψ(φ) termwise from the series.
- Im ψ(0) = (Im ψ)′(0) = (Im ψ)″(0) = 0, because ψ′(0) = 0 and ψ″(0) = −B is real. Also |(Im ψ)‴| ≤ |ψ‴| ≤ κ₃.
- (Re ψ)″(0) = −B, (Re ψ)‴(0) = Re(i³κ₃) = 0, and |(Re ψ)⁗| ≤ κ₄. □

**Machine sanity (ws_close_rep.out R5).** At 4 instances × 6 random φ ∈ (−π, π), mpmath numerical derivatives satisfy |ψ^{(r)}| ≤ κ_r (r = 2..6),
|H′| ≤ c₁, |H^{(r)}| ≤ c_r (r = 2..5) and Re H ≤ 0. Also ψ^{(r)}(0) = i^rκ_r, ψ′(0) = 0 and e^ψ = χ. The identities at 0 are also exact sympy
(ws_mech_lemA B1/B2).

---

## 3. The [ELEM] one-liners of ws_mech_lemA.py — trust item 3

Here t = |φ|. All functions are evaluated at real φ.

**3.1 D5: Re ω ≤ κ₄t⁴/24 + κ₆t⁶/720.**
1. Apply Taylor with N = 5 to ψ: ψ(φ) = Σ_{r=2}^{5} i^rκ_rφ^r/r! + R₆, with |R₆| ≤ κ₆|φ|⁶/720 (Lemma 2.1, sup|ψ^{(6)}| ≤ κ₆).
2. Take real parts. Re(i²κ₂φ²/2) = −Bφ²/2, Re(i⁴κ₄φ⁴/24) = κ₄φ⁴/24, and the r = 3 and r = 5 terms are purely imaginary.
3. Hence Re ω = Re ψ + Bφ²/2 = κ₄φ⁴/24 + Re R₆ ≤ κ₄t⁴/24 + κ₆t⁶/720. □

(Numeric sanity: ws_close_elem X3.)

**3.2 E2: the three pieces of T_H.**
- (i) |H(φ)| ≤ c₁t. This follows from H(φ) = ∫₀^φ H′ and |H′| ≤ c₁.
- (ii) For Re z ≤ 0: |e^z − 1 − z − z²/2| ≤ |z|³/6.
  The identity e^z − 1 − z − z²/2 = z³∫₀¹(1−s)²/2·e^{sz} ds holds; it is sympy-checked for z > 0 in ws_close_elem X1, and extends to all
  z ∈ ℂ because both sides are entire. For Re z ≤ 0, |e^{sz}| = e^{s Re z} ≤ 1, and ∫₀¹(1−s)²/2 ds = 1/6.
  Apply this with z = H, which has Re H ≤ 0, and use (i): ≤ c₁³t³/6.
- (iii) |H − ΣH_r| ≤ c₅t⁵/120. Here ΣH_r = H₁ + … + H₄ is the degree-4 Taylor polynomial of H at 0 (Lemma 2.2). Apply Taylor with N = 4
  and sup|H^{(5)}| ≤ c₅.
- (iv) |H² − H₁²|/2 ≤ c₁c₂t³/2. Write H² − H₁² = (H − H₁)(H + H₁). Then |H − H₁| ≤ c₂t²/2 (Taylor with N = 1, |H″| ≤ c₂), and
  |H + H₁| ≤ |H| + c₁t ≤ 2c₁t by (i).

With the exact decomposition e^H − P_H = (e^H−1−H−H²/2) + (H−ΣH_r) + (H²−H₁²)/2 (sympy, ws_mech_lemA E1), the sum of the three bounds is
T_H = c₅t⁵/120 + (c₁c₂/2 + c₁³/6)t³ (E3). □

**3.3 F4: L − L₂ = log(R/(g₀M₂)).** L₂ = −θ − log(1−q) + log M₂ = log(e^{−θ}/(1−q)) + log M₂ = log g₀ + log M₂, because g₀ = x/(1−q) with
x = e^{−θ}. Subtract this from L = log R. □

**3.4 F5: the minor-arc bounds.** The minor arcs φ₀ ≤ |φ| ≤ π have total length 2(π − φ₀) ≤ 2π. On them |e^H − 1| ≤ |e^H| + 1 ≤ 2, using Re H ≤ 0.
Hence |∫_{minor} χ(e^H − 1)| ≤ 2·2π·sup|χ| and |∫_{minor} χ| ≤ 2π·sup|χ|. □

**3.5 G3: t²S(t) ≥ F₁(Kt) − t.**
1. t²S(t) = t²Σ_{i≤k} i/(e^{it}−1) = tΣ_{i≤k} f₁(it), since f₁(u) = u/(e^u−1).
2. f₁ is decreasing (E1), so tf₁(it) ≥ ∫_{it}^{(i+1)t} f₁.
3. Summing over i ≤ k gives t²S ≥ ∫_t^{Kt} f₁ = F₁(Kt) − F₁(t).
4. F₁(t) = ∫₀^t f₁ ≤ t, because 0 < f₁ ≤ 1 (f₁(0⁺) = 1 and f₁ is decreasing).

The companion upper bound t²S(t) ≤ ζ(2) follows the same way from tf₁(it) ≤ ∫_{(i−1)t}^{it} f₁ and ∫₀^∞ f₁ = ζ(2). □

(Numeric sanity: X3.)

---

## 4. Lemma P1.1 and the closed form of F_r — trust item 4

**Lemma P1.1.** Let f: [0, ∞) → ℝ have total variation TV(f) < ∞ (here f = f_r, which is continuous). Then for θ > 0 and k ≥ 1,

  |θΣ_{i≤k} f(iθ) − ∫₀^{kθ} f| ≤ θ·Var_{[0,kθ]} f ≤ θ·TV(f).

*Proof.*
1. For each i, θf(iθ) − ∫_{(i−1)θ}^{iθ} f = ∫_{(i−1)θ}^{iθ}(f(iθ) − f(u))du.
2. For u in that interval, |f(iθ) − f(u)| ≤ Var_{[(i−1)θ, iθ]} f, by the definition of variation with the two-point partition {u, iθ}.
3. So the i-th difference has modulus ≤ θ·Var_{[(i−1)θ,iθ]} f.
4. Sum over i = 1..k. The intervals [(i−1)θ, iθ] are non-overlapping, and variation is additive over adjacent intervals,
   so Σ_i Var_{[(i−1)θ,iθ]} f = Var_{[0,kθ]} f ≤ TV(f). □

The constants TV(f_r) are certified by ws_mech_tv.py (MECH §3a, referee-audited as sound).

**Corollary P1.3.**
- s_r = θ^{r+1}κ_r = θΣ_{i≤k}(iθ)^r Li_{1−r}(e^{−iθ}) = θΣ f_r(iθ).
- So s_r ∈ F_r(kθ) ± θ·TV_r, with kθ = ν − θ.
- Since f_r ≥ 0, F_r is nondecreasing, and F_r(∞) = r!ζ(2) (Lemma 4.1 below). □

**Lemma 4.1 (closed form).** For r ≥ 1 and v > 0, F_r(v) = r!(ζ(2) − Σ_{j=0}^{r} v^j/j!·Li_{2−j}(e^{−v})).
Here Li₂ is the dilogarithm, Li₁(y) = −log(1−y), and Li_{−s}(y) = yP_s(y)/(1−y)^{s+1} (Eulerian).

*Proof.* Call the right side G_r(v).
1. **Derivative.** d/dv Li_s(e^{−v}) = −Li_{s−1}(e^{−v}). This follows from y d/dy Li_s(y) = Li_{s−1}(y), termwise from Li_s(y) = Σ y^l/l^s.
   The Leibniz rule then gives d/dv[v^j/j!·Li_{2−j}] = v^{j−1}/(j−1)!·Li_{2−j} − v^j/j!·Li_{1−j}, where the first term is absent for j = 0.
   The sum over j = 0..r telescopes to −v^r/r!·Li_{1−r}. Hence G_r′(v) = v^r Li_{1−r}(e^{−v}) = f_r(v).
   (sympy: ws_close_elem F1 for r = 1..8, with Li_s as formal symbols and the chain rule above; F2 checks the chain rule for sympy's polylog.)
2. **Value at 0⁺.**
   - Li₂(e^{−v}) → Li₂(1) = ζ(2).
   - For j = 1, v·Li₁(e^{−v}) = −v log(1−e^{−v}) → 0.
   - For j ≥ 2, Li_{2−j}(e^{−v}) = O(v^{−(j−1)}), so v^j Li_{2−j}(e^{−v}) = O(v) → 0.

   These limits are sympy-computed for j ≤ 8 in F2. Hence G_r(0⁺) = r!(ζ(2) − ζ(2)) = 0 = F_r(0).
3. **Conclusion.** f_r is continuous on [0, v] (it extends continuously to u = 0 with f_r(0) = (r−1)!), so F_r is C¹ with F_r′ = f_r.
   Two C¹ functions on (0, v] with equal derivatives and equal limits at 0⁺ coincide.
4. **F_r(∞).** As v → ∞ every Li_{2−j}(e^{−v}) is O(e^{−v}), so v^j Li_{2−j}(e^{−v}) → 0 and F_r(∞) = r!ζ(2). □

**Rigorous cross-check (ws_close_elem F3).** Arb's certified quadrature (acb.integral) of f_r over [0.25,1], [0.25,3], [1,9] and [3,30]
overlaps F_fun(r,b) − F_fun(r,a) for r = 1..6. This shows that the project's code F_fun implements Lemma 4.1. P1 checks P1.1 at three (θ, k)
instances with the certified TV.

---

## 5. Lemma B: B.1–B.4, the T₂ bound, the hull mean-value step, and the Region II ingredients — trust item 5

**Notation.** S_k(t) := Σ_{i≤k} i/(e^{it}−1), b(k,t) := tΣ_{i≤k} f₂(it), and B_k(t) := −S_k′(t).
- S_k is strictly decreasing, because each term is.
- S_k(θ) = m and S_{k+1}(θ′) = m−1, where θ′ is the saddle at k+1. Then m′ = n−k−1 = m−1.
- B_k(t) = Σ i²e^{it}/(e^{it}−1)² = t^{−2}Σ f₂(it) = b(k,t)/t³, with f₂(u) = u²e^u/(e^u−1)² = (u/2)²/sinh²(u/2).
  (sympy, ws_close_elem B1, k ≤ 6; for general k it is the same identity termwise.)
- S_{k+1}(t) − S_k(t) = K/(e^{Kt}−1). At t = θ this equals Kq/(1−q) = c₁ − 1 (B1).

Fix τ = 0.01 and b_min := F₂(ν−θ) − θ(1+τ)TV₂.

**Claim 5.0.** b(k,t) ≥ b_min for every t ∈ [θ, (1+τ)θ].
*Proof.* P1.3 gives b(k,t) ≥ F₂(kt) − t·TV₂. Since F₂ is nondecreasing and kt ≥ kθ = ν−θ, F₂(kt) ≥ F₂(ν−θ). Also t·TV₂ ≤ (1+τ)θ·TV₂. □

**B.1.** If τ·b_min > (1+τ)³c₁θ², then θ < θ′ ≤ (1+τ)θ.
*Proof.*
1. S_{k+1}(θ) = m + c₁ − 1 > m − 1 = S_{k+1}(θ′). S_{k+1} is strictly decreasing, so θ′ > θ.
2. On [θ, (1+τ)θ], −S_{k+1}′(t) = B_{k+1}(t) ≥ B_k(t), since the extra term is positive.
   And B_k(t) = b(k,t)/t³ ≥ b_min/((1+τ)θ)³ (Claim 5.0).
3. Integrating, S_{k+1}((1+τ)θ) ≤ S_{k+1}(θ) − τθ·b_min/((1+τ)θ)³ = m + c₁ − 1 − τb_min/((1+τ)³θ²).
4. By the hypothesis, τb_min/((1+τ)³θ²) > c₁, so S_{k+1}((1+τ)θ) < m − 1 = S_{k+1}(θ′). Monotonicity then gives θ′ < (1+τ)θ. □

**B.2.** Under B.1, dθ := θ′ − θ ≤ (1+τ)³c₁θ³/b_min.
*Proof.* By the mean value theorem, c₁ = S_{k+1}(θ) − S_{k+1}(θ′) = B_{k+1}(ξ)·dθ for some ξ ∈ (θ, θ′) ⊂ [θ, (1+τ)θ].
By step 2 of B.1, B_{k+1}(ξ) ≥ b_min/((1+τ)θ)³. Divide. □

**B.3.** dν := ν_{k+1} − ν_k = (K+1)θ′ − Kθ = θ′ + K·dθ. This is > 0, and ≤ (1+τ)θ + (ν/θ)dθ ≤ (1+τ)θ + (ν/θ + 1)dθ.
The algebra is sympy-checked (B5) and uses K = ν/θ. □

**B.4.** For r = 2..5, |s_r(k+1,θ′) − s_r(k,θ)| ≤ dθ·((1+r)s̄_r + s̄_{r+1})/θ + (1+τ)θ·sup_{u∈[ν, ν+dν]} f_r(u),
where s̄_r := r!ζ(2) + (1+τ)θTV_r.
*Proof.*
1. s_r(k+1,θ′) = s_r(k,θ′) + θ′f_r(Kθ′) (sympy B4).
2. ∂_t s_r(k,t) = ((1+r)s_r(k,t) − s_{r+1}(k,t))/t. This uses u f_r′ = r f_r − f_{r+1} (sympy B4 and E14).
3. For t ∈ [θ, θ′]: t ≥ θ, and 0 ≤ s_r(k,t) ≤ F_r(kt) + tTV_r ≤ r!ζ(2) + (1+τ)θTV_r = s̄_r (P1.3, Lemma 4.1).
   So |∂_t s_r| ≤ ((1+r)s̄_r + s̄_{r+1})/θ, and |s_r(k,θ′) − s_r(k,θ)| ≤ dθ·that.
4. The remaining term: θ′ ≤ (1+τ)θ, and Kθ′ ∈ [Kθ, (K+1)θ′] = [ν, ν_{k+1}] ⊂ [ν, ν+dν]. □

**T₂ bound.** L₂,k − L₂,k+1 = (θ′ − θ) + log((1−q′)/(1−q)) − (log M₂′ − log M₂).
1. The first term is > 0 by B.1.
2. q′ = e^{−(K+1)θ′} < e^{−(K+1)θ} = qe^{−θ}.
3. So (1−q′)/(1−q) > (1 − qe^{−θ})/(1−q) = 1 + q(1−e^{−θ})/(1−q). The identity is sympy B5. Take logs.

Hence step ≥ T₂ − |Δ log M₂|. □

**Hull mean-value step (Region I).**
1. M₂ is an explicit rational-exponential function of v = (θ, ν, b, s₃, s₄, s₅); see ws_b_core.M2_formula, which equals LemmaA §3 by
   ws_mech_lemA C2/C3. Let P₀ = v(k) and P₁ = v(k+1).
2. ws_b_core.hull_point builds a product box H ⊂ ℝ⁶ with the following coordinates:
   - θ ∈ [θ, (1+τ)θ] ∋ θ, θ′ (B.1);
   - ν ∈ [ν, ν + dν] ∋ ν, ν_{k+1} (B.3);
   - s_r ∈ [F_r(ν−θ) − (1+τ)θTV_r, F_r(ν+dν) + (1+τ)θTV_r].
3. The s_r interval contains s_r(k,θ) ∈ F_r(ν−θ) ± θTV_r, and s_r(k+1,θ′) ∈ F_r(Kθ′) ± θ′TV_r with Kθ′ ∈ [ν, ν+dν] and F_r nondecreasing.
   So P₀, P₁ ∈ H.
4. A box is convex, so the segment [P₀, P₁] ⊂ H. The script asserts M₂ > 0 on H, so log M₂ is C¹ on a neighbourhood of the segment.
5. By the mean value theorem applied to φ(λ) = log M₂(P₀ + λ(P₁−P₀)), there is a ξ on the segment with
   |log M₂(P₁) − log M₂(P₀)| = |Σ_j ∂_j log M₂(ξ)(P₁−P₀)_j| ≤ Σ_j sup_H|∂_jM₂/M₂|·|Δv_j|.
6. The sup is enclosed by forward-mode AD over Arb balls. That is the chain rule evaluated in interval arithmetic, so it encloses every
   value of the derivative on H. |Δv_j| is bounded by B.2–B.4. □

**Validation at real saddles (ws_close_elem B6).** At (n,k) = (10⁵, 870), (10⁵, 1358), (10⁵, 1850), (10⁶, 3700), (10⁶, 5191) and (10⁶, 6700),
θ_k and θ_{k+1} are enclosed by Arb bisection. The script then checks:
- B.1, B.2, B.3 and B.4 at these points;
- that both endpoints lie in the hull box;
- that the true |Δ log M₂| is below the MVT bound.

**Region II ingredients (ws_b_asym.py), each with proof.**
- **(a)** sup θ′^e over both points is (1+τ)^eθ^e if e > 0 and θ^e if e < 0, because θ ≤ θ′ ≤ (1+τ)θ.
- **(b)** q′ ≤ q, by step 2 of the T₂ bound.
- **(c)** ν_{k+1} ≤ ρν with ρ = 1 + (1+τ)θ₁/ν_min + (1+τ)³(c₁θ²)_max(1 + θ₁/ν_min)/b_min.
  By B.3 and B.2, ν_{k+1}/ν ≤ 1 + (1+τ)θ/ν + (1/θ + 1/ν)(1+τ)³c₁θ³/b_min, which equals the expression with θ in place of θ₁. Now use
  θ ≤ θ₁ and ν ≥ ν_min. The bound (c₁θ²)_max follows from:
  - c₁θ² = θ² + θ²νe^{−ν′}/(1−q);
  - d/dθ[θ²ν] = θ(2ν−1) > 0 at fixed ν′ (sympy B5), so θ²ν ≤ θ₁²(ν′_hi + log(1/θ₁));
  - e^{−ν′} ≤ e^{−v₀};
  - q ≤ θ₁e^{−v₀}.
- **(d)** f_r(u) = u^r e^{−u}P_{r−1}(e^{−u})/(1−e^{−u})^r (Eulerian, sympy B5). For u ∈ [ν, ν_{k+1}], u ≤ ρν, e^{−u} ≤ q ≤ y := e^{−ν_min},
  and P_{r−1}(y)/(1−y)^r is increasing in y. So f_r(u) ≤ (ρν)^r q P_{r−1}(y)/(1−y)^r.
- **(e)** |g_r′(u)| ≤ (r g_r + g_{r+1})/u, since g_r = f_r and u f_r′ = r f_r − f_{r+1} with f ≥ 0. Then u ≥ ν.
- **(f)** c₁ = 1 + h(ν)/θ with h(u) = ue^{−u}/(1−e^{−u}).
  - h′ = (e^u − 1 − ue^u)/(e^u−1)². The numerator is ≤ 0: it is 0 at u = 0 and has derivative −ue^u.
  - So |h′| = (ue^u − e^u + 1)/(e^u−1)² ≤ (u+1)e^u/(e^u−1)² = (u+1)e^{−u}/(1−e^{−u})², since (u+1)e^u − (ue^u − e^u + 1) = 2e^u − 1 > 0.
    All of this is sympy B5.
  - Δc₁ = (h(ν_{k+1}) − h(ν))/θ′ + h(ν)(1/θ′ − 1/θ). The first part is ≤ dν·sup|h′|/θ. The second has modulus ≤ h(ν)dθ/θ², because
    θ′ ≥ θ and h(ν) ≤ ρνq/(1−q₁). This is exactly the script's dc1.
- **(g)** Products and quotients are handled with
  - |Δ(xy)| ≤ |Δx|·sup|y| + sup|x|·|Δy|;
  - |Δ(1/N)| ≤ |ΔN|/N_min²;
  - |Δ b^{−r/2}| ≤ (r/2)b_min^{−r/2−1}|Δb|;
  - |Δθ^e| ≤ e·sup θ^{e−1}·dθ for the exponents used, all e > 0;
  - |Δ log M₂| ≤ |ΔM₂|/min M₂ (MVT for log).
  Each is the mean value theorem, or the telescoping Δ(xy) = Δx·y′ + x·Δy.
- **(h)** T₂/(θq) ≥ (1−θ/2)(1 − θq/(2(1−q))). Write x = q(1−e^{−θ})/(1−q). Then:
  - x ≥ q(θ − θ²/2), since 1 − e^{−θ} ≥ θ − θ²/2;
  - x ≤ θq/(1−q);
  - log(1+x) ≥ x − x²/2 for x ≥ 0.

  Then use θ ≤ θ₁ and q ≤ q₁.
- **(i)** Every majorant is a nonnegative combination of t^p q^{b_q} ν^j. After dividing by t²q, each monomial is t^Pν^je^{−β_qν′}.
  - At fixed ν′, d/dt[t^Pν^j] = t^{P−1}ν^{j−1}(Pν − 2j) (sympy E7). This is ≥ 0 when ν ≥ ν_min > 2j/P, which the script asserts.
    For j < 0 it holds trivially.
  - Each monomial is therefore maximised at t = t₁ = √θ₁, with ν bounded by v₁ + log(1/θ₁) when j ≥ 0 and by ν_min when j < 0.
  - The factor e^{−β_qν′} is bounded at the appropriate end of the ν′ cell.

I re-read ws_b_asym.py line by line against (a)–(i), and every term matches.

---

## 6. Lemma C0, C1 and the CL-a minor-arc terms — trust item 6

**Lemma C0.** Let φ₀ = c₀θ and suppose Re ψ(φ) ≤ −(1−δ)Bφ²/2 on |φ| ≤ φ₀ with 0 ≤ δ < 1. Define T₁ = c₁κ₃/(2B²), T₂ = (c₂+c₁²)/(2B),
A₃ = κ₃B^{−3/2}, w₀ = φ₀√B and M = √(2πB)·sup_{φ₀≤|φ|≤π}|χ|. If the denominator below is positive, then

  |R_k/g₀ − 1| ≤ ρ := (T₁/(1−δ)^{5/2} + T₂/(1−δ)^{3/2} + 2M) / (1 − erfc(w₀/√2) − √(2/π)A₃/(3(1−δ)²) − M).

*Proof.* By (2), R_k/g₀ − 1 = Ñ/D. All integrals below are taken in units of √(2π/B).

**(i) Central numerator.**
1. Taylor with N = 1 for e^H: e^H − 1 = ic₁φ + E(φ), since (e^H)′(0) = H′(0) = ic₁.
2. (e^H)″ = (H″ + H′²)e^H, so |(e^H)″| ≤ (c₂ + c₁²)·|e^H| ≤ c₂ + c₁² (Re H ≤ 0). Hence |E| ≤ (c₂+c₁²)φ²/2.
3. Since χ(−φ) = conj χ(φ), ∫_{−φ₀}^{φ₀}χφ dφ = ∫₀^{φ₀}φ(χ(φ) − conj χ(φ))dφ = 2i∫₀^{φ₀}φ·Im χ.
4. Im χ = e^{Re ψ}sin(Im ψ), so |Im χ| ≤ e^{Re ψ}|Im ψ| ≤ e^{−(1−δ)Bφ²/2}κ₃|φ|³/6. The last step is Taylor with N = 2 for Im ψ
   (Corollary 2.3).
5. Hence |c₁∫χφ| ≤ 2c₁∫₀^{φ₀}e^{−(1−δ)Bφ²/2}κ₃φ⁴/6 ≤ (c₁κ₃/6)∫_ℝ φ⁴e^{−(1−δ)Bφ²/2} = √(2π/B)·T₁/(1−δ)^{5/2}.
   Likewise ∫|χ||E| ≤ ((c₂+c₁²)/2)∫_ℝ φ²e^{−(1−δ)Bφ²/2} = √(2π/B)·T₂/(1−δ)^{3/2}.
   Both Gaussian integrals are sympy C0a, with 1−δ a positive symbol.

**(ii) Minor numerator.** By §3.4, the minor arcs contribute ≤ 2·2π sup|χ| = √(2π/B)·2M (sympy C0c).

**(iii) Denominator.**
1. D is real and positive (§1.4), so D ≥ Re∫_{central}χ − |∫_{minor}χ|.
2. Write χ = e^{−Bφ²/2}e^ω. By Corollary 2.3 and Taylor with N = 2, |ω| ≤ κ₃|φ|³/6.
3. e^ω − 1 = ω∫₀¹e^{sω}ds (X1), so |e^ω − 1| ≤ |ω|·max(1, e^{Re ω}). Re ω = Re ψ + Bφ²/2 ≤ δBφ²/2 and δ ≥ 0, so the max is ≤ e^{δBφ²/2}.
4. Hence |χ − e^{−Bφ²/2}| ≤ e^{−(1−δ)Bφ²/2}κ₃|φ|³/6, and its integral over ℝ is √(2π/B)·√(2/π)A₃/(3(1−δ)²) (C0a).
5. ∫_{|φ|≤φ₀}e^{−Bφ²/2} = √(2π/B)(1 − erfc(w₀/√2)) (C0b). The minor part is ≤ 2π sup|χ| = √(2π/B)M.

**(iv) Conclusion.** |Ñ|/D ≤ ρ whenever the denominator lower bound is positive. □

**Admissible δ.**
- **δ_T.** Re ψ is even, with Re ψ(0) = 0, (Re ψ)″(0) = −B, odd derivatives zero at 0, and |(Re ψ)⁗| ≤ κ₄ (Corollary 2.3).
  Taylor with N = 3 gives Re ψ ≤ −Bφ²/2 + κ₄φ⁴/24 ≤ −Bφ²/2 + κ₄φ₀²φ²/24 = −(1 − δ_T)Bφ²/2, with δ_T = κ₄φ₀²/(12B).
- **δ_C** (when c₀v < 2).
  1. |χ|² = Π_i (1−y_i)²/|1−y_ie^{iiφ}|². Each factor is (1 + sin²(iφ/2)/sinh²(iθ/2))^{−1} (sympy C0d), so
     −2Re ψ = Σ_i log(1 + sin²(iφ/2)/sinh²(iθ/2)).
  2. For |φ| ≤ c₀θ and i ≤ k, the angle x = i|φ|/2 ≤ c₀v/2 < 1.
  3. sin x ≥ x − x³/6 ≥ 0 (E10), and (x − x³/6)² = x²(1 − x²/3) + x⁶/36 (sympy C0d). So sin²x ≥ x²(1 − x²/3) ≥ λx²
     with λ := 1 − (c₀v)²/12 ∈ (0,1].
  4. With a_i := (iφ/2)²/sinh²(iθ/2) ≤ (φ/θ)² ≤ c₀² (sinh y ≥ y): log(1 + sin²(iφ/2)/sinh²(iθ/2)) ≥ log(1 + λa_i) ≥ λlog(1 + a_i).
     The last step is E4 (concavity).
  5. log(1 + a) ≥ a·log(1+c₀²)/c₀² for 0 ≤ a ≤ c₀², by E4 again (chord).
  6. Σa_i = φ²Σ_i i²e^{iθ}/(e^{iθ}−1)² = Bφ² (sympy C0d).
  7. So −2Re ψ ≥ λ(log(1+c₀²)/c₀²)Bφ², which is the claim with δ_C = 1 − λlog(1+c₀²)/c₀². This lies in [0,1), since log(1+c²) < c².

  ws_cd_common.X_cell uses δ_C only when c₀v₂ < 2, and with v₂ ≥ v; the formula is monotone in v.

**Consequences.**
- log R = −θ − log(1−q) + log(R/g₀).
- −log(1−q) ≤ q/(1−q) and −log(1−q) ≥ q. Both follow from the series −log(1−q) = Σq^l/l.
- log(1+ρ) ≤ ρ.
- For ρ < 1, log(1−ρ) ≥ −ρ/(1−ρ). This is the series again: −log(1−ρ) = Σρ^l/l ≤ Σρ^l.
- 1 − q = 1 − e^{−ν} ≤ ν.

These give the three displayed inequalities of LemmaCD §1. □

**Lemma C1 (θ-uniform cell bound).** The ingredient bounds are as follows. Each is proved, and the scaled identities are sympy C1a.
- s_r ≤ F_r(v₂) + θTV_r: P1.3, with F_r nondecreasing.
- s_r ≤ kθ·sup f_r = v·sup f_r, since each term θf_r(iθ) ≤ θ sup f_r.
  - The sup values come from ws_cd_common.sup_f. On [a,b] it uses f_r = g^r h_r with g(u) = u/(1−e^{−u}) increasing, so g ≤ g(b).
  - h_r(u) = e^{−u}P_{r−1}(e^{−u}) is decreasing, because P has positive coefficients, so h_r ≤ h_r(a).
  - On [60, ∞), f_r is decreasing (E6).
- b ≥ F₂(v₁) − θ (§3.5 argument with f₂ in place of f₁, f₂ decreasing (E2), 0 < f₂ ≤ 1). Also b ≥ kθf₂(kθ) ≥ v₁f₂(v₂), and b ≤ F₂(kθ) ≤ F₂(v₂)
  (right-endpoint Riemann sum of a decreasing function).
- c₁θ = θ + f₁(ν) ≤ θ + f₁(v₁) and c₂θ² = f₂(ν) ≤ f₂(v₁), since ν = v + θ ≥ v₁ and f₁, f₂ are decreasing (sympy C1a).
- M ≤ √(2πb/θ³)e^{−Ψ/(2θ)}, because sup|χ| = e^{−inf S/2} ≤ e^{−Ψ/(2θ)}.

Monotonicity in θ:
1. ρ is increasing in T₁, T₂, M, δ, A₃ and erfc(w₀/√2), and the numerator is increasing and the denominator decreasing in each (sympy C1b).
2. After substituting the bounds, T₁/θ, T₂/θ, A₃ and δ are nondecreasing in θ, and w₀ is nonincreasing in θ. Each is an explicit expression
   in θ through θ, F₂(v₁)−θ and θTV_r, with signed derivative (ws_mech_mono).
3. M/θ = √(2πb)θ^{−5/2}e^{−Ψ/(2θ)} has d/dθ log = (Ψ − 5θ)/(2θ²) (C1a). This is > 0 when Ψ > 5θ̄ ≥ 5θ, which X_cell asserts.

Therefore ρ/θ is bounded by its value with every ingredient taken at θ̄, and that value is X_cell. □

**CL-a minor-arc terms (ws_mech_cla.py).** Here N := √(B/2π)∫_{φ₀≤|φ|≤π}|χ| replaces M. In C0, M entered only through
∫_{minor}|χ| ≤ 2π sup|χ|, and N bounds the same integrals directly. The ranges are φ ∈ [c₀θ, 10θ], [10θ, π/k] and [π/k, π].
They cover the minor arc, since 10θ = 10v/k ≤ 1/k < π/k.
- **N₁, c = φ/θ ∈ [c_j, c_{j+1}].**
  1. iφ/2 ≤ c_{j+1}kθ/2 ≤ c_{j+1}V₀/2 = T_j < π/2 (asserted), so sin(iφ/2) ≥ s_j·iφ/2 with s_j = sinT_j/T_j (E3).
  2. So the summand is ≥ log(1 + s_j²c²a′_i), where a′_i = f₂(iθ) ≤ 1.
  3. By E4 (chord on [0, s_j²c_{j+1}²]) this is ≥ c²a′_iλ_j with λ_j = log(1+s_j²c_{j+1}²)/c_{j+1}².
  4. Σc²a′_i = Bφ², so |χ| ≤ e^{−λ_jBφ²/2}.
  5. 2√(B/2π)∫_{c_jθ}^∞e^{−λ_jBφ²/2} = erfc(c_jθ√(λ_jB/2))/√λ_j (sympy CLa). Here θ√B = √(b/θ) ≥ √(kf₂(v)) ≥ √(K₀f₂(V₀)), and erfc is decreasing.
- **N₂, φ ∈ [10θ, π/k].**
  1. iφ/2 ≤ π/2, so sin(iφ/2) ≥ (2/π)(iφ/2).
  2. The summand is ≥ log(s²c²f₂(iθ)) ≥ log(s²c²f₂(V₀)), using f₂ decreasing and iθ ≤ v ≤ V₀. So |χ| ≤ (s²f₂(V₀))^{−k/2}c^{−k}.
  3. N₂ ≤ √(B/2π)·2θ∫_{10}^∞(…)dc = √(2/π)(θ√B)(s²f₂(V₀))^{−k/2}10^{1−k}/(k−1) (sympy CLa).
  4. θ√B = √(b/θ) ≤ √(v/θ) = √k, since b ≤ v.
- **N₃₄, φ ∈ [π/k, π], d = φ/2 ∈ [π/(2k), π/2].**
  1. The summand is ≥ sin²(id)w_i with w_i = log(1 + 1/sinh²(iθ/2)) (E4), and w_i ≥ w_k (decreasing).
  2. So S ≥ w_kΣ_{i≤k}sin²(id) = w_k(k/2 − sin(kd)cos((k+1)d)/(2 sin d)) ≥ w_k(k/2 − 1/(2 sin d)) (E9).
  3. sin d ≥ sin(π/(2k)) = (π/(2k))sinc(π/(2k)) ≥ (π/(2k))sinc(π/(2K₀)) (E3, k ≥ K₀). So the sum is ≥ k(½ − 1/(π sinc(π/2K₀))) = 2gk, with g > 0 asserted.
  4. w_k ≥ −2 log sinh(v/2) ≥ −2 log(vσ₀/2), since sinh(v/2)/(v/2) ≤ σ₀ for v ≤ V₀ (E2).
  5. So |χ| = e^{−S/2} ≤ (vσ₀/2)^{2gk}, and N₃₄ ≤ √(B/2π)·2π·(…) = √(2πB)(…) ≤ √(2π)k^{3/2}v^{−1}(vσ₀/2)^{2gk}. This uses B = b/θ³ ≤ v·k³/v³.
  6. The bound is increasing in v (exponent 2gK₀ − 1 ≥ 0 asserted, base factor σ₀/2 fixed) and decreasing in k (ratio < 1 asserted).
     So it is evaluated at (K₀, V₀).

Finally log R_k ≥ −θ − log ν + log(1−ρ) > 0 iff 1 − ρ > e^θν. This holds because ν = v(1 + 1/k) ≤ V₀(1 + 1/K₀) and θ ≤ V₀/K₀. □

---

## 7. Lemma D (D1–D3), including "Ψ(s) increasing" — trust item 7 (PRIORITY: now proved analytically and mechanically)

**D1.** C(n−1,k−1)/k! ≤ p(n,k) ≤ C(n−1+D,k−1)/k! with D = k(k−1)/2.
*Proof.*
- *Lower bound.* Sorting a composition of n into k positive parts gives a partition of n into k parts. This map is onto, and each fibre has
  at most k! elements (the orderings of the parts). So C(n−1,k−1) = #compositions ≤ k!·p(n,k).
- *Upper bound.* For λ₁ ≥ … ≥ λ_k ≥ 1 with Σλ_i = n, put μ_i = λ_i + (k−i). Then μ₁ > μ₂ > … > μ_k ≥ 1 and Σμ_i = n + D. The map λ ↦ μ is
  injective (λ_i = μ_i − (k−i)). A partition of n+D into k distinct parts corresponds to exactly k! compositions (all orderings distinct).
  So k!·p(n,k) ≤ k!·#{distinct-part partitions} ≤ #compositions of n+D into k parts = C(n+D−1,k−1). □

(Exact check, all 1 ≤ k ≤ n ≤ 250: ws_close_rep D1.)

**D2.** If Φ(k,n) := log(n−k) − log(k(k+1)) − k(k−1)²/(2(n−k+1)) > 0, then p(n,k+1) > p(n,k) (for 1 ≤ k ≤ n−1).
*Proof.*
1. By D1 it suffices that C(n−1,k)/(k+1)! > C(N−1,k−1)/k! with N = n+D, that is, C(n−1,k)/C(N−1,k−1) > k+1.
2. Write C(n−1,k)/C(N−1,k−1) = [C(n−1,k)/C(n−1,k−1)]·[C(n−1,k−1)/C(N−1,k−1)]. The first factor is (n−k)/k. The second is
   Π_{j=0}^{k−2}(n−1−j)/(n−1−j+D) = Π(1 + D/(n−1−j))^{−1}. Both are sympy D2a.
3. Since 1 + x ≤ e^x (D2b) and n−1−j ≥ n−k+1 for j ≤ k−2, the product is ≥ exp(−(k−1)D/(n−k+1)) = exp(−k(k−1)²/(2(n−k+1))).
4. Taking logs, the sufficient condition becomes log(n−k) − log k − k(k−1)²/(2(n−k+1)) > log(k+1), which is Φ > 0. □

**D3.** For n ≥ 10⁵ and every integer 1 ≤ k ≤ 1.7n^{1/3}, Φ(k,n) > 0.
*Proof.* Put γ = 1.7 and s = n^{1/3} ≥ s_A := 10^{5/3} = 46.4158….
1. **Φ is decreasing in real k ∈ [1, n−1].** dΦ/dk = −1/(n−k) − 1/k − 1/(k+1) − (k−1)[(3k−1)(n−k+1) + k(k−1)]/(2(n−k+1)²) (sympy D3a).
   Every term is negative or zero for 1 ≤ k < n. Since γs ≤ n−1 (1.7·46.5 ≪ 10⁵), Φ(k,n) ≥ Φ(γs, s³) for k ≤ γs.
2. **Φ(γs, s³) ≥ Ψ(s) := log((s²−γ)/(γ(γs+1))) − γ³/(2(1 − γ/s²)).**
   - The log terms agree exactly: (s³ − γs)/(γs(γs+1)) = (s²−γ)/(γ(γs+1)) (sympy D3b).
   - For the cubic term, k(k−1)²/(n−k+1) ≤ k³/(n−k), because k³(n−k+1) − k(k−1)²(n−k) = k[(n−k)(2k−1) + k²] ≥ 0 (sympy D3b).
     At k = γs, n = s³, the right side is γ³/(1 − γ/s²).
3. **Ψ is increasing on s > √γ.** Exactly (sympy D3c):

   **Ψ′(s) = (γs² + 2s + γ²)/((s² − γ)(γs + 1)) + γ⁴s/(s² − γ)².**

   For s > √γ = 1.30…, every factor is positive, so Ψ′ > 0. This replaces the hand argument of REFEREE A.3 with an exact closed form.
   As a redundant check, ws_close_D D3d also certifies Ψ′ > 0 by an Arb sweep on [s_A, 10⁴]; the closed form covers s > 10⁴.
4. **Value.** Ψ(s_A) ∈ [0.30456254 ± 4·10⁻⁹] > 0 (Arb, with s_A enclosed as the ball 10^{5/3}). Hence Ψ(s) ≥ Ψ(s_A) > 0 for all s ≥ s_A, that is,
   for all n ≥ 10⁵. □

**Corollary.** For n ≥ 10⁵ and 1 ≤ k ≤ k_D(n) = ⌊1.7n^{1/3}⌋, p(n,k+1) > p(n,k).
k_D(n) = max{k : 1000k³ ≤ 4913n} is nondecreasing, and k_D(10⁵) = 78 (exact integers, ws_close_D D3e).
Exact sanity: a_{k+1} > a_k for all k ≤ k_D at n = 10⁵ (k_D = 78) and at n = 2·10⁵ (k_D = 99), from big-integer columns.

---

## 8. Region E — trust item 8

**Lemma E.** If 2k ≥ n, then p(n,k) = p(n−k). For ⌈n/2⌉ ≤ k ≤ n−1, a_{k+1} ≤ a_k, with equality only at k = n−1.
*Proof.*
1. Subtract 1 from each of the k parts of a partition of n into exactly k parts, and delete the zeros. The result is a partition of n−k into at
   most k parts. The inverse map pads with zeros to length k and adds 1 to each entry, so this is a bijection.
2. Every partition of n−k has at most n−k ≤ k parts, so "at most k parts" is no restriction. Hence p(n,k) = p(n−k).
3. If k ≥ ⌈n/2⌉, then also k+1 ≥ n/2, so a_{k+1} = p(n−k−1) and a_k = p(n−k).
4. Appending a part 1 is an injection from partitions of j into partitions of j+1, so p(j) ≤ p(j+1). For j ≥ 1 it is not onto, because the
   partition (j+1) is missed, so p(j) < p(j+1).
5. With j = n−k−1: a_{k+1} = p(j) ≤ p(j+1) = a_k, with strict inequality iff j ≥ 1, that is, k ≤ n−2. □

(Exact check, all n ≤ 400: ws_close_rep E1, E2.)

---

## 9. Threshold audit: every hypothesis holds for every n ≥ 10⁵ (ws_close_thresh.py → .out)

### 9.1 Parameter constants for all n
Every constant below is certified on all β ≥ β_A by ws_mech_window.py. It uses 5486 Arb boxes on ℓ ∈ [ℓ_A, 60] plus one tail box for
ℓ ≥ 60 in the atoms x = 1/β, ℓx, ℓ²x and 1/ℓ. Those atoms are monotone for ℓ ≥ 60 by E6. No monotonicity in β is assumed.
- θ ≤ θ_A = 0.0040939
- ν′ ∈ [−2.5101980, 2.0839964]
- ν ≥ 3.1158184
- ν_CL ≤ −1.9536125
- ν_CR ≥ 1.0430522
- ν′_L ≤ −1.9066744 (printed form)
- G ≥ 1.0092

### 9.2 Lemma by lemma

| Lemma | Hypothesis | Holds for all n ≥ 10⁵ because |
|---|---|---|
| P2 (Lemma A minor arcs) | θ ≤ 0.0042, (k+1)θ ≥ 3.0, φ ≥ 0.4θ | θ ≤ θ_A, ν ≥ 3.1158; c_* = 0.15973 ≥ 0.15926 (value used); c_*/3.5 > 0.002 (Region II tail monotonicity) |
| Thm A-unif, Region I | θ ∈ [0.002, θ_A], ν′ ∈ [−2.52, 2.10], ν ≥ 3.1158 | 9.1; the box cover runs to 2.11 |
| Thm A-unif, Region II | θ ≤ 0.002, ν′ ∈ [−2.52, 2.10], ν ≥ v₀ + log 500 ≥ 3.69 | 9.1 |
| Lemma B, Regions I/II | the same ranges at k; B.1's τ-step asserted per box/cell | 9.1; quantities at k+1 come from B.1–B.4, not from Lemma W at k+1 |
| Corollary of A-unif | Σ_min/1.01 > 2Γ_max | 0.91918/1.01 = 0.91008 > 0.8564; Γ certified ≤ 0.42024 ≤ 0.4282 |
| (B3) | ν′_L ≤ −1.906674, ν′_R ≥ 1.0430522, \|log M₂\|/θ ≤ 0.14236 | 9.1; the certified sup is 0.14215 |
| CR | θ ≤ θ_CR := √(2ζ(2)/(n−1)) | m ≥ ⌊n/2⌋ ≥ (n−1)/2 and θ²m ≤ ζ(2) (§3.5). θ_CR is decreasing in n, ≤ 0.0057358 at n = 10⁵ |
| CR | v ≥ 6.198 | v = ν′ − log θ − θ ≥ 1.0430522 + log(1/θ_CR) − θ_CR ≥ 6.19835, and this increases with n |
| CR | Ψ > 5θ̄ | 0.33808 > 0.0287 |
| CL-a | k ≥ 47 and kθ ≤ 0.1 | independent of n |
| CL-b | k ≥ 79, θ ≤ θ_A, ν′ ≤ ν_CL | k > k_D ≥ 78; k ≤ μ−2β; 9.1 |
| D | 1 ≤ k ≤ 1.7n^{1/3}, n ≥ 10⁵ | §7 (Ψ increasing) |
| Representation (1)/(2) | 2 ≤ k ≤ n−2 | CL, W and CR have 79 ≤ k ≤ ⌈n/2⌉ ≤ n−2 |

### 9.3 Cover and assembly
The regions are:
- D = [1, k_D];
- CL = (k_D, μ−2β];
- W = (μ−2β, μ+2β);
- CR = [μ+2β, ⌈n/2⌉];
- E = [⌈n/2⌉, n−1].

**These cover [1, n−1].** This needs μ+2β = β(ℓ+2) ≤ n/2. The ratio β(ℓ+2)/(n/2) = (2/ζ(2))(ℓ+2)e^{−ℓ} is decreasing for ℓ > −1, and is < 0.04
at n = 10⁵ (Arb). If CL is empty (k_D ≥ μ−2β), nothing changes.

**Sign pattern.**
- R_k > 1 on D ∪ CL, R_k < 1 on CR, and R_k ≤ 1 on E.
- On W, every k except at most one has sign L_k = sign L₂,k (Corollary of A-unif). L₂ is strictly decreasing on W (B1).
- So the determined signs on W read +…+ −…−, with at most one sign change. The single undetermined k lies in the crossing pair, adjacent to the
  change, so whatever its sign, the pattern on W is still +…+ −…−.
- Concatenating, the whole row reads +…+ −…−, which is weak unimodality.

**Remark.** (B3) (the edge signs) is not needed for this conclusion. Any pattern +…+ −…− on W, including all + or all −, concatenates with the
+ of CL and the − of CR/E to a single change. (B3) only locates the mode inside W. It is a correct, certified statement, but it is not
load-bearing.

### 9.4 Finite part
n ≤ 2·10⁵ is verified by two independent programs (FiniteCheck*.md). The ranges n ≤ 2·10⁵ and n ≥ 10⁵ overlap.

---

## 10. Minor-arc lemma (ws_mech_minor.py / ws_mech_cd.py) — written proof, for completeness

This item was not on the trust list, because the scripts carry the derivation in their docstrings. It is written out here so that this file is
self-contained.

**Claim.** Let θ ≤ θ_max, (k+1)θ ≥ u_top and φ ∈ [c₀θ, π]. Then θS(φ) ≥ PSI, where S = −2log|χ| = Σ_{i≤k}log(1 + sin²(iφ/2)/sinh²(iθ/2)).
Each summand is ≥ 0 (§6, δ_C step 1).

**(a) c = φ/θ in a cell [c_lo, c_hi] ⊂ [c₀, c_a].**
1. Take T ∈ (0, π/2] and s = sinT/T. For i ≤ 2T/φ we have iφ/2 ≤ T, so sin(iφ/2) ≥ s·iφ/2 (E3).
2. The summand is then ≥ log(1 + s²c_lo²f₂(iθ)) =: l(iθ), using c ≥ c_lo. The function l is decreasing (E2).
3. Drop the other summands, which are ≥ 0, and take I = min(k, ⌊2T/φ⌋). Then θS ≥ θΣ_{i≤I}l(iθ) ≥ ∫_θ^{(I+1)θ}l.
4. (I+1)θ ≥ min((k+1)θ, 2T/c) ≥ min(u_top, 2T/c_hi), and θ ≤ θ_max. Since l ≥ 0, ∫_θ^{(I+1)θ}l ≥ ∫_{θ_max}^{min(u_top, 2T/c_hi)}l.
5. A right Riemann sum of a decreasing function is a lower bound for its integral.

**(D) φ ∈ [c_aθ, π], d = φ/2.**
1. By E4, the summand is ≥ sin²(id)·w_i with w_i = 2 log coth(iθ/2), which is decreasing (E8).
2. By Abel summation (E16), Σ sin²(id)w_i = Σ_{I≤k}G_I(w_I − w_{I+1}), with w_{k+1} := 0 and G_I = Σ_{i≤I}sin²(id).
3. By E9, G_I ≥ I/2 − 1/(2 sin d). Also G_I ≥ 0. So G_I ≥ (I − N₀)₊/2 with N₀ = ⌈1/sin d⌉.
4. The coefficients w_I − w_{I+1} are ≥ 0, so S ≥ ½Σ_{i>N₀}w_i ≥ ½θ^{−1}∫_{(N₀+1)θ}^{(k+1)θ}w.
5. (N₀+1)θ ≤ θ/sin d + 2θ. Also sin d ≥ sin(c_aθ/2) ≥ (c_aθ/2)sinc(c_aθ_max/2), because sin is increasing on [0, π/2], sinc is decreasing (E3),
   and c_aθ_max/2 ≤ π/2 is asserted. So (N₀+1)θ ≤ L := 2θ_max + 2/(c_a·sinc(c_aθ_max/2)).
6. Hence θS ≥ ∫_L^{u_top}log coth(u/2)du. If L ≥ u_top the bound is 0, which is still valid.

**Result.** PSI is the min of the (a) cells and (D). It depends on θ only through θ ≤ θ_max and (k+1)θ ≥ u_top. □

---

## 11. Findings of this pass

- **No genuine gap was found.** Every derivation on the trust list was re-derived and holds as stated.
- Two improvements over the previous text:
  - **(i)** "Ψ(s) increasing" (D3) now has an exact closed-form derivative, Ψ′ = (γs² + 2s + γ²)/((s²−γ)(γs+1)) + γ⁴s/(s²−γ)² > 0. It is no
    longer a hand argument; the old wording "2s(γs+1) > γ(s²−γ)" is the same fact.
  - **(ii)** (B3) is shown not to be load-bearing (§9.3).
- Small text corrections, none of which affects the proof:
  - LemmaB §4 writes ν′ ≤ ρν for ν_{k+1} ≤ ρν_k. The prime there means the k+1 point, not ν′ = ν + log θ.
  - In B.3 the "+1" in (ν/θ + 1)dθ is slack.
  - LemmaCD's statement "Ψ is increasing in s, because 2s(γs+1) > γ(s²−γ)" omitted that the second term of Ψ is also increasing; §7 covers both.
- **Remaining trust base:** python-flint/Arb 0.9.0, sympy 1.14, mpmath (sanity tests only), gcc 13.3 and the two C finite-check programs. Also the
  classical facts that are quoted, not re-proved:
  - Taylor's theorem with integral remainder;
  - the mean value theorem;
  - term-by-term differentiation of uniformly convergent series;
  - additivity of total variation;
  - the identity theorem for entire functions.
