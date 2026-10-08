# Finite check (Workstream 4): ACG unimodality for all n ≤ 200 000

**Result.** For every 1 ≤ n ≤ 200 000 the row p(n,1),…,p(n,n) is weakly unimodal. PASS.
Covers the analytic threshold N_A = 100 000 with a 2× margin. Zero undetermined pairs, so no
exact resolution was needed.

## 1. Reduction to k ≤ ⌈n/2⌉ + 1
**Lemma E.** If 2k ≥ n then p(n,k) = p(n−k).
*Proof.* Subtracting 1 from each part is a bijection from partitions of n into exactly k parts onto
partitions of n−k into at most k parts. Since n−k ≤ k, every partition of n−k has at most n−k ≤ k
parts, so the constraint does nothing and the count is p(n−k). ∎

So for c = ⌈n/2⌉ and k ≥ c, p(n,k) = p(n−k) ≥ p(n−k−1) = p(n,k+1), because p is nondecreasing
(adding a part 1 injects partitions of j into partitions of j+1). The tail k ≥ c is nonincreasing.
The checker therefore certifies the signs of d_k = p(n,k+1) − p(n,k) for 1 ≤ k ≤ c, and the row is
unimodal iff no '+' comes after a '−' among those k (the tail adds only '−'/'0').

## 2. Algorithm (ws_fc_check.c)
v[m] = p_{≤k}(m), updated in place for k = 1,2,…: v[m] += v[m−k] for m = k..N−k (increasing m).
After step k, p(n,k) = v[n−k] for every n ≤ N at once. Only m ≤ N−k is updated at step k. This is
enough because later steps only read indices ≤ N−k' < N−k. Time O(N²/2), memory O(N).
Per n, the previous value p(n,k−1) is kept and the pair is classified. k runs to ⌈N/2⌉+1.

## 3. Rigorous error argument
* Arithmetic is x87 `long double`: 64-bit significand, exponent up to 2^16383. log2 p(200000) ≈ 1660,
  so **no scaling is needed** and nothing overflows or underflows. All values are ≥ 1 or exactly 0.
* Two copies are kept. lo[] is computed by a thread whose rounding mode is `fesetround(FE_DOWNWARD)`,
  hi[] by a thread with `FE_UPWARD`. glibc's fesetround sets both the x87 control word and MXCSR.
  The build uses `-frounding-math`, and the array is accessed through a `volatile` pointer so the
  operations are not reordered or folded.
* Induction: v[0] = 1 is exact. Every update is a single addition of two nonnegative numbers.
  With a ≥ A ≥ 0 and b ≥ B ≥ 0, fl_down(a+b) ≤ A+B. The upward case is symmetric. Hence
  lo[m] ≤ p_{≤k}(m) ≤ hi[m] at every step. No subtraction is ever performed, so there is no cancellation.
* Classification of d_k = a_{k+1} − a_k, using intervals [l_k, h_k] ∋ a_k:
  '+' if l_{k+1} > h_k; '−' if h_{k+1} < l_k; '0' only if both intervals are the same single point
  (so the values are exactly known and equal); otherwise '?' (written out for exact resolution).
  Comparisons are exact, so they do not depend on the rounding mode.
* A row is flagged as a violation if a certified '+' follows a certified '−'. '?' pairs would have been
  resolved by `ws_fc_ref.py resolve` (exact Python big integers). There were none.
* Reported mode = the largest k with a certified increase into k (1 if there is none). Under unimodality
  it is the first index of the maximum.

## 4. Validation
`ws_fc_ref.py ref 2000` computes p(n,k) exactly with Python big integers via p(n,k) = p(n−1,k−1) + p(n−k,k).
It also asserts Lemma E and the nonincreasing tail for all n ≤ 2000. Its per-n output (n, mode, violation,
#undetermined) is **byte-identical** to the C log for N = 2000. Sign counts agree: +151295 and 8 exact
ties. The '−' count differs by 1 because the trivial n=1 pair (k=1, c=1, p(1,2)=0) is counted only by the
reference.

## 5. Results (2-core x86_64, gcc 13.3, OpenMP; GMP not installed and not needed)
| N | runtime | '+' | '−' | exact ties | '?' | violations |
|---|---|---|---|---|---|---|
| 2 000 | 0.02 s | 151 295 | 849 696 | 8 | 0 | 0 |
| 100 000 | 25 s | 85 086 401 | 2 414 963 590 | 8 | 0 | 0 |
| 200 000 | 103 s | 256 705 214 | 9 743 394 777 | 8 | 0 | 0 |

The 8 ties all come from small n, where both values are certified exact points.
Sample modes: n=40000 → 789, 80000 → 1191, 120000 → 1513, 160000 → 1792, 200000 → 2042.
The full per-n mode is in log*.txt (columns: n mode violation #undetermined).

SHA-256 (also in SHA256SUMS):
```
675740a5c82f111a986517bef17ce56ed9b9129665fbc98d6fdff696c4d89e2d  log100k.txt
45ebb662ffbcd5a46af32fe1b1133c78a30e9e5167080a9a148d67c3685e3132  log200k.txt
67e92fa044d5bcc472b87d94cfd213d040ba9a7641cb77b4c4d9986060b4936d  ws_fc_check.c
95fb9bd1882ee69fe9679e338134842ec18e6c3abafcdbafae0b72ed684a5d4c  ws_fc_ref.py
```

## 6. Commands
```
gcc -O2 -fopenmp -frounding-math -o ws_fc_check ws_fc_check.c -lm
./ws_fc_check 2000 und2k.txt log2k.txt;  python3 ws_fc_ref.py ref 2000 ref2k.txt;  cmp ref2k.txt log2k.txt
./ws_fc_check 100000 und100k.txt log100k.txt
./ws_fc_check 200000 und200k.txt log200k.txt
python3 ws_fc_ref.py resolve und200k.txt   # empty: nothing to resolve
```
Caveat (trust base): the proof relies on the correctness of x86 x87 directed rounding and of glibc
fesetround in OpenMP threads. The rounding mode is set per thread, and a hard-coded check (§4) cross-validates
the output. An independent second implementation, e.g. MPFR or exact-integer for N up to about 2e4, would
further reduce that trust base.
