# Finite check, independent exact verifier: ACG unimodality for all n ≤ 200 000 (no floating point)

**Result: PASS.** For every 1 ≤ n ≤ 200 000 the row p(n,1),…,p(n,n) is weakly unimodal. Every
comparison is exact integer arithmetic. The program contains no floating point, no rounding modes and no threads.

## Design (`ws_fc_exact.c`)
* Bignums: fixed-width slots of L 64-bit limbs per index, with a per-entry used-length. L = 27 for N = 2e5,
  from p(N) < exp(π√(2N/3)) plus slack. A carry past L aborts with exit code 2, so the bound is never trusted silently.
  Addition uses `__builtin_add_overflow`. Comparison checks the length first, then limbs from the top.
* Recurrence: this is the same column recurrence as `ws_fc_check`. After step k, v[m] = p_{≤k}(m) = p(m+k,k).
  At step k, m runs upward and v[m] += v[m−k]. Right after v[m] is updated, v[m+1] has not been updated yet,
  so it still holds p(n,k−1) with n = m+k. The comparison p(n,k) vs p(n,k−1) therefore needs no second array.
  Entries with m < k do not change at step k and are compared directly.
* Range: the pairs (k−1,k) with 2 ≤ k and k−1 ≤ ⌈n/2⌉ are checked, which is the reduction of FiniteCheck.md §1.
  A violation is a '+' after a '−'. The mode is the last '+' index (1 if there is none).
* Output: the log has the same per-n format as `ws_fc_check` (`n mode viol #und`, with #und ≡ 0).
* Memory: (N+2)·L·8 B ≈ 43 MB at N = 2e5. The program is single-threaded. Time is about N²·(avg limbs)/2.

## Validation
* N = 2000: the log is byte-identical to `ref2k.txt` (Python big ints, `ws_fc_ref.py`) and to `log2k.txt`.
* N = 20000, 100000, 200000: the logs are byte-identical to the x87 interval-arithmetic logs
  (`exlog100k.txt` ≡ `log100k.txt`, `exlog200k.txt` ≡ `log200k.txt`, with the same SHA-256 values).
  The sign counts agree exactly with FiniteCheck.md §5.

## Results (2-core x86_64, gcc 13.3, single thread)
| N | limbs | time | '+' | '−' | ties | violations |
|---|---|---|---|---|---|---|
| 2 000 | 4 | 0.01 s | 151 295 | 849 696 | 8 | 0 |
| 20 000 | 10 | 1.5 s | 6 437 567 | 93 572 424 | 8 | 0 |
| 100 000 | 20 | 71 s | 85 086 401 | 2 414 963 590 | 8 | 0 |
| 200 000 | 27 | 8 m 18 s | 256 705 214 | 9 743 394 777 | 8 | 0 |

SHA-256:
```
6f0a4aad6d1c7b01f3a3a80fa174afd292e4eebd8012517e3e32d953a0815e68  ws_fc_exact.c
675740a5c82f111a986517bef17ce56ed9b9129665fbc98d6fdff696c4d89e2d  exlog100k.txt
45ebb662ffbcd5a46af32fe1b1133c78a30e9e5167080a9a148d67c3685e3132  exlog200k.txt
```
The SHA-256 values of the exact logs equal those of `log100k.txt` and `log200k.txt`. The two verifiers are independent
(one uses exact integers, the other x87 directed-rounding intervals) and agree bit-for-bit on every n's mode and verdict.

## Commands
```
gcc -O3 -march=native -Wall -o ws_fc_exact ws_fc_exact.c
./ws_fc_exact 2000 ex2000.txt && cmp ex2000.txt ref2k.txt
./ws_fc_exact 100000 exlog100k.txt && cmp exlog100k.txt log100k.txt
./ws_fc_exact 200000 exlog200k.txt && cmp exlog200k.txt log200k.txt
```
Exit code: 0 means PASS, 3 means a violation was found, 2 means limb overflow.
Trust base: gcc integer code generation and the shared reduction lemma (FiniteCheck.md §1).
