"""Resurgence test 1: large-order growth of the integer Bohm-Gross coefficients a_n is
predicted, with no fitting, by the one-instanton (Landau-damping) series Phi1:

    a_n  ~  s * (2 e^{-3/2}/sqrt(pi)) * sum_h f_h 2^{n-h} Gamma(n + 3/2 - h)

(s = overall orientation sign, fixed by the leading term), with relative error
O(2^{-n}) from the two-instanton Borel singularity at t = 1.
"""
import sys, time
import mpmath as mp
from series import bohm_gross, phi1

N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
t0 = time.time()
a = bohm_gross(N + 1)
f = phi1(a, N // 2 + 2)
print(f"series: N={N}, {time.time()-t0:.1f}s", flush=True)

mp.mp.dps = int(0.6 * N) + 60
C = 2 * mp.e ** mp.mpf(-1.5) / mp.sqrt(mp.pi)
fq = [mp.mpf(int(v.p)) / int(v.q) for v in f]


def predict(n, H):
    return C * mp.fsum(fq[h] * mp.mpf(2) ** (n - h) * mp.gamma(n + mp.mpf(1.5) - h)
                       for h in range(H))


s = mp.sign(mp.mpf(int(a[N])) / predict(N, 1))
print("orientation sign s =", int(s))
for n in [N // 8, N // 4, N // 2, N]:
    exact = mp.mpf(int(a[n]))
    best = None
    for H in range(1, n // 2 + 2):
        rel = abs(exact / (s * predict(n, H)) - 1)
        if best is None or rel < best[1]:
            best = (H, rel)
    lead = abs(exact / (s * predict(n, 1)) - 1)
    print(f"n={n:5d}  a_n has {len(str(int(a[n])))} digits | leading-term rel.err {mp.nstr(lead,3)}"
          f" | optimal H={best[0]} rel.err {mp.nstr(best[1],3)}  (2^-n = {mp.nstr(mp.mpf(2)**-n,3)})")
