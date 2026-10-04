"""All transseries sectors of the Landau root, exactly, and the resurgence tower test.

Equation (g = k^2, W = omega^2): G(W) - sigma K(W) = 1 with
    G(W) = P(g/W)/g,  P(x) = sum (2n+1)!! x^{n+1}             (asymptotic part of Z'/2g)
    K(W) = i sqrt(pi) (W/(2g))^{1/2} e^{-W/(2g)} / g          (Stokes term of Z'/2g)
Put W = W0 + g Z and eps' = sigma K(W0)/g. Then Z solves Z = eps' phi(Z) with
    phi(Z) = kappa(Z) / sum_{j>=1} G_j g^{j-1} Z^{j-1},   kappa(Z) = (1 + gZ/W0)^{1/2} e^{-Z/2},
    G_j = [Y^j] G(W0 + Y),
and every coefficient is a power series in g with rational coefficients.  Lagrange:
    z_m = [eps'^m] Z = (1/m) [Z^{m-1}] phi(Z)^m.
Sector m of the transseries is
    u_m = sigma^m c_m T^m g^{1 - 5m/2} Psi_m(g),   c_m = (i sqrt pi)^m e^{-3m/2} 2^{-m/2},  T = e^{-1/(2g)},
    Psi_m = W0^{m/2} E^m z_m,   E = exp(-(W0 - 1 - 3g)/(2g))          (rational power series)
Psi_1 = Phi1 of series.py.

Stokes automorphism (one jump of Z, sigma -> sigma + 2):  S_+[W(sigma)] = S_-[W(sigma + 2)], so
    Disc u_m = sum_{j>=1} C(m+j, m) 2^j u_{m+j}.
Leading term => large-order law for the coefficients psi_{m,h} of Psi_m (m >= 1):
    psi_{m,h} ~ s (m+1) (4 e^{-3/2}/sqrt(pi)) sum_j psi_{m+1,j} 2^{h-j} Gamma(h + 5/2 - j)
and for the Bohm-Gross integers (m = 0):
    a_n       ~ s       (2 e^{-3/2}/sqrt(pi)) sum_j psi_{1,j}   2^{n-j} Gamma(n + 3/2 - j).
The bridge multipliers 1, 2, 3, ... (times the common 2) are fixed by the single Stokes jump of Z;
nothing is fitted.
"""
import sys, time
import flint
from flint import fmpq, fmpq_series
from series import bohm_gross, double_factorials


def _ser(c, L):
    return fmpq_series([fmpq(v) for v in c[:L]], prec=L)


def _coeffs(s, L):
    c = s.coeffs()
    return [c[i] if i < len(c) else fmpq(0) for i in range(L)]


def sectors(L, M):
    """Return (a, [Psi_1..Psi_M]) with L g-coefficients each (exact rationals)."""
    flint.ctx.cap = L + 2
    a = bohm_gross(L + 2)
    W0 = _ser(a, L)
    invW = W0.inv()
    p = double_factorials(L + 1)                  # P(x) = sum p_n x^{n+1}
    # G_j g^{j-1} = (1/g^... ) : G(W0+Y) = P(x/(1+Y/W0))/g.
    #  [Y^j] = (1/g) sum_n p_n x^{n+1} binom(-(n+1), j) W0^{-j}
    # We need A_j := G_j g^{j-1} = g^{j-2} W0^{-j} sum_n p_n binom(-(n+1), j) x^{n+1}.
    # x^{n+1} = g^{n+1} * (1/W0)^{n+1}; so the sum is g * S_j with S_j a power series.
    A = []
    xg = invW                                     # x/g
    for j in range(1, M + 1):
        acc = [fmpq(0)] * L
        powx = fmpq_series([1], prec=L)
        for n in range(L):
            powx = powx * xg if n else xg          # (x/g)^{n+1}
            b = fmpq(1)
            for i in range(j):                     # binom(-(n+1), j)
                b *= fmpq(-(n + 1) - i, i + 1)
            cs = _coeffs(powx, L)
            # term = p_n * b * g^{n+1} * (x/g)^{n+1}
            for t in range(L - (n + 1)):
                acc[t + n + 1] += p[n] * b * cs[t]
        # sum = sum_t acc[t] g^t (valuation >= 1) ; A_j = g^{j-2} W0^{-j} * sum
        # valuation of sum: 1  => A_j has valuation j-1 >= 0
        shift = j - 2
        if shift >= 0:
            ser = [fmpq(0)] * shift + acc[:L - shift]
        else:
            ser = acc[1:] + [fmpq(0)]
        A.append(_ser(ser, L) * invW ** j)
    # denominator D(Z) = sum_{j>=1} A_j Z^{j-1};  A_1 = G_1 = -1/R
    # kappa(Z) = (1 + gZ/W0)^{1/2} e^{-Z/2} = sum_k kap_k Z^k
    gW = fmpq_series([0] + _coeffs(invW, L - 1), prec=L)     # g/W0
    kap = []
    for k in range(M + 1):
        # [Z^k] (1+gZ/W0)^{1/2} e^{-Z/2} = sum_{i+l=k} binom(1/2,i)(g/W0)^i (-1/2)^l/l!
        s = fmpq_series([0], prec=L)
        for i in range(k + 1):
            l = k - i
            b = fmpq(1)
            for r in range(i):
                b *= fmpq(1, 2) - r
                b /= r + 1
            fl = fmpq(1)
            for r in range(1, l + 1):
                fl *= r
            s = s + gW ** i * (b * fmpq(-1, 2) ** l / fl)
        kap.append(s)
    # phi = kappa / D as a polynomial in Z with series coefficients, degree < M
    Dinv = [A[0].inv()]
    for k in range(1, M):                          # inverse of D(Z) in Z, up to Z^{M-1}
        s = fmpq_series([0], prec=L)
        for j in range(1, k + 1):
            if j < len(A):
                s = s + A[j] * Dinv[k - j]
        Dinv.append(-s * Dinv[0])
    phi = []
    for k in range(M):
        s = fmpq_series([0], prec=L)
        for i in range(k + 1):
            s = s + kap[i] * Dinv[k - i]
        phi.append(s)

    def zmul(u, v):
        out = [fmpq_series([0], prec=L) for _ in range(M)]
        for i in range(M):
            for j in range(M - i):
                out[i + j] = out[i + j] + u[i] * v[j]
        return out

    E = fmpq_series([fmpq(0)] + [-fmpq(a[n + 1]) / 2 for n in range(1, L)], prec=L).exp()
    sqW = W0.sqrt()
    Psi, pw = [], [fmpq_series([1], prec=L)] + [fmpq_series([0], prec=L)] * (M - 1)
    for m in range(1, M + 1):
        pw = zmul(pw, phi)
        z_m = pw[m - 1] * fmpq(1, m)
        Psi.append(_coeffs(sqW ** m * E ** m * z_m, L))
    return a, Psi


if __name__ == "__main__":
    import mpmath as mp
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    t0 = time.time()
    a, Psi = sectors(L, M)
    print(f"sectors: L={L}, M={M}, {time.time()-t0:.1f}s")
    for m, P in enumerate(Psi, 1):
        print(f"Psi_{m} = " + ", ".join(str(c) for c in P[:5]) + ", ...")
    mp.mp.dps = int(0.4 * L) + 60
    q = lambda c: mp.mpf(int(c.p)) / int(c.q)
    e32 = mp.e ** mp.mpf(-1.5) / mp.sqrt(mp.pi)
    levels = [("a_n  (m=0)", [mp.mpf(int(v)) for v in a[:L]], [q(c) for c in Psi[0]], 2 * e32, mp.mpf(1.5))]
    for m in range(1, M):
        levels.append((f"Psi_{m} (m={m})", [q(c) for c in Psi[m - 1]], [q(c) for c in Psi[m]],
                       (m + 1) * 4 * e32, mp.mpf(2.5)))
    for name, tgt, nxt, C, beta in levels:
        def pred(n, H):
            return -C * mp.fsum(nxt[j] * mp.mpf(2) ** (n - j) * mp.gamma(n + beta - j) for j in range(H))
        row = []
        for n in [L // 8, L // 4, L // 2 - 1]:
            best = min((abs(tgt[n] / pred(n, H) - 1), H) for H in range(1, n // 2 + 2))
            row.append(f"n={n}: {mp.nstr(best[0], 3)} (2^-n {mp.nstr(mp.mpf(2) ** -n, 2)})")
        # negative control: the same law with the bridge multiplier shifted by one
        n = L // 2 - 1
        m = levels.index((name, tgt, nxt, C, beta))
        Cw = C * (m + 2) / (m + 1) if m else C * 2
        wrong = min(abs(tgt[n] / (Cw / C * pred(n, H)) - 1) for H in range(1, n // 2 + 2))
        print(f"{name:16s} <- sector {m+1}, Stokes multiplier {2*(m+1)}:  " + " | ".join(row)
              + f" || control, multiplier shifted: {mp.nstr(wrong, 3)}", flush=True)
