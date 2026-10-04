"""Exact small-k series for the Maxwellian Langmuir (Landau) root.

Units: omega_p = lambda_D = v_th = 1, g = k^2, W = omega^2, zeta = omega/(sqrt(2) k).
Dispersion relation: eps = 1 - Z'(zeta)/(2 g) = 0, Z the plasma dispersion function.

Perturbative sector.  With x = g/W the asymptotic part of Z' gives  P(x) = g,
P(x) = sum_{n>=0} (2n+1)!! x^{n+1}.  So x(g) = P^{-1}(g) (series reversion, integer
coefficients) and W0(g) = g/x(g) = 1 + 3g + 6g^2 + ...  (integers).

One-instanton sector (closed form, derived in NOTES.md):
    W = W0 + sigma * i*sqrt(pi/2) * e^{-3/2} * g^{-3/2} e^{-1/(2g)} * Phi1(g) + O(e^{-1/g})
    Phi1 = -sqrt(W0) * (W0 - g W0') * exp(-(W0 - 1 - 3g)/(2g))        (rational series)
"""
import flint
from flint import fmpq, fmpq_series, fmpz


def double_factorials(N):
    """(2n+1)!! for n = 0..N-1."""
    out, c = [], fmpz(1)
    for n in range(N):
        c = c * (2 * n + 1) if n else fmpz(1)
        out.append(c)
    return out


def bohm_gross(N):
    """Integer coefficients a_0..a_{N-1} of W0(g) = omega^2 = sum a_n g^n."""
    flint.ctx.cap = N + 2
    df = double_factorials(N + 1)
    P = fmpq_series([0] + [fmpq(d) for d in df[:N + 1]], prec=N + 2)
    x = P.reversion()                                  # x(g) = g - 3g^2 + ...
    xg = fmpq_series(x.coeffs()[1:], prec=N + 1)       # x/g = 1 - 3g + ...
    W = xg.inv()
    a = [W.coeffs()[n] if n < len(W.coeffs()) else fmpq(0) for n in range(N)]
    assert all(c.q == 1 for c in a), "Bohm-Gross coefficients must be integers"
    return [c.p for c in a]


def phi1(a, M):
    """Rational coefficients f_0..f_{M-1} of Phi1(g) from W0 coefficients a (len >= M+2)."""
    flint.ctx.cap = M + 1
    W = fmpq_series([fmpq(c) for c in a[:M + 2]], prec=M + 1)
    Wp = W.derivative()
    gWp = fmpq_series([0] + Wp.coeffs()[:M], prec=M + 1)
    # (W - 1 - 3g)/(2g): drop first two coefficients, shift down by 2... careful: (W-1-3g)/g
    c = [fmpq(a[n + 1]) / 2 for n in range(M + 1)]   # (W-1)/(2g) = 3/2 + 3g + ...
    c[0] = fmpq(0)                                    # subtract the 3/2 (pulled out as e^{-3/2})
    E = fmpq_series([-v for v in c], prec=M + 1).exp()
    F = -(W.sqrt()) * (W - gWp) * E
    co = F.coeffs()
    return [co[n] if n < len(co) else fmpq(0) for n in range(M)]


if __name__ == "__main__":
    a = bohm_gross(16)
    print("W0 coefficients:", [int(v) for v in a])
    f = phi1(bohm_gross(30), 10)
    print("Phi1 coefficients:", [str(v) for v in f])
