"""Certified Landau root of eps(k, omega) = 1 + (1 + zeta Z(zeta))/k^2 = 0, Maxwellian.

Z(zeta) = i sqrt(pi) exp(-zeta^2) erfc(-i zeta)  (entire; the Landau continuation).
Newton in acb, then a Krawczyk test proves a unique zero inside the returned ball.
"""
import flint
from flint import acb, arb


def _f(z, g):
    Z = acb.pi().sqrt() * acb(0, 1) * (-z * z).exp() * (acb(0, -1) * z).erfc()
    u = 1 + z * Z                       # 1 + zeta Z
    f = 1 + u / g
    fp = (Z - 2 * z * u) / g            # d/dzeta (1 + zeta Z) = Z + zeta Z',  Z' = -2u
    return f, fp


def landau_root(k, prec_bits=400, omega_guess=None):
    """Return (omega ball, certified flag) for the least-damped root at wavenumber k."""
    old = flint.ctx.prec
    kf = float(arb(k).mid()) if isinstance(k, arb) else float(k)
    # guard bits: exp(-zeta^2) * erfc(-i zeta) cancels ~ zeta^2/ln2 bits
    flint.ctx.prec = prec_bits + 64 + int(2.0 / (kf * kf))
    try:
        k = arb(k)
        g = k * k
        if omega_guess is None:
            w0 = (1 + 3 * g).sqrt()
            omega_guess = acb(w0, 0)
        z = acb(omega_guess) / (k * arb(2).sqrt())
        for _ in range(200):
            f, fp = _f(acb(z.mid()), g)
            dz = f / fp
            z = acb((z - dz).mid())
            if dz.abs_upper() < arb(2) ** (-prec_bits - 8):
                break
        # Krawczyk: K(B) = m - f(m)/f'(m) + (1 - f'(B)/f'(m)) (B - m) subset of int(B)
        r = arb(2) ** (-prec_bits + 10)
        B = acb(arb(z.real.mid(), r), arb(z.imag.mid(), r))
        fm, fpm = _f(z, g)
        _, fpB = _f(B, g)
        K = z - fm / fpm + (1 - fpB / fpm) * (B - z)
        cert = (B.real.contains_interior(K.real) and B.imag.contains_interior(K.imag))
        omega = K * k * arb(2).sqrt()
        return omega, cert
    finally:
        flint.ctx.prec = old


def landau_root_cont(k, prec_bits=400):
    """Least-damped root at k, reached by continuation from k = 0.3 (Bohm-Gross start) in
    steps of 0.05 at low precision, then refined and Krawczyk-certified at prec_bits."""
    kf = float(k)
    guess, kk = None, 0.3
    while kk < kf - 1e-12:
        guess, _ = landau_root(repr(kk), 120, guess)
        kk = min(kk + 0.05, kf)
    return landau_root(k, prec_bits, guess)


if __name__ == "__main__":
    guess = None
    for k in ["0.1", "0.2", "0.3", "0.5", "0.7", "0.85", "1", "1.5", "2"]:
        w, c = landau_root(k, 200, guess)
        guess = w
        print(k, c, w.str(30))
