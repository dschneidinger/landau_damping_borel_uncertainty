"""Resurgence test 3: the all-orders Landau damping series.

    gamma(k) = Im omega = -sqrt(pi/8) k^{-3} exp(-1/(2k^2) - 3/2) * Psi(k^2) + O(e^{-3/(2k^2)}),
    Psi(g)   = -Phi1(g)/sqrt(W0(g)) = (W0 - g W0') exp(-(W0 - 1 - 3g)/(2g)) / ... (rational series)

Psi is itself Gevrey-1 (Borel singularity at 1/2); optimally truncated it should match the
certified root to relative accuracy ~ e^{-1/(2k^2)}.
"""
import flint
from flint import fmpq, fmpq_series
import mpmath as mp
from series import bohm_gross, phi1
from root import landau_root


def psi_coeffs(M):
    a = bohm_gross(M + 3)
    f = phi1(a, M + 1)
    flint.ctx.cap = M + 1
    F = fmpq_series(f[:M + 1], prec=M + 1)
    W = fmpq_series([fmpq(c) for c in a[:M + 1]], prec=M + 1)
    P = -F * W.rsqrt()
    co = P.coeffs()
    return [co[n] if n < len(co) else fmpq(0) for n in range(M)]


if __name__ == "__main__":
    M = 120
    d = psi_coeffs(M)
    print("Psi(g) = " + " + ".join(f"({c}) g^{n}" for n, c in enumerate(d[:8])) + " + ...")
    mp.mp.dps = 60
    dq = [mp.mpf(int(c.p)) / int(c.q) for c in d]
    for k in ["0.08", "0.1", "0.12", "0.15", "0.2"]:
        om, cert = landau_root(k, 260)
        flint.ctx.prec = 300
        gam = mp.mpf(om.imag.mid().str(70, radius=False))
        kk = mp.mpf(k); g = kk * kk
        pref = -mp.sqrt(mp.pi / 8) / kk ** 3 * mp.exp(-1 / (2 * g) - mp.mpf(1.5))
        terms = [dq[n] * g ** n for n in range(M)]
        nopt = min(range(1, M), key=lambda n: abs(terms[n]))
        approx = pref * mp.fsum(terms[:nopt])
        landau = pref                                   # textbook formula with the -3/2
        print(f"k={k}: gamma={mp.nstr(gam, 15)} cert={cert} | textbook rel.err {mp.nstr(abs(landau/gam-1), 3)}"
              f" | all-orders sigma=1 (n_opt={nopt}) rel.err {mp.nstr(abs(approx/gam-1), 3)}"
              f" | sigma=2 ('factor 2') rel.err {mp.nstr(abs(2*approx/gam-1), 3)}"
              f" | e^(-1/2k^2)={mp.nstr(mp.exp(-1/(2*g)), 3)}")
