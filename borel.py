"""Resurgence test 2: the lateral Borel sum of the integer Bohm-Gross series IS the
complex Landau root.  Claim:  omega_L(k)^2 = S_{theta}[W0](k^2) for any ray theta in (-pi/2, 0)
(the damped root); rays in (0, pi/2) give its complex conjugate.  Checked against the
Krawczyk-certified Arb root (root.py).

Borel transform B(t) = sum c_n t^n / n!  (singularities at t = 1/2, 1, ... on R+),
Pade-approximated in the conformal variable w with t = 2w/(1+w)^2 (cut [1/2,inf) -> |w|=1),
Laplace integral along the ray arg s = theta.
"""
import sys
import mpmath as mp


def borel_pade(c):
    """Return lateral(g, theta) for the formal series sum c[n] g^n (mpf coefficients)."""
    N = len(c)
    b = [c[n] / mp.factorial(n) for n in range(N)]
    tw = mp.taylor(lambda w: 2 * w / (1 + w) ** 2, 0, N - 1)
    Bw = [mp.mpf(0)] * N
    pw = [mp.mpf(1)] + [mp.mpf(0)] * (N - 1)          # t(w)^n
    for n in range(N):
        for j in range(n, N):
            Bw[j] += b[n] * pw[j]
        pw = [mp.mpf(0)] + [mp.fsum(pw[i] * tw[j - i] for i in range(j)) for j in range(1, N)]
    L = (N - 1) // 2
    P, Q = mp.pade(Bw, L, N - 1 - L)
    P, Q = P[::-1], Q[::-1]

    def B(t):
        s = mp.sqrt(1 - 2 * t)
        w = (1 - s) / (1 + s)
        return mp.polyval(P, w) / mp.polyval(Q, w)

    def lateral(g, theta):
        e = mp.expjpi(mp.mpf(theta) / mp.pi)
        return mp.quad(lambda r: mp.exp(-r * e) * B(g * r * e) * e, [0, 1, 4, 16, 64, mp.inf])

    return lateral


def landau_W(k, dps, guess=None):
    """Certified omega_L^2 at wavenumber k as an mpc (plus the omega ball for continuation)."""
    import flint
    from root import landau_root_cont
    om, cert = landau_root_cont(k, int(3.4 * dps) + 64)
    assert cert, f"Krawczyk certification failed at k={k}"
    flint.ctx.prec = 4 * dps
    W = om * om
    return mp.mpc(mp.mpf(W.real.mid().str(dps + 10, radius=False)),
                  mp.mpf(W.imag.mid().str(dps + 10, radius=False))), om


if __name__ == "__main__":
    from series import bohm_gross
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    thetas = sys.argv[3].split(",") if len(sys.argv) > 3 else ["-0.25", "-0.7", "-1.2"]
    ks = sys.argv[4].split(",") if len(sys.argv) > 4 else ["0.2", "0.3", "0.5", "0.7", "1", "1.5", "2"]
    lateral = borel_pade([mp.mpf(int(v)) for v in bohm_gross(N)])
    print(f"N={N} terms, dps={mp.mp.dps}, rays theta={thetas}")
    guess = None
    for k in ks:
        Wc, guess = landau_W(k, mp.mp.dps, guess)
        g = mp.mpf(k) ** 2
        errs = [abs(lateral(g, th) - Wc) for th in thetas]
        print(f"k={k:5s} omega_L^2={mp.nstr(Wc, 20)}  |S_theta[W0]-omega_L^2|: "
              + "  ".join(f"{th}:{mp.nstr(e, 3)}" for th, e in zip(thetas, errs)), flush=True)
