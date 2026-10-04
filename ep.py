"""Certified exceptional points (EPs) of the Maxwellian Langmuir dispersion relation, and
whether the Landau root itself runs into them.

eps(zeta; g) = 1 - Z'(zeta)/(2g).  A double root (EP) needs eps = d eps/d zeta = 0, i.e.
    Z''(zeta*) = 0,   g* = k*^2 = Z'(zeta*)/2,   omega*^2 = 2 g* zeta*^2.
Z' = -2(1 + zeta Z),  Z'' = -2(Z + zeta Z'),  Z''' = -2(2Z' + zeta Z'').
Each zeta* is certified unique in its ball by a Krawczyk test on Z''.
Sheet test: continue the least-damped root from g = 0.01 along the real axis to |g*|, then
along the arc |g| = |g*| to arg g* (and to arg g* +/- 0.02), and compare with zeta*.
"""
import json
import flint
from flint import acb, arb


def Zs(z):
    Z = acb.pi().sqrt() * acb(0, 1) * (-z * z).exp() * (acb(0, -1) * z).erfc()
    Z1 = -2 * (1 + z * Z)
    Z2 = -2 * (Z + z * Z1)
    Z3 = -2 * (2 * Z1 + z * Z2)
    return Z, Z1, Z2, Z3


def certify(z0, bits=300):
    flint.ctx.prec = bits + 64
    z = acb(z0)
    for _ in range(100):
        _, _, Z2, Z3 = Zs(acb(z.mid()))
        z = acb((z - Z2 / Z3).mid())
    r = arb(2) ** (-bits + 20)
    B = acb(arb(z.real.mid(), r), arb(z.imag.mid(), r))
    _, _, Z2m, Z3m = Zs(z)
    _, _, _, Z3B = Zs(B)
    K = z - Z2m / Z3m + (1 - Z3B / Z3m) * (B - z)
    ok = B.real.contains_interior(K.real) and B.imag.contains_interior(K.imag)
    _, Z1, _, _ = Zs(K)
    g = Z1 / 2
    return K, g, ok


def asym_seed(n, it=60):
    """EP asymptotics: Z'' = 0 balances -2/zeta^3 against sigma*i*sqrt(pi)*4 zeta^2 e^{-zeta^2}
    with sigma = 2 (lower half zeta-plane):  zeta_n^2 = ln(4 i sqrt(pi) zeta_n^5) - 2 pi i n."""
    import cmath, math
    z2 = complex(math.log(4 * math.sqrt(math.pi)), -2 * math.pi * n) + 0j
    for _ in range(it):
        z = cmath.sqrt(z2)
        z2 = cmath.log(4j * math.sqrt(math.pi) * z ** 5) - 2j * math.pi * n
    z = cmath.sqrt(z2)
    return z if z.real > 0 else -z


def track(gpath, z, bits=120):
    """Follow a root of eps(zeta; g) = 0 along a list of g values (Newton)."""
    flint.ctx.prec = bits
    for g in gpath:
        for _ in range(60):
            Z, Z1, Z2, _ = Zs(z)
            f = 1 - Z1 / (2 * g)
            fp = -Z2 / (2 * g)
            dz = f / fp
            z = acb((z - dz).mid())
            if dz.abs_upper() < arb(2) ** -90:
                break
    return z


def landau_at(g):
    """Least-damped root at small real g via the Bohm-Gross start."""
    w = (1 + 3 * g).sqrt()
    return track([g], acb(w) / (2 * g).sqrt())


if __name__ == "__main__":
    import sys, math, cmath
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    out = []
    for n in range(1, NMAX + 1):
        seed = asym_seed(n)
        K, g, ok = certify(seed)
        zc = complex(float(K.real.mid()), float(K.imag.mid()))
        gm = complex(float(g.real.mid()), float(g.imag.mid()))
        r, ph = abs(gm), cmath.phase(gm)
        bits = 160 + int(2.2 * abs(zc) ** 2)
        flint.ctx.prec = bits
        g0 = min(0.01, r * 0.9) if r < 0.01 else 0.01
        z_start = landau_at(acb(arb(repr(g0))))
        steps = 300
        path = [acb(arb(repr(g0 + (r - g0) * t / steps))) for t in range(1, steps + 1)]
        path += [acb(arb(repr(r))) * acb(arb(0), arb(repr(ph * t / steps))).exp() for t in range(1, steps + 1)]
        zt = track(path, z_start, bits)
        d = abs(complex(float(zt.real.mid()), float(zt.imag.mid())) - zc)
        row = {"n": n, "zeta*": K.str(25), "g*": g.str(25), "certified": ok,
               "asymptotic_seed": [seed.real, seed.imag], "rel_err_asymptotic": abs(seed - zc) / abs(zc),
               "abs_g*": r, "arg_g*_over_pi": ph / math.pi, "4*pi*n*|g*|": 4 * math.pi * n * r,
               "dist_tracked_landau_root": d}
        out.append(row)
        print(f"n={n:3d} cert={ok} g*={g.str(12)} |g*|={r:.5f} arg/pi={ph/math.pi:.4f} 4pi n|g*|={4*math.pi*n*r:.4f}"
              f" | asym rel.err {row['rel_err_asymptotic']:.1e} | |tracked Landau root - zeta*| {d:.1e}", flush=True)
    json.dump(out, open("results/exceptional_points.json", "w"), indent=1)
