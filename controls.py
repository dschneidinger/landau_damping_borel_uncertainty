"""Controls and convergence for the Borel-sum identity (borel.py).

  planted  : same pipeline on the asymptotic series of Z' itself, D(x) = sum (2n+1)!! x^n,
             whose lateral sum is known exactly: D = Z'(zeta)/(2x), zeta = (2x)^(-1/2).
             Must recover the exponentially small i*sqrt(pi) term from the series alone.
  sheet    : rays theta > 0 must give conj(omega_L^2), NOT the Landau root (check can fail).
  corrupt  : multiply the integers a_n, n >= 40, by (1 + 1e-6): the match must break.
  converge : error vs N at fixed k must fall (not a coincidence of one Pade).
Usage: python controls.py {planted|sheet|corrupt|converge} [N] [dps]
"""
import sys, json
import mpmath as mp
from borel import borel_pade, landau_W
from series import bohm_gross, double_factorials

which = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 120
mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 80
out = {"control": which, "N": N, "dps": mp.mp.dps, "rows": []}

if which == "planted":
    lat = borel_pade([mp.mpf(int(v)) for v in double_factorials(N)])
    Z = lambda z: 1j * mp.sqrt(mp.pi) * mp.exp(-z * z) * mp.erfc(-1j * z)
    for x in ["0.02", "0.05", "0.1", "0.3", "1"]:
        x = mp.mpf(x); z = 1 / mp.sqrt(2 * x)
        exact = -2 * (1 + z * Z(z)) / (2 * x)          # Z'(z)/(2x), Z' = -2(1+zZ)
        for th in ["-1.2", "1.2"]:
            e1 = abs(lat(x, th) - exact); e2 = abs(lat(x, th) - mp.conj(exact))
            out["rows"].append({"x": str(x), "theta": th, "err_vs_Z": float(e1), "err_vs_conj": float(e2),
                                "Im_exact": float(exact.imag)})
            print(f"x={mp.nstr(x,3)} theta={th}: |S-Z'/2x|={mp.nstr(e1,3)}  |S-conj|={mp.nstr(e2,3)}"
                  f"  (Im exact = {mp.nstr(exact.imag,3)})", flush=True)
else:
    a = [mp.mpf(int(v)) for v in bohm_gross(N)]
    if which == "corrupt":
        a = [v * (1 + mp.mpf("1e-6")) if n >= 40 else v for n, v in enumerate(a)]
    if which == "converge":
        Ns = [40, 60, 80, 100, 120, 160, 200]
        ks = ["0.5", "1", "2"]
        Ws, guess = {}, None
        for k in ks:
            Ws[k], guess = landau_W(k, mp.mp.dps, guess)
        for n in Ns:
            lat = borel_pade(a[:n])
            errs = {k: float(abs(lat(mp.mpf(k) ** 2, "-1.2") - Ws[k])) for k in ks}
            out["rows"].append({"N": n, "err": errs})
            print(f"N={n:4d} " + "  ".join(f"k={k}: {e:.3e}" for k, e in errs.items()), flush=True)
    else:
        lat = borel_pade(a)
        guess = None
        for k in ["0.3", "0.5", "1", "2"]:
            W, guess = landau_W(k, mp.mp.dps, guess)
            g = mp.mpf(k) ** 2
            th = "1.2" if which == "sheet" else "-1.2"
            e = abs(lat(g, th) - W); ec = abs(lat(g, th) - mp.conj(W))
            out["rows"].append({"k": k, "theta": th, "err_vs_landau": float(e), "err_vs_conj": float(ec)})
            print(f"k={k} theta={th}: |S-omega_L^2|={mp.nstr(e,3)}  |S-conj(omega_L^2)|={mp.nstr(ec,3)}", flush=True)

json.dump(out, open(f"results/control_{which}.json", "w"), indent=1)
