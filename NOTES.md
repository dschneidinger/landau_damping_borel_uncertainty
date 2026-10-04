# Landau damping is the Borel ambiguity of an integer sequence

**Result.** For the Maxwellian Langmuir wave, the full complex Landau root ω_L(k), damping
included, is the directional Borel sum of the small-k Bohm–Gross series

    ω² = 1 + 3k² + 6k⁴ + 24k⁶ + 180k⁸ + 1872k¹⁰ + 24192k¹² + …        (all coefficients integers)

taken along any Laplace ray strictly below the real axis. Rays above the axis give the complex
conjugate. Concretely:

- Landau damping, the exponentially small Im ω, is exactly half the Borel ambiguity of a series
  of integers.
- Landau's causal contour is the choice of lateral resummation.

Consequences. Items 1–3 are checked against Krawczyk-certified Arb roots of ε(k, ω) = 0. Items 4–5
are exact-rational and interval-certified computations respectively:

1. **Large-order law, with no fitted parameters.** The integers aₙ grow as

       aₙ ≈ (2e^{-3/2}/√π) · Σ_h (−f_h) 2^{n−h} Γ(n + 3/2 − h),

   where f_h are the exact rational coefficients of the damping series. At n = 1000, a
   2870-digit integer is reproduced to **299 digits**, at the predicted 2^{−n}
   two-instanton floor.
2. **All-orders Landau damping formula** with exact rational coefficients:

       γ(k) = −√(π/8) k⁻³ e^{−1/(2k²) − 3/2} · Ψ(k²),
       Ψ(g) = 1 − 3g − (27/2)g² − (177/2)g³ − (7965/8)g⁴ − (557361/40)g⁵ − …

   It also has a closed form, Ψ = ω(ω − k·dω/dk) e^{−(ω² − 1 − 3k²)/(2k²)}, evaluated on the
   Bohm–Gross series. Optimally truncated, it matches the certified γ to 2×10⁻³¹ relative at
   k = 0.08 and 5×10⁻²⁰ at k = 0.1. The textbook formula is off by 2.0% and 3.3% there.
3. **The Stokes multiplier is σ = 1, not 2.** Kowalenko (Phys. Lett. A 337, 408 (2005))
   states that "a proper treatment of the Stokes phenomenon reveals a factor of 2 missing in
   the standard result for Landau damping". In the field-free Maxwellian case, doubling the
   standard rate is off by 100% against the certified roots at every k tested. The resurgent
   reading explains where a 2 can appear:
   - σ = 2 is correct only when paired with the *lateral* sum of the real series, which
     carries an imaginary part of −½ the jump.
   - Pairing σ = 2 with the real (truncated or median) series double counts.
   - I read only the abstract, not the full text. Whether its claim is meant to cover the
     field-free limit is unconfirmed.

4. **The whole resurgent tower closes on one Stokes jump** (`transseries.py`).
   - Every sector of the transseries is computed exactly as a rational power series Ψ_m(g).
     The method is Lagrange inversion in the rescaled variable Z = (W − W₀)/g.
   - The single Stokes jump of Z (σ → σ + 2) fixes all bridge relations,
     Disc u_m = Σ_j C(m+j, m) 2^j u_{m+j}.
   - The coefficients of each sector m = 0…4 are therefore predicted by sector m+1, with Stokes
     multipliers 2, 4, 6, 8, 10 and no fitting.
   - With 800 terms per sector (`results/transseries_L800.txt`), every level matches to between
     3×10⁻¹¹⁸ and 10⁻¹¹⁶ at n = 399, near the 2^{−n} floor of 7.7×10⁻¹²¹. Shifting any multiplier
     by one gives errors of 0.17 to 0.5.
5. **The Landau root has infinitely many exceptional points accumulating at k = 0** (`ep.py`).
   - At an exceptional point (EP) the Langmuir root merges with another kinetic root:
     ε = ∂ε/∂ζ = 0, equivalently Z″(ζ*) = 0 with g* = Z′(ζ*)/2.
   - The first 30 EPs are each certified by a Krawczyk test.
   - Continued from real k (along the real axis, then the arc |g| = |g*ₙ|), the physical Landau
     root reaches every one of them, to the ~10⁻⁸ floor of Newton's method at a double root.
   - They accumulate at k² = 0 along the imaginary axis: 4πn|g*ₙ| → 1 and arg g*ₙ → π/2, so
     g*ₙ ≈ i/(4πn). The first is at k² = 0.036389 + 0.035523i.
   - Balancing the two parts of Z″ with the off-line Stokes multiplier σ = 2 gives
     ζₙ² = ln(4i√π ζₙ⁵) − 2πi n. This predicts the certified EPs with relative error falling
     from 2×10⁻² (n = 1) to 5×10⁻⁵ (n = 30).
   - So the Landau root is not analytic in any neighbourhood of k = 0. All the EPs lie in the
     upper half of the k² plane, which is consistent with only the lower lateral sum representing
     the physical root. I have not proved that link.

Units throughout: ω_p = λ_D = v_th = 1, g = k², W = ω², ζ = ω/(√2 k).

## Derivation

The dispersion relation is ε = 1 − Z′(ζ)/(2g) = 0. For large ζ,

    Z′(ζ) = Σₙ (2n+1)!!/(2ⁿ ζ^{2n+2}) − 2iσ√π ζ e^{−ζ²}.

**Perturbative sector.** Put x = g/W = 1/(2ζ²). The asymptotic part reads
P(x) := Σ (2n+1)!! x^{n+1} = g. So x(g) is the series reversion of P, and W₀ = g/x(g). Its
coefficients are integers because P has integer coefficients and leading coefficient 1
(`series.py`).

**One-instanton sector.** Write G(W, g) = W⁻¹ D(g/W). The inverse-function rule gives
∂_W G = −1/(W² x′(g)), which simplifies to 1/∂_W G = −(W − gW′). Then

    δW = −iσ√π (W/2)^{1/2} g^{−3/2} e^{−W/(2g)} (W − gW′),
    δω = −iσ√(π/8) k⁻³ ω(ω − k v_g) e^{−ω²/(2k²)},

with v_g = dω/dk on the Bohm–Gross series. This is the classic weak-damping formula. The new
point is that it is the exact one-instanton sector of a resurgent transseries, so all its
rational coefficients are meaningful.

**Stokes structure.**
- Z = i√π e^{−ζ²} erfc(−iζ) is entire, and its asymptotic series has a single Stokes line on
  real ζ, where the jump is 2i√π e^{−ζ²}.
- Analytically continuing the upper-half-plane lateral Borel sum gives exactly Z.
- Lateral Borel summation commutes with the algebra and the implicit-function step. The
  1-summable series form a differential algebra closed under composition and implicit
  functions; see Balser, *Formal Power Series and Linear Systems of Meromorphic ODEs*.
- So the damped root is the lateral Borel sum of W₀, and Re ω_L² and Im ω_L² are the half-sum
  and half-difference of the two lateral sums.
- Large-order behaviour then follows from a Cauchy/dispersion integral over the discontinuity
  2τW₁. This gives law (1).

## Evidence (all files in `results/`)

### Large order (`large_order.py`)

The integers are predicted from Ψ alone: s = −1 orientation, optimal truncation in h.

| n | digits of aₙ | leading-term rel. err | resurgent rel. err | 2^{−n} |
|---|---|---|---|---|
| 125 | 248 | 6.2e-3 | 2.7e-35 | 2.4e-38 |
| 250 | 569 | 3.1e-3 | 6.8e-74 | 5.5e-76 |
| 500 | 1286 | 1.5e-3 | 7.4e-149 | 3.1e-151 |
| 1000 | 2870 | 7.5e-4 | **4.5e-299** | 9.3e-302 |

### Borel sum vs certified root (`borel.py`, `controls.py converge`)

The table gives |S_θ[W₀](k²) − ω_L²| for the θ = −1.2 ray. Every ω_L carries a Krawczyk
uniqueness certificate.

Convergence run at dps = 110 (`results/control_converge.txt`):

| k | ω_L (certified, 20 digits) | N=40 | N=60 | N=80 | N=100 | N=120 | N=160 | N=200 |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 1.4156618886045364298 − 0.15335946690960483409 i | 2.6e-8 | 5.6e-11 | 7.7e-13 | 2.3e-15 | 2.1e-17 | 2.2e-19 | 5.7e-19 |
| 1 | 2.0459048656906245943 − 0.85133045869205615163 i | 1.0e-5 | 7.2e-8 | 6.3e-9 | 4.4e-11 | 9.1e-13 | 2.5e-14 | 4.9e-14 |
| 2 | 3.1891361929982971943 − 2.8272002686707789326 i | 1.8e-3 | 2.7e-5 | 5.9e-6 | 1.2e-7 | 5.6e-9 | 3.7e-10 | 5.7e-10 |

Single-N runs (`results/borel_largek.txt` and the first `borel.py` run):

| k | ω_L (certified, 20 digits) | N | dps | error |
|---|---|---|---|---|
| 0.2 | 1.0639843406877026103 − 5.5107414768628414033e-5 i | 120 | 80 | 3.0e-27 |
| 3 | 4.2449483969497490581 − 5.1730183051984860268 i | 160 | 100 | 3.6e-8 (\|ω²\| ≈ 45) |
| 5 | 6.2184893169398914267 − 10.458500590204565806 i | 160 | 100 | 6.4e-6 (\|ω²\| ≈ 148) |
| 8 | ω² = −289.85887832747460511 − 346.07795255120472009 i | 160 | 100 | 4.0e-4 (\|ω²\| ≈ 451) |

From N = 40 to N = 160 the error falls steadily at every k. From N = 160 to N = 200 it levels
off and rises slightly at all three k. The likely cause is the Padé approximant running into the
110-digit working precision, but I have not tested that.

Near-axis rays (θ = −0.25) fail for k ≳ 0.27,
while rotated rays work. This is consistent with complex-k² exceptional points of the root,
the nearest at k² ≈ 0.0364 + 0.0355i, where a lateral sum's half-plane of convergence cannot
reach. Those exceptional points are located but not certified.

### All-orders damping (`damping.py`)

| k | certified γ | textbook rel. err | Ψ (optimal) rel. err | σ = 2 rel. err |
|---|---|---|---|---|
| 0.08 | −3.15054949672458e-32 | 2.0e-2 | 2.2e-31 | 1.0 |
| 0.10 | −2.61207782362829e-20 | 3.3e-2 | 5.1e-20 | 1.0 |
| 0.12 | −6.42360194962512e-14 | 4.9e-2 | 2.8e-13 | 1.0 |
| 0.15 | −8.5531013073366e-9 | 8.2e-2 | 2.1e-8 | 1.0 |
| 0.20 | −5.51074147686284e-5 | 1.8e-1 | 4.8e-5 | 1.0 |

### Controls (`controls.py`)

- **Planted truth.** The identical pipeline, applied to the bare asymptotic series of Z′,
  recovers the exact Z′(ζ)/(2x) to 1e-79 to 3e-80, including its exponentially small
  imaginary part, at x = 0.02, 0.05 and 0.1.
- **Wrong sheet.** Rays with θ > 0 return conj(ω_L²) and fail against ω_L². For example, at
  k = 0.3 the error is 0.059 against ω_L² and 3e-20 against the conjugate.
- **Corruption.** Multiplying aₙ (n ≥ 40) by 1 + 10⁻⁶ destroys the match (error ≥ 10¹⁷).

## What is and is not established

- **Rigorous.** The integers aₙ and the rational coefficients of every sector Ψ_m are exact. The
  roots ω_L(k) and the 30 exceptional points carry interval-Newton certificates.
- **Also numerical.** Two steps in the extensions are numerical rather than certified:
  - That the Landau root reaches each EP comes from Newton tracking along a path at 160+ bits.
    Near a double root Newton resolves the root only to about the square root of its tolerance,
    which explains the ~10⁻⁸ match. The first 6 EPs were also
    checked with endpoints rotated by ±0.02 in phase, which moved the root off by ~0.1.
  - The tower's large-order agreements are high-precision evaluations of exact quantities, not
    interval bounds.
- **Numerical, not certified.** The Borel–Padé sums are uncertified numerics. Their agreement
  with certified roots improves steadily with N up to N = 160, reaching 10⁻⁶ to 10⁻²⁷ depending
  on k, then levels off.
- **Theorem status.** The all-k identity is a theorem modulo the standard summability-closure
  results cited above. I have written a proof sketch, not a full proof.
- **Novelty.** The check was limited:
  - The OEIS has no match for 1, 3, 6, 24, 180, 1872, ….
  - Web searches found the Kowalenko Stokes/Landau paper but no resurgent treatment of the
    Landau root, no large-order law, and no all-orders Ψ.
  - Exceptional points of Vlasov dispersion relations are a known general concept; for example, a
    2026 preprint (arXiv:2608.19624) studies them in magnetized plasmas. I found no prior
    account of the infinite EP sequence of the Maxwellian Langmuir root, its i/(4πn)
    accumulation law, or the exact resurgence tower.
  - A proper literature audit is still needed before any priority claim.
- **Scope.** This is a structural result about linear Landau damping of the Maxwellian
  Langmuir wave. It changes no practical damping rate, since the textbook formula was already
  adequate where the damping is small.

## Reproduce

    ../.venv/bin/python series.py              # integers and Phi1
    ../.venv/bin/python large_order.py 1000    # ~20 s
    ../.venv/bin/python damping.py
    ../.venv/bin/python borel.py 120 80 "-0.25,-0.7,-1.2,-1.5"
    ../.venv/bin/python controls.py {planted|sheet|corrupt|converge}
    ../.venv/bin/python transseries.py 800 5   # resurgence tower, m = 0..4 (~3 min)
    ../.venv/bin/python ep.py 30               # certified exceptional points
