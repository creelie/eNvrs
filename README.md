# eNvrs: blowup above the exponent one third in dyadic models

Code, data and paper for

> Deep Bhattacharjee, *Blowup above the exponent one third in viscous
> Katz–Pavlović dyadic models with small shell ratios*.

Archived on Zenodo: [doi:10.5281/zenodo.23247974](https://doi.org/10.5281/zenodo.23247974)
(all versions; v1.0.0 is [doi:10.5281/zenodo.23247975](https://doi.org/10.5281/zenodo.23247975)).

## What is proved

The viscous Katz–Pavlović dyadic model is the infinite system

    du_n/dt + ν Λ^{2αn} u_n = Λ^{n-1} u_{n-1}^2 − Λ^n u_n u_{n+1},   n ≥ 0,   u_{-1} = 0.

For Λ > 2^{3/2}, global regularity holds exactly for α ≥ 1/3 (Cheskidov 2008
for α < 1/3; Zhang 2026 for α ≥ 1/3). The paper proves that this threshold is
not universal: for

| Λ     | blowup for every α ≤ | growth factor κ per shell        |
|-------|----------------------|----------------------------------|
| 6/5   | 0.377                | [1.01628489, 1.01628532]         |
| 13/10 | 0.378                | [1.02401568, 1.02401602]         |
| 3/2   | 0.365                | [1.02642827, 1.02642851]         |
| 17/10 | 0.348                | [1.01599737, 1.01599757]         |
| 9/5   | 0.339                | [1.00747733, 1.00747753]         |

explicit nonnegative finitely supported data lose regularity in finite time,
for every viscosity ν > 0. The blowup is carried by a self-similar front
whose amplitude grows by the factor κ > 1 per shell.

The proof is computer-assisted. The reduction to one renormalized step of the
front, the control of the infinitely many shells outside a finite window, the
truncation error and the iteration are proved by hand (Sections 2–4 of the
paper). What remains is a finite list of inequalities about one system of
ordinary differential equations, which `cap/verify_step.cpp` checks in
interval arithmetic with the [CAPD](https://github.com/CAPDGroup/CAPD)
library. The proof is rigorous provided CAPD and this program are correct.

## What is not claimed

This is a theorem about a dyadic model, not about the Navier–Stokes
equations, and it does not touch the Clay Millennium problem. It is proved at
the five ratios above, not for all Λ < Λ_c ≈ 1.8754; the conjectured
threshold α\*(Λ) = 1/3 + log κ(Λ)/(2 log Λ) (Section 6) is supported by
computation only. For Λ = 13/10 the range 0.378 < α < 1/2 remains open.

## Layout

| Path | Contents |
|------|----------|
| `paper/` | `main.tex` (amsart) and the TikZ figures with their data |
| `cap/verify_step.cpp`, `cap/MyHOE.h` | the interval verifier (CAPD 6.1.0) |
| `cap/cases/` | the profiles `w*` and the box radii of the five cases |
| `cap/logs/` | the verifier's output for each case, with all enclosures |
| `cap/run_all.sh` | builds the verifier and runs the five cases |
| `cap/make_case.py` | how the profiles and radii were produced (not part of the proof) |
| `crosscheck/recheck.py` | 60-digit interval re-check of the scalar inequalities from the logs |
| `crosscheck/front_steps.c` | independent long-double iteration of the full model |
| `crosscheck/kappa.jl` | independent Julia computation of κ(Λ) |
| `numerics/` | fixed-point computations, the κ table and the simulations of Section 6 |
| `lean/` | Lean 4 proofs (Mathlib) of eight elementary lemmas, with an axiom audit |
| `scripts/` | `build_capd.sh`, `fast_checks.sh`, `build_paper.sh` |

## Reproducing

    scripts/build_capd.sh capd-install            # CAPD at commit 2f06098, about 10 minutes
    cap/run_all.sh capd-install/bin/capd-config   # the five cases, about 20 seconds
    scripts/fast_checks.sh                        # interval re-check, C and Julia checks
    (cd lean && lake exe cache get && lake build && lake env lean CheckDyadic.lean)
    scripts/build_paper.sh                        # dist/: PDF, tex.zip, arXiv tarball

Each case prints `PASS` only if every inequality (V1)–(V9) of Section 5 holds.
The CI workflow `.github/workflows/verify.yml` runs all of this on every pull
request.

## Licence and citation

MIT licence, © 2026 Deep Bhattacharjee. See `CITATION.cff` for how to cite the
code and the paper.
