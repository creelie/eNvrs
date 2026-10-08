**Blowup above the exponent one third in viscous Katz–Pavlović dyadic models with small shell ratios**, by Deep Bhattacharjee.

For the viscous Katz–Pavlović dyadic model with shell ratio Λ > 2^{3/2}, global regularity holds exactly for α ≥ 1/3. This release proves that the threshold is not universal: for Λ = 6/5, 13/10, 3/2, 17/10 and 9/5, explicit nonnegative finitely supported data lose regularity in finite time for every α up to 0.377, 0.378, 0.365, 0.348 and 0.339 respectively, for every viscosity. The blowup is carried by a self-similar front whose amplitude grows by a factor κ > 1 per shell.

The proof is computer-assisted: the reductions are proved by hand, and one renormalized step of the front is verified in interval arithmetic with CAPD 6.1.0. It is rigorous provided CAPD and the verifier are correct. The result concerns the dyadic model, not the Navier–Stokes equations; the threshold for general Λ is a conjecture supported by computation.

Files:

- `dyadic-blowup.pdf`: the paper;
- `dyadic-blowup-tex.zip`: LaTeX source with the figures as PNG and their TikZ sources;
- `dyadic-blowup-arxiv.tar.gz`: LaTeX source with the figures as PDF.

The repository at this tag also holds the verifier, its output for the five cases, the independent Python, C and Julia checks, and the Lean 4 proofs of eight elementary lemmas.

Version 1.1.0 adds the DOI of the Zenodo archive (concept DOI [10.5281/zenodo.23247974](https://doi.org/10.5281/zenodo.23247974)) to the paper, the README and the citation file. The mathematics is unchanged from version 1.0.0 ([10.5281/zenodo.23247975](https://doi.org/10.5281/zenodo.23247975)).
