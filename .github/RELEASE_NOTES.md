**Blowup above the exponent one third in viscous Katz–Pavlović dyadic models with small shell ratios**, by Deep Bhattacharjee.

For the viscous Katz–Pavlović dyadic model with shell ratio Λ > 2^{3/2}, global regularity holds exactly for α ≥ 1/3. This release proves that the threshold is not universal: for Λ = 6/5, 13/10, 3/2, 17/10 and 9/5, explicit nonnegative finitely supported data lose regularity in finite time for every α up to 0.377, 0.378, 0.365, 0.348 and 0.339 respectively, for every viscosity. The blowup is carried by a self-similar front whose amplitude grows by a factor κ > 1 per shell.

The proof is computer-assisted: the reductions are proved by hand, and one renormalized step of the front is verified in interval arithmetic with CAPD 6.1.0. It is rigorous provided CAPD and the verifier are correct. The result concerns the dyadic model, not the Navier–Stokes equations; the threshold for general Λ is a conjecture supported by computation.

Files:

- `dyadic-blowup.pdf`: the paper;
- `dyadic-blowup-tex.zip`: LaTeX source with the figures as PNG and their TikZ sources;
- `dyadic-blowup-arxiv.tar.gz`: LaTeX source with the figures as PDF;
- `dyadic-blowup-physica-d.pdf`: the same paper in Elsevier's elsarticle class, formatted for Physica D: Nonlinear Phenomena;
- `dyadic-blowup-physica-d-source.zip`: its LaTeX source, one self-contained `main.tex` with the figures as PNG.

The repository at this tag also holds the verifier, its output for the five cases, the independent Python, C and Julia checks, and the Lean 4 proofs of eight elementary lemmas.

Version 1.2.0 adds the Physica D version of the paper. It is assembled from the same source by `scripts/build_physica_d.py`, with Elsevier's front matter, declarations and numbered references; the text of the theorems and proofs is identical. One wide display in Proposition 3.1 is split over two lines. The mathematics is unchanged from versions 1.0.0 ([10.5281/zenodo.23247975](https://doi.org/10.5281/zenodo.23247975)) and 1.1.0 ([10.5281/zenodo.23248055](https://doi.org/10.5281/zenodo.23248055)); the concept DOI [10.5281/zenodo.23247974](https://doi.org/10.5281/zenodo.23247974) covers all versions.
