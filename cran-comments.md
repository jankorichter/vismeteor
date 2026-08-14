## Release summary

This is a minor release (3.1.0) adding regression interfaces for the two
magnitude models:

* `vmgeom_glm()` and `vmideal_glm()` fit the geometric model and the ideal
  distribution as generalized linear models, so the population index `r` and
  the location parameter `psi` can depend on covariates.
* `vmgeom_vst_lm()` fits the variance-stabilized magnitudes of the geometric
  model with `lm()`, with `predict()` returning `r`, `1/r` or `log(r)`.

The release also recalibrates both variance-stabilizing transformations. The
transformed magnitudes therefore live on a different scale than in 3.0.1, and
`dmideal()`, `pmideal()` and `qmideal()` now reject missing values instead of
returning `NA`/`NaN`. Both changes are documented under "Breaking changes" in
NEWS.md.

## Test environments

* local macOS (Darwin 25.6), R 4.6.1 — R CMD check --as-cran
* GitHub Actions: ubuntu-latest, windows-latest, macOS-latest;
  R-release and R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes on GitHub Actions
(ubuntu-latest, windows-latest, macOS-latest; R-release and R-devel).

Local `R CMD check --as-cran` on macOS produces 1 NOTE that is unrelated to
the package: HTML manual validation is skipped because the system `tidy` is
older than the version required by `R CMD check`, and math rendering is
skipped because `V8` is not installed. Both are local-toolchain issues and
do not occur on CRAN's check servers.

## Reverse dependencies

There are currently no reverse dependencies for this package.
