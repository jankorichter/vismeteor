## Release summary

This is a major release (3.0.0) with breaking changes:

* All exported functions, parameters, and `load_vmdb_*()` columns were
  renamed to snake_case to align with tidyverse / r-lib conventions and
  with the imo-vmdb API (which itself uses snake_case). A migration sed
  snippet for downstream scripts is included in NEWS.md.
* `vmperception()` parameter renamed from `m` to `dm` to reflect that it
  is the difference (limiting magnitude − meteor magnitude); positional
  calls remain compatible.
* `perception_fun` now defaults to `vmperception` instead of `NULL` with
  an internal fallback.
* New `select_knots()` performs forward/backward stepwise selection of
  spline knots with a user-supplied scoring function (see
  `vignette("select_knots")`).

## Test environments

* local macOS (Darwin 26.4), R 4.6.0 — R CMD check --as-cran
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
