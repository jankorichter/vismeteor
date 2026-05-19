## Release summary

This is a patch release (3.0.1) with one bug fix and no API changes:

* `load_vmdb_rates()` and `load_vmdb_magnitudes()` now correctly request
  the per-magnitude-class frequencies from the imo-vmdb REST API.  The
  previous release sent `include=magnitudes`, which on `/magnitudes`
  returned `HTTP 400` and on `/rates` produced a response shape the
  internal parser could not handle.  The functions now send
  `include=magnitude_details`, restoring the documented behaviour.
  This requires imo-vmdb >= 2.0.0 on the server side; older servers reject the new include combination with `HTTP 400`.

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
