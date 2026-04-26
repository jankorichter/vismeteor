## Release summary

This is a minor release (2.1.0) with breaking changes.
`load_vmdb_rates()` and `load_vmdb_magnitudes()` now connect to an
[imo-vmdb](https://pypi.org/project/imo-vmdb/) REST API instead of a
direct database connection. Minimum R version raised to 4.1.0.

## Test environments

* local macOS, R 4.5.2
* win-builder (R-devel)
* GitHub Actions (ubuntu-latest, windows-latest, macOS-latest), R-release + R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no reverse dependencies for this package.
