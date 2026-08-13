# vismeteor unreleased

## New features

- `vmgeom_glm()` fits the geometric model of visual meteor magnitudes as a
  generalized linear model, allowing the population index `r` to be modelled
  as a function of covariates.  The response is given as a
  two-column matrix `cbind(magn, lim_magn)`, so that the limiting magnitude
  is subsetted together with the remaining data.  The link function is
  `logit(1/r)`, which enforces `r > 1` for every finite linear predictor.
  Magnitude counts can be supplied as `weights`, which is equivalent to
  repeating each row accordingly; the dispersion is fixed at `1`, so the
  aggregated and the replicated form also agree in their standard errors.
- `predict()` on such a fit returns the population index (`type = "r"`, the
  default), `1/r`, `log(r)`, or the linear predictor.  With `se.fit = TRUE` it
  also returns standard errors and confidence limits, which are formed on the
  link scale and therefore never fall below `1`.
- `vignette("vmgeom")` gained a section demonstrating the new model.
- `vmideal_glm()` fits the ideal distribution of visual meteor magnitudes as a
  generalized linear model, allowing the location parameter `psi` to be
  modelled as a function of covariates.  The response is given as a two-column
  matrix `cbind(magn, lim_magn)` and magnitude counts can be supplied as
  `weights`, both as for `vmgeom_glm()`.  The link is the identity, since `psi`
  is unrestricted on the real line.  Unlike the geometric model, the mean is
  not a sufficient statistic for `psi`, so the fit is a quasi-likelihood
  estimate: it is consistent, but not fully efficient.  How closely it tracks
  the maximum likelihood estimate depends on the sample size and on how far
  `psi` lies above the limiting magnitude — from a few hundredths for a few
  thousand meteors at a moderate `psi` to a whole magnitude for small samples
  near the limit of identifiability.  See `?vmideal_glm` for details.
- `predict()` on such a fit returns `psi` (`type = "psi"`, the default) or the
  linear predictor.  With `se.fit = TRUE` it also returns standard errors and
  symmetric confidence limits.  Where the mean magnitude has reached the value
  the ideal distribution converges to, no data can distinguish a larger `psi`
  from a smaller one and `type = "psi"` returns `Inf`.  This is a result rather
  than an error — the magnitudes are then geometric with `r = 10^0.4` — and it
  also covers data too faint for the model at any `psi`.  The linear predictor
  stays finite, so `type = "link"`, `summary()` and `anova()` remain usable.
  Close to that point the mean barely moves with `psi` and the iteration can
  run out of steps; such a fit is reported with a warning rather than returned
  as if it had converged.
- `vignette("vmideal")` gained a section demonstrating the new model.

## Changes

- `dvmideal()`, `pvmideal()`, `qvmideal()` and `rvmideal()` accept a `psi` of
  `Inf`, which denotes the geometric limit the ideal distribution converges to,
  with a population index of `r = 10^0.4`.  The `Inf` that
  `predict(type = "psi")` returns can therefore be passed on unchanged instead
  of having to be intercepted by the caller.  The result is the geometric one
  exactly, not an approximation of it, and a large finite `psi` was already
  evaluated that way.  A `psi` of `-Inf` remains an error, as the distribution
  has no limit there.
- `dvmideal()`, `pvmideal()` and `qvmideal()` keep a `perception_fun` given to
  them when the location parameter is large enough for the geometric model to
  answer.  Previously they fell back to `vmperception()` in that case.
- `vmideal_vst_to_psi()` returns `Inf` for a `tm` below `0.001`, where it
  previously returned `NA`.  A vanishing mean of the transformed magnitudes is
  what an unbounded `psi` produces, and the estimate can now say so instead of
  reporting the parameter as unknown.  The result may be passed straight to
  `dvmideal()`.  The upper end above `8.22` stays `NA`, as do the derivatives
  used for the delta method wherever `psi` is not finite.
- `vmgeom_vst_from_magn()` and `vmgeom_vst_to_r()` use a new algorithm.  The
  transformation is now derived from the rate-based statistic
  `f(lm - m - 1) / f(lm - m)`, whose expectation is `1/r`, and stabilises its
  variance with a power transformation.  This replaces the interpolated table
  of fitted parameters.  Both functions keep their arguments and their usage
  unchanged; `vignette("vmgeom")` explains the relation to the rate-based
  method.
- The scale of the transformed magnitudes changed accordingly: `1.4 <= r <= 4`
  now corresponds to `tm` between about `1.59` and `3.53`, previously `3.96` and
  `5.74` over the narrower range `1.4 <= r <= 3.5`.  Transformed values stored
  with earlier versions are not comparable.
- `vmgeom_vst_to_r()` maps the mean of the transformed magnitudes back onto
  `r`.  Its constants are obtained by regressing `log(E[g])` on `log(E[t])`
  under the model.  Inverting the power transformation algebraically would be
  correct for a single meteor only: the estimator averages transformed
  magnitudes, and the mean of a power is not the power of the mean.
- `vmgeom_vst_to_r()` no longer restricts its argument to a calibrated window.
  It extrapolates monotonically instead, and returns `NA` only for values that
  cannot occur, that is negative ones and those above `4.52`.  This is aimed at
  predictive models, where sparse data can yield an implausible `r`: such a
  value is now returned inflated but ordered, instead of dropping out as `NA`.
  The recovered `r` is strictly monotone in `tm` throughout, so the estimator
  never saturates or folds back.
- The back-transformation carries a small systematic deviation that does not
  shrink as the sample grows.  Over `1.7 <= r <= 3.3`, the range met in
  practice, it stays below `1.5%` and is negligible against the random error;
  it shrinks towards smaller `r` and grows towards larger ones.  It stems from
  a slight dependence of the transformed mean on the limiting magnitude; where
  it matters, use the rate-based estimator or `vmgeom_glm()`.
- The variance stabilization, unlike the above, remains two-sided and is what
  the range `1.4 <= r <= 4` now describes: there the variance stays within
  `0.19` of `1.0`, compared with `0.018` over `1.4 <= r <= 3.5` before.  It
  bounds the use of the transformed magnitudes as a response in linear models,
  not the reading of a single `r` off their mean.

## Performance

- Fitting a model whose population index depends on covariates is roughly ten
  times faster.  The perception probabilities depend on the fractional part of
  the limiting magnitude alone, and are now evaluated once per distinct value
  instead of once per observation.  The results are unchanged.

# vismeteor 3.0.1

## Bug fixes

- `load_vmdb_rates()` and `load_vmdb_magnitudes()` now correctly request
  the per-magnitude-class frequencies from the imo-vmdb REST API.
  Previously `with_magnitudes = TRUE` sent `include=magnitudes`, which on
  `/magnitudes` returned `HTTP 400` and on `/rates` produced a response
  shape that `.parse_magnitudes()` could not handle.  The functions now
  send `include=magnitude_details`, restoring the documented behaviour.
- The `period_start` and `period_end` columns of the `observations`
  data.frame returned by `load_vmdb_rates()` and `load_vmdb_magnitudes()`
  are now `POSIXct` (UTC), matching the behaviour of vismeteor 2.x.
  They were inadvertently returned as `character` after the switch to
  `httr2` in 3.0.0.
- The `period` argument of `load_vmdb_rates()` and
  `load_vmdb_magnitudes()` is now serialised to the strict ISO 8601
  form (`YYYY-MM-DDTHH:MM:SS`, UTC) that imo-vmdb 2.0 expects.  `Date`,
  `POSIXct`, and character (full datetime or date-only) inputs are all
  accepted; date-only inputs are expanded to midnight (lower) and
  23:59:59 (upper) of the given day.

## Compatibility

- Requires **imo-vmdb ≥ 2.0.0** on the server side.  2.0 introduced
  the `include=magnitude_details` parameter and the strict ISO 8601
  wire format (`T` separator, no timezone marker) for `period_start` /
  `period_end`; both are required by this release.

# vismeteor 3.0.0

## New features

- `select_knots()` performs forward/backward stepwise selection of spline
  knots from a candidate set, scoring each fit with a user-supplied function
  (e.g. AIC/BIC). Backward selection supports a "bulk removal" mode, and
  scoring can run in parallel via the `parallel` package. See
  `?select_knots` and `vignette("select_knots")`.

## Breaking changes

### Naming overhaul (snake_case)

All public identifiers were migrated from dotted / camelCase to `snake_case` to
match the tidyverse / r-lib convention and to remove visual confusion with S3
dispatch (`print.foo`). This affects exported functions, function parameters,
and `data.frame` columns returned by `load_vmdb_*()` and the example datasets
`PER_2015_rates` / `PER_2015_magn` (which were regenerated).

The d/p/q/r distribution prefixes (`dvmgeom`, `pvmideal`, ...) are unchanged,
and `lower.tail` is preserved to match base R conventions.

Renamed exported functions:

| 2.1.0                | 3.0.0                    |
|----------------------|--------------------------|
| `freq.quantile`      | `freq_quantile`          |
| `vmgeomVstFromMagn`  | `vmgeom_vst_from_magn`   |
| `vmgeomVstToR`       | `vmgeom_vst_to_r`        |
| `vmidealVstFromMagn` | `vmideal_vst_from_magn`  |
| `vmidealVstToPsi`    | `vmideal_vst_to_psi`     |

Renamed parameters: `lim.magn`, `magn.id`, `rate.id`, `session.id`,
`perception.fun`, `sun.alt.max`, `moon.alt.max`, `deriv.degree`, `withSessions`,
`withMagnitudes` → their snake_case forms.

### Closer alignment with the imo-vmdb API

The imo-vmdb API already uses snake_case, so column remapping is reduced to a
few semantic renames: `id` → `rate_id` / `magn_id` / `session_id` (so foreign
keys in the same data frame stay unambiguous) and `mean` → `magn_mean` (to
avoid shadowing base R `mean()`). Three columns that previously carried
R-specific names now pass through unchanged from the API:

| 2.1.0           | 3.0.0       |
|-----------------|-------------|
| `shower.code`   | `shower`    |
| `radiant.alt`   | `rad_alt`   |
| `radiant.az`    | `rad_az`    |

### Other breaking changes

- `perception_fun` now defaults to `vmperception` instead of `NULL` with an
  internal fallback. Callers that explicitly passed `perception.fun = NULL`
  must drop the argument or supply a function.
- `vmperception(m)` is now `vmperception(dm)` — the parameter is the
  difference between the limiting magnitude and the meteor magnitude, and
  the new name reflects that. Positional calls remain compatible.

### Migration

A `sed` sweep of your scripts is usually sufficient:

```sh
sed -i '' -E '
  s/\bfreq\.quantile\b/freq_quantile/g;
  s/\bvmgeomVstFromMagn\b/vmgeom_vst_from_magn/g;
  s/\bvmgeomVstToR\b/vmgeom_vst_to_r/g;
  s/\bvmidealVstFromMagn\b/vmideal_vst_from_magn/g;
  s/\bvmidealVstToPsi\b/vmideal_vst_to_psi/g;
  s/\blim\.magn\b/lim_magn/g;
  s/\bmagn\.id\b/magn_id/g;
  s/\brate\.id\b/rate_id/g;
  s/\bsession\.id\b/session_id/g;
  s/\bperception\.fun\b/perception_fun/g;
  s/\bderiv\.degree\b/deriv_degree/g;
  s/\bwithSessions\b/with_sessions/g;
  s/\bwithMagnitudes\b/with_magnitudes/g;
  s/\bshower\.code\b/shower/g;
  s/\bperiod\.start\b/period_start/g;
  s/\bperiod\.end\b/period_end/g;
  s/\bsl\.start\b/sl_start/g;
  s/\bsl\.end\b/sl_end/g;
  s/\bt\.eff\b/t_eff/g;
  s/\btime\.sidereal\b/sidereal_time/g;
  s/\bsun\.alt\b/sun_alt/g;
  s/\bsun\.az\b/sun_az/g;
  s/\bmoon\.alt\b/moon_alt/g;
  s/\bmoon\.az\b/moon_az/g;
  s/\bmoon\.illum\b/moon_illum/g;
  s/\bfield\.alt\b/field_alt/g;
  s/\bfield\.az\b/field_az/g;
  s/\bradiant\.alt\b/rad_alt/g;
  s/\bradiant\.az\b/rad_az/g;
  s/\bmagn\.mean\b/magn_mean/g;
  s/\blocation\.name\b/location_name/g;
  s/\bobserver\.id\b/observer_id/g;
  s/\bobserver\.name\b/observer_name/g;
' your-script.R
```

# vismeteor 2.1.0

## Breaking changes

- `load_vmdb_rates()` and `load_vmdb_magnitudes()` now connect to an
  [imo-vmdb](https://pypi.org/project/imo-vmdb/) REST API instead of a
  direct database connection. 
- Multi-range filters for `period`, `sl`, and `lim_magn` (previously a
  matrix with multiple rows that were OR-joined) are no longer supported.
  Each filter is collapsed to a single bounding min/max. Use multiple calls
  combined with `rbind()` if disjoint ranges are needed.
- Minimum R version raised to **4.1.0** (was 3.5.0) due to use of the native pipe operator `|>` and `\()`.

# vismeteor 2.0.2

## Changes

- Clarified documentation and vignettes for `vmgeom` and `vmideal` models, including improved mathematical notation and explanations.
- Hidden internal parameters from exported function documentation to reduce user confusion.

# vismeteor 2.0.1

## Changes

- Updated README guidance, vignettes, and roxygen references to clarify data sources, models, and to correct documentation errors.
- Improved the performance of `vmgeom_vst_from_magn()` and `vmideal_vst_from_magn()` by tightening the interpolation steps used during the variance-stabilising transforms.
- Streamlined the `vmtable()` rounding routine to reduce allocations.
- Added a GitHub Actions workflow to run package checks automatically.
- Documented the derivation scripts in `inst/derivation/` to improve reproducibility.
- Added `load_vmdb()` regression tests.

# vismeteor 2.0.0

## Highlights
This release introduces variance-stabilizing transformations for the ideal distribution of visual meteor magnitudes (`vmideal_vst_from_magn()`) as well as for visual meteor magnitudes under a geometric distribution (`vmgeom_vst_from_magn()`).

## Other changes
The function `vmperception()` now better matches the perception probabilities of *Koschack & Rendtel (1990b)*.  
The argument `deriv_degree` has been removed, as it was only intended for internal testing and had no practical relevance for regular use.

Laplace-transformed perception probabilities have been replaced by variance-stabilizing transformations, which also means that the function `vmperception.l()` has been removed.

**Note:** This release includes breaking changes and is **not fully backward compatible** due to the removal of parameters and functions.

# vismeteor 1.8.5
- Initial CRAN release
