#' @name vmgeom_vst
#' @aliases vmgeom_vst_from_magn
#' @aliases vmgeom_vst_to_r
#' @title Variance-Stabilizing Transformation for Geometric Visual Meteor Magnitudes
#'
#' @description
#' Applies a variance-stabilizing transformation to visual meteor magnitudes
#' under the geometric model.
#'
#' @param m integer; meteor magnitude.
#' @param lm numeric; limiting magnitude.
#' @param tm numeric; transformed magnitude.
#' @param log logical; if `TRUE`, the logarithm of the population index `r` is returned.
#' @param deriv_degree integer; the order of the derivative at `tm` to return
#'   instead of `r` or `log(r)`. Must be `0`, `1`, or `2`.
#'
#' @details
#' `vmgeom_vst_from_magn` maps visual meteor magnitudes onto a scale on which
#' their variance no longer depends on the [population index][vismeteor::vmgeom]
#' `r`, provided the magnitudes follow the geometric model. The limiting
#' magnitude `lm` enters the mapping, so it needs no further attention in the
#' subsequent analysis. The variance of the transformed magnitudes is close to
#' `1.0`, and their mean estimates `r` through `vmgeom_vst_to_r`.
#'
#' The transformation is a monotone rescaling of the underlying quantity used
#' by the rate-based estimator, which is formed from the perception
#' probabilities of two magnitudes one class apart. `vmgeom_vst_to_r` maps the
#' mean of the transformed magnitudes back onto `r`. See `vignette("vmgeom")`
#' for the derivation, the properties of the resulting estimator and a
#' comparison with the other methods the package offers.
#'
#' Two consequences are relevant when calling these functions. Since that
#' quantity lies between `0` and `1`, the transformed magnitudes are bounded as
#' well, and the upper bound corresponds to \eqn{r = 1}. And once both
#' perception probabilities are close to `1.0`, it no longer distinguishes how
#' bright a meteor actually was, so bright meteors are pooled into a single
#' class.
#'
#' The variance is close to `1.0` over \eqn{1.4 \le r \le 4.0}, which covers the
#' \eqn{1.7 \le r \le 3.3} met in practice, and degrades outside that window.
#' This matters when the transformed magnitudes are used as a response in a
#' linear model; it does not affect reading a single `r` off their mean.
#'
#' The back-transformation carries a small systematic deviation that does not
#' shrink as the sample grows. In the practical range it stays below `1.5%` of
#' `r`, and is negligible against the random error, but it grows for larger `r`
#' and comes to dominate in very large samples. Where this matters, the
#' rate-based estimator shown in `vignette("vmgeom")` is unbiased for `1/r`,
#' and `vmgeom_glm` fits the exact likelihood.
#'
#' The numerical form of the transformation is version-specific and may change
#' substantially in future releases. Do not rely on equality of transformed
#' values across package versions.
#'
#' @return
#' * `vmgeom_vst_from_magn`: numeric value, the transformed meteor magnitude.
#' * `vmgeom_vst_to_r`: numeric value of the population index `r`, derived from
#'   the mean of `tm`.
#'
#' The argument `deriv_degree` can be used to apply the delta method.
#' If `log = TRUE`, the logarithm of `r` is returned.
#'
#' `vmgeom_vst_to_r` is defined for every `tm` the transformation can produce
#' and returns `NA` only for values it cannot, that is negative ones and those
#' above the upper bound; `tm = 0` yields `Inf`. The recovered `r` is strictly
#' monotone in `tm`, so an implausible estimate — as sparse data may produce in
#' a predictive model — comes back too large, but never as `NA`.
#'
#' @note
#' The transformation is based on the perception probabilities produced by
#' [vismeteor::vmperception].
#'
#' The transformation uses the perception probability one magnitude class below
#' `m`, which vanishes for `lm - m <= 0.5`. Such magnitudes are transformed to
#' `0`, including the faintest class still inside the support of the model.
#' This is the correct value — the ratio is genuinely zero there — and those
#' meteors belong in the mean; do not exclude them. Magnitudes outside the
#' support, that is `lm - m <= -0.5`, likewise yield `0` and should as usual be
#' removed beforehand.
#'
#' @seealso [vismeteor::vmgeom] [vismeteor::vmperception] [vismeteor::vmgeom_glm]
#'   `vignette("vmgeom")`
#' @examples
#' N <- 100
#' r <- 2.0
#' limmag <- 6.3
#'
#' # Simulate magnitudes
#' m <- rvmgeom(N, limmag, r)
#'
#' # Variance-stabilizing transformation
#' tm <- vmgeom_vst_from_magn(m, limmag)
#' tm_mean <- mean(tm)
#' tm_var <- var(tm)
#'
#' # Estimator for r from the transformed mean
#' r_hat <- vmgeom_vst_to_r(tm_mean)
#'
#' # Derivative dr/d(tm) at tm_mean (needed for the delta method)
#' dr_dtm <- vmgeom_vst_to_r(tm_mean, deriv_degree = 1L)
#'
#' # Variance of the sample mean of tm
#' var_tm.mean <- tm_var / N
#'
#' # Delta method: variance and standard error of r_hat
#' var_r.hat <- (dr_dtm^2) * var_tm.mean
#' se_r.hat <- sqrt(var_r.hat)
#'
#' # Results
#' print(r_hat)
#' print(se_r.hat)
#'
#' # The transformation depends on the magnitude only through the difference
#' # to the limiting magnitude. Bright meteors run into the upper bound,
#' # where the transformed value no longer resolves their brightness.
#' old_par <- par(mfrow = c(1, 1))
#' plot(
#'     function(dm) vmgeom_vst_from_magn(0.0, dm),
#'     -0.5, 10.0,
#'     main = paste(
#'         "variance-stabilizing transformation of",
#'         "visual meteor magnitudes"
#'     ),
#'     col = "blue",
#'     xlab = "dm",
#'     ylab = "tm"
#' )
#' abline(h = vmgeom_vst_from_magn(0.0, 100.0), lty = "dashed")
#'
#' par(old_par)

#' @rdname vmgeom_vst
#' @export
vmgeom_vst_from_magn <- function(m, lm) {
    a <- .vmgeom_vst_params$a
    b <- .vmgeom_vst_params$b

    # ratio of the perception probabilities one magnitude apart;
    # it is zero where the numerator or the denominator vanishes
    dm <- lm - m
    g <- rep(NA_real_, length(dm))
    known <- !is.na(dm)
    if (any(known)) {
        den <- vmperception(dm[known])
        g[known] <- ifelse(den > 0.0, vmperception(dm[known] - 1.0) / den, 0.0)
    }

    a * g^b
}

#' @rdname vmgeom_vst
#' @export
vmgeom_vst_to_r <- function(tm, log = FALSE, deriv_degree = 0L) {
    c <- .vmgeom_vst_params$c
    d <- .vmgeom_vst_params$d

    # `g` is a ratio of a monotonically increasing function and therefore
    # bounded by 0 and 1, which bounds `tm` by 0 and `a`.
    # Both limits follow from the model, not from a calibration range.
    tm[tm < 0.0 | tm > .vmgeom_vst_params$a] <- NA_real_

    if (deriv_degree > 2L) {
        stop(paste("deriv_degree", deriv_degree, "not implemented!"))
    }

    # `log` masks base::log() here, hence the qualified calls below
    if (log) {
        if (2L == deriv_degree) {
            -d / tm^2
        } else if (1L == deriv_degree) {
            d / tm
        } else {
            c + d * base::log(tm)
        }
    } else {
        if (2L == deriv_degree) {
            d * (d - 1.0) * exp(c) * tm^(d - 2.0)
        } else if (1L == deriv_degree) {
            d * exp(c) * tm^(d - 1.0)
        } else {
            exp(c) * tm^d
        }
    }
}

#' Parameters of the variance-stabilizing transformation.
#'
#' `a` and `b` define the transformation, `c` and `d` the back-transformation.
#' Since `g` cannot exceed 1, the factor `a` is at the same time the upper
#' bound of the transformed magnitudes.
#'
#' Derivation of `c` and `d`
#'
#' The transformation is built on the ratio of the perception probabilities of
#' two magnitude classes one apart,
#'
#'     g = f(lm - m - 1) / f(lm - m),   with mean(g) = 1/r,
#'
#' which is the statistic of the rate-based estimator. It is unbiased for
#' `1/r` and, unlike the magnitudes themselves, no longer depends on the
#' limiting magnitude. The transformation is the power transformation
#'
#'     t = a * g^b,
#'
#' whose sole purpose is to make the variance independent of `r`.
#'
#' Inverting `t = a * g^b` algebraically would only be correct per meteor.
#' Applied to a mean it is not, because `b` is not 1 and therefore
#' `mean(g^b) != mean(g)^b`. Since the estimator averages the transformed
#' magnitudes, the back-transformation has to map `mean(t)` onto `mean(g)`,
#' and that relation is obtained by regression rather than by rearranging:
#'
#'     lm(log(mean(g)) ~ log(mean(t)))  =>  log(g) = k0 + k1 * log(t),
#'
#' evaluated over a grid of `r` and `lm` under the geometric model. With
#' `log(r) = -log(g)` this gives
#'
#'     log(r) = -k0 - k1 * log(t) = c + d * log(t),
#'
#' hence `c = -k0` and `d = -k1`. See `inst/derivation/vmgeom_vst.R`.
#'
#' @noRd
.vmgeom_vst_params <- (\() {
    list(
        a = 4.52,
        b = 0.68,
        c = 2.034928,
        d = -1.339972
    )
})()
