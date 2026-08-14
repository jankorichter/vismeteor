#' @name vmgeom_vst_lm
#' @aliases predict.vmgeom_vst_lm
#' @title Linear Model of Variance-Stabilized Visual Meteor Magnitudes
#' @description
#' Fits visual meteor magnitudes as an ordinary linear model after applying the
#' [variance-stabilizing transformation][vismeteor::vmgeom_vst] of the
#' [geometric model][vismeteor::vmgeom], and returns the estimates on the scale
#' of the population index `r`.
#' @param formula an object of class [stats::formula]; the response must be a
#'     two-column matrix `cbind(magn, lim_magn)`.
#' @param data a data.frame containing the variables of the model.
#' @param object an object returned by `vmgeom_vst_lm`.
#' @param newdata optional data.frame in which to look for variables with which
#'     to predict. If omitted, the fitted values are used.
#' @param type character; the scale on which to return predictions.
#'     One of `"r"`, `"inv_r"`, `"log_r"` or `"tm"`.
#' @param bias_correction logical; if `TRUE`, the second-order term of the
#'     delta method is applied to correct the bias of the non-linear
#'     back-transformation.
#' @param ... further arguments passed to [stats::lm] respectively
#'     [stats::predict.lm].
#' @details
#' The [variance-stabilizing transformation][vismeteor::vmgeom_vst] maps visual
#' meteor magnitudes onto a scale on which their variance no longer depends on
#' the population index `r`. The transformed magnitudes are therefore
#' homoscedastic by construction, which is precisely what an ordinary linear
#' model assumes, and `vmgeom_vst_lm` fits them as such.
#'
#' The response must be given as a two-column matrix, where the first column
#' contains the integer meteor magnitudes and the second column the
#' corresponding limiting magnitudes:
#'
#' ```
#' vmgeom_vst_lm(cbind(magn, lim_magn) ~ x, data = magnitudes)
#' ```
#'
#' Both are needed because the transformation depends on their difference. The
#' limiting magnitude is a fixed property of each observation and is never
#' estimated. Passing it as part of the response also ensures that it is
#' subsetted together with the remaining data whenever [stats::lm] drops rows.
#'
#' Estimating a single global `r` is the intercept-only case:
#'
#' ```
#' vmgeom_vst_lm(cbind(magn, lim_magn) ~ 1, data = magnitudes, weights = Freq)
#' ```
#'
#' Since the back-transformation is a power relation, effects combine additively
#' on the logarithmic scale, and a coefficient translates into a change of
#' \eqn{\log r} up to a constant factor.
#'
#' Visual magnitude observations are usually reported as counts per magnitude
#' class rather than as individual meteors. Such data are passed in aggregated
#' form by supplying the counts as `weights`. One consequence deserves
#' attention: [stats::lm] treats weights as precision weights and normalizes the
#' residual variance by the number of rows, whereas the counts state how many
#' meteors each row stands for. `predict` corrects for this by deriving the
#' residual scale from `sum(weights)` instead, so that standard errors and
#' confidence limits refer to the number of meteors. Without that correction
#' they come out too wide. Unlike [vismeteor::vmgeom_glm], fractional counts
#' need not be rounded beforehand.
#'
#' Predictions are obtained on the transformed scale and mapped back afterwards
#' with [vismeteor::vmgeom_vst_to_r]. The confidence limits are the transformed
#' limits of that scale. Since every transformation offered here is monotonic,
#' these limits are exact rather than an approximation and are asymmetric about
#' the estimate. They are formed by [stats::predict.lm] and therefore rest on
#' the t-distribution, which is exact for a linear model, rather than on the
#' normal approximation used by [vismeteor::vmgeom_glm].
#'
#' The reported `se.fit` is a first-order delta-method value,
#' \eqn{|g'(t)|\,\mathrm{se}(t)}, matching the behaviour of
#' [stats::predict.glm] and of [vismeteor::vmgeom_glm]. Setting
#' `bias_correction = TRUE` additionally applies the second-order term, which
#' corresponds to `deriv_degree = 2L` of [vismeteor::vmgeom_vst_to_r]. It
#' corrects the curvature of the back-transformation and is small at the sample
#' sizes met in practice. It does *not* address the systematic deviation of the
#' back-transformation itself, which does not shrink as the sample grows; see
#' `vignette("vmgeom")`. Where that deviation matters,
#' [vismeteor::vmgeom_glm] fits the exact likelihood.
#'
#' @return
#' * `vmgeom_vst_lm`: an object of class `vmgeom_vst_lm` inheriting from
#'   [stats::lm].
#' * `predict`: a numeric vector of predictions, or the list returned by
#'   [stats::predict.lm] with its components transformed, if `se.fit = TRUE` or
#'   `interval` is given.
#' @seealso [vismeteor::vmgeom_vst] [vismeteor::vmgeom] [vismeteor::vmgeom_glm]
#'   [stats::lm] [stats::predict.lm] `vignette("vmgeom")`
#' @examples
#' # Simulate magnitudes with a constant population index
#' set.seed(1)
#' lim_magn <- rep(c(5.8, 6.1, 6.4), each = 200)
#' magn <- mapply(\(lim_magn) rvmgeom(1L, lim_magn, r = 2.0), lim_magn)
#' obs <- data.frame(lim_magn = lim_magn, magn = magn)
#'
#' # Fit the model
#' fit <- vmgeom_vst_lm(cbind(magn, lim_magn) ~ 1, data = obs)
#'
#' # The model has no covariates, so any single row predicts the global r-value
#' newdata <- data.frame(row.names = "")
#' predict(fit, newdata) # r
#' predict(fit, newdata, type = "inv_r") # 1/r
#' predict(fit, newdata, type = "log_r") # log(r)
#'
#' # ... including confidence limits
#' predict(fit, newdata, interval = "confidence")
NULL

#' @rdname vmgeom_vst_lm
#' @export
vmgeom_vst_lm <- function(formula, data, ...) {
    cl <- match.call()

    # Build the stats::lm() call rather than forwarding `...`, so that
    # non-standard arguments such as `subset` and `weights` are evaluated in
    # the caller's frame.
    fit_call <- cl
    fit_call[[1L]] <- quote(stats::lm)
    fit_call$formula <- .vmgeom_vst_formula(formula)

    # The transformation is applied to the response before `lm()` sees it. The
    # formula therefore names an ordinary variable, which keeps it free of a
    # reference into this namespace: the formula of the returned fit carries
    # the environment of the caller, where such a reference could not be
    # resolved. The environment below supplies the variable instead.
    env <- new.env(parent = environment(formula))
    assign(
        ".vmgeom_vst_tm",
        .vmgeom_vst_response(eval(formula[[2L]], data, environment(formula))),
        envir = env
    )
    environment(fit_call$formula) <- env

    fit <- eval(fit_call, parent.frame())
    fit$call <- cl
    class(fit) <- c("vmgeom_vst_lm", class(fit))

    fit
}

#' @rdname vmgeom_vst_lm
#' @export
predict.vmgeom_vst_lm <- function(
    object,
    newdata,
    type = c("r", "inv_r", "log_r", "tm"),
    bias_correction = FALSE,
    ...
) {
    type <- match.arg(type)

    # `scale` and `df` carry the correction for aggregated counts and are
    # therefore derived here; passing them again would be silently ignored by
    # the rescaling or abort with R's "matched by multiple actual arguments".
    dots <- list(...)
    if (any(c("scale", "df") %in% names(dots))) {
        stop("scale and df are determined by the model and must not be given!")
    }

    # `lm()` normalizes the residual variance by the number of rows, while the
    # weights count meteors. The residual scale is rescaled accordingly, so
    # that standard errors and limits refer to the number of meteors.
    w <- stats::weights(object)
    if (is.null(w)) {
        w <- rep(1.0, length(stats::residuals(object)))
    }
    n_eff <- sum(w)
    df_residual <- n_eff - length(stats::coef(object))
    scale <- sqrt(sum(w * stats::residuals(object)^2) / df_residual)

    # The bias correction needs the variance of the estimate, which
    # `predict.lm` reports only when asked for it.
    dots$se.fit <- isTRUE(dots$se.fit) || bias_correction
    args <- c(list(object), dots, list(scale = scale, df = df_residual))
    if (!missing(newdata) && !is.null(newdata)) {
        args$newdata <- newdata
    }
    pred <- do.call(stats::predict.lm, args)

    trans <- .vmgeom_vst_transform(type)

    # Without `se.fit` and `interval` predict.lm() returns a bare vector, and
    # with `interval` a matrix. Both are still on the transformed scale.
    if (!is.list(pred)) {
        return(.vmgeom_vst_apply(pred, trans, bias_correction, 0.0))
    }

    # The delta method expands around the estimate, so the derivatives are
    # evaluated on the transformed scale, before the back-transformation.
    tm <- .vmgeom_vst_estimate(pred$fit)
    tm_var <- pred$se.fit^2

    pred$fit <- .vmgeom_vst_apply(pred$fit, trans, bias_correction, tm_var)
    pred$se.fit <- if (bias_correction) {
        # Second-order variance, matching the corrected estimate.
        sqrt(trans$deriv(tm)^2 * tm_var + 0.5 * trans$deriv2(tm)^2 * tm_var^2)
    } else {
        abs(trans$deriv(tm)) * pred$se.fit
    }
    pred$residual.scale <- NULL

    # `se.fit` was requested only to obtain the bias correction.
    if (!isTRUE(list(...)$se.fit)) {
        return(pred$fit)
    }

    pred
}

#' apply the back-transformation to a prediction
#'
#' Handles both the vector returned without `interval` and the matrix returned
#' with it. The delta method expands around the estimate itself, so the
#' bias correction applies to the `fit` column only; the limits are transformed
#' as they are, which keeps them exact. Since the transformations are monotone
#' but not all increasing, the limits are reordered afterwards.
#'
#' @noRd
.vmgeom_vst_apply <- function(x, trans, bias_correction, tm_var) {
    if (is.matrix(x) && all(c("fit", "lwr", "upr") %in% colnames(x))) {
        out <- x
        out[, "fit"] <- trans$fun(x[, "fit"])
        if (bias_correction) {
            out[, "fit"] <- out[, "fit"] -
                0.5 * trans$deriv2(x[, "fit"]) * tm_var
        }
        lwr <- trans$fun(x[, "lwr"])
        upr <- trans$fun(x[, "upr"])
        out[, "lwr"] <- pmin(lwr, upr)
        out[, "upr"] <- pmax(lwr, upr)

        return(out)
    }

    out <- trans$fun(x)
    if (bias_correction) {
        out <- out - 0.5 * trans$deriv2(x) * tm_var
    }

    out
}

#' the estimate a standard error belongs to
#'
#' `predict.lm` returns `fit` either as a vector or, with `interval`, as a
#' matrix. `se.fit` is a vector in both cases and refers to the estimate
#' itself, which is the `fit` column of the matrix.
#'
#' @noRd
.vmgeom_vst_estimate <- function(fit) {
    if (is.matrix(fit)) {
        return(fit[, "fit"])
    }

    fit
}

#' transformations of the variance-stabilized magnitudes
#'
#' Every scale is a power of the transformed magnitude, so the derivatives are
#' available in closed form. `vmgeom_vst_to_r` provides them for `r` and
#' `log_r`; `inv_r` is the reciprocal of `r` and is differentiated here.
#'
#' @noRd
.vmgeom_vst_transform <- function(type) {
    switch(type,
        "tm" = list(
            fun = \(tm) tm,
            deriv = \(tm) rep(1.0, length(tm)),
            deriv2 = \(tm) rep(0.0, length(tm))
        ),
        "r" = list(
            fun = \(tm) vmgeom_vst_to_r(tm),
            deriv = \(tm) vmgeom_vst_to_r(tm, deriv_degree = 1L),
            deriv2 = \(tm) vmgeom_vst_to_r(tm, deriv_degree = 2L)
        ),
        "log_r" = list(
            fun = \(tm) vmgeom_vst_to_r(tm, log = TRUE),
            deriv = \(tm) vmgeom_vst_to_r(tm, log = TRUE, deriv_degree = 1L),
            deriv2 = \(tm) vmgeom_vst_to_r(tm, log = TRUE, deriv_degree = 2L)
        ),
        # 1/r = exp(-c) * tm^(-d), the reciprocal of the power relation
        "inv_r" = list(
            fun = \(tm) 1.0 / vmgeom_vst_to_r(tm),
            deriv = \(tm) {
                d <- .vmgeom_vst_params$d
                -d * exp(-.vmgeom_vst_params$c) * tm^(-d - 1.0)
            },
            deriv2 = \(tm) {
                d <- .vmgeom_vst_params$d
                d * (d + 1.0) * exp(-.vmgeom_vst_params$c) * tm^(-d - 2.0)
            }
        )
    )
}

#' formula with the transformed response
#'
#' Replaces `cbind(magn, lim_magn)` by the transformed magnitudes, so that
#' `stats::lm()` fits an ordinary univariate response while the caller states
#' the model in terms of the observed quantities.
#'
#' The response becomes a plain variable rather than a call, because the
#' formula of the returned fit carries the environment of the caller. A call
#' into this namespace could not be resolved there, and `stats::predict.lm()`
#' does not evaluate the response for `newdata` anyway.
#'
#' @noRd
.vmgeom_vst_formula <- function(formula) {
    if (length(formula) < 3L) {
        stop("formula must have a response of the form cbind(magn, lim_magn)!")
    }

    formula[[2L]] <- quote(.vmgeom_vst_tm)

    formula
}

#' transformed response
#'
#' Applies the transformation to the two-column response. Kept as a function of
#' the response rather than of its columns, so that the limiting magnitude is
#' subsetted together with the magnitudes whenever rows are dropped.
#'
#' @noRd
.vmgeom_vst_response <- function(response) {
    response <- as.matrix(response)
    if (2L != ncol(response)) {
        stop("the response must be a two-column matrix cbind(magn, lim_magn)!")
    }

    vmgeom_vst_from_magn(response[, 1L], response[, 2L])
}
