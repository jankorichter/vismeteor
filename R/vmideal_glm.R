#' @name vmideal_glm
#' @aliases predict.vmideal_glm
#' @title Generalized Linear Model for the Ideal Distribution of Visual Meteor Magnitudes
#' @description
#' Fits the [ideal distribution][vismeteor::vmideal] of visual meteor magnitudes
#' as a generalized linear model, so that the location parameter `psi` can be
#' estimated as a function of covariates.
#' @param perception_fun function; perception probability function (optional).
#'     Default is [vismeteor::vmperception].
#' @param psi_range numeric; the range of location parameters considered when
#'     inverting the mean. The default leaves `psi` unrestricted.
#' @param formula an object of class [stats::formula]; the response must be a
#'     two-column matrix `cbind(magn, lim_magn)`.
#' @param data a data.frame containing the variables of the model.
#' @param object an object returned by `vmideal_glm`.
#' @param newdata optional data.frame in which to look for variables with which
#'     to predict. If omitted, the fitted linear predictors are used.
#' @param type character; the scale on which to return predictions. `"psi"`
#'     (the default) returns the location parameter, `"link"` the linear
#'     predictor. The link is the identity, so both agree, except where the
#'     data no longer determine `psi`: `"psi"` is `Inf` there, while `"link"`
#'     remains the finite value the fit stopped at.
#' @param se.fit logical; if `TRUE`, standard errors and confidence limits are
#'     returned in addition to the predictions.
#' @param level numeric; the confidence level of the returned interval.
#' @param ... further arguments passed to [stats::glm] respectively
#'     [stats::predict.glm].
#' @details
#' In the [ideal distribution][vismeteor::vmideal] of visual meteor magnitudes,
#' the location parameter `psi` is usually estimated as a single global value.
#' `vmideal_glm` allows `psi` to depend on covariates instead, such as the solar
#' longitude or the limiting magnitude.
#'
#' The response must be given as a two-column matrix, where the first column
#' contains the integer meteor magnitudes and the second column the
#' corresponding limiting magnitudes:
#'
#' ```
#' vmideal_glm(cbind(magn, lim_magn) ~ x, data = magnitudes)
#' ```
#'
#' The limiting magnitude is a fixed property of each observation and is never
#' estimated. Passing it as part of the response ensures that it is subsetted
#' together with the remaining data whenever [stats::glm] drops rows, for
#' example through `subset` or `na.action`.
#'
#' The link is the identity, since `psi` is unrestricted on the real line. The
#' fit is an ordinary [stats::glm] object, so `summary`, `anova`, `AIC`/`BIC`
#' and [vismeteor::select_knots] apply unchanged.
#'
#' Visual magnitude observations are usually reported as counts per magnitude
#' class rather than as individual meteors. Such data can be passed in
#' aggregated form by supplying the counts as `weights`:
#'
#' ```
#' vmideal_glm(cbind(magn, lim_magn) ~ x, data = magnitudes, weights = Freq)
#' ```
#'
#' This is equivalent to repeating every row `Freq` times, and yields identical
#' coefficients, standard errors, deviance and `AIC`. Fractional counts should
#' be rounded with [vismeteor::vmtable] beforehand.
#'
#' Unlike [vismeteor::vmgeom_glm], the fit is not an exact maximum likelihood
#' estimate but a quasi-likelihood one: the mean of the ideal distribution is
#' not a sufficient statistic for `psi`. It is consistent, and for a few
#' thousand meteors with `psi` at or below the limiting magnitude it differs
#' from the maximum likelihood estimate by a few hundredths. For small samples,
#' or a `psi` well above the limiting magnitude, the difference can reach a
#' whole magnitude, and the standard errors are then no longer conservative.
#' Use [stats::optim] on [dvmideal][vismeteor::vmideal] directly when an exact
#' maximum likelihood estimate of a single global `psi` is required.
#' `vignette("vmideal")` sets out the reason and the size of the effect.
#'
#' @section Range of psi:
#' `psi` is unrestricted, and the default `psi_range` leaves it that way. What
#' the data can resolve is governed by the distance between `psi` and the
#' limiting magnitude: `psi` is estimated reliably as long as it does not
#' exceed the limiting magnitude by more than about two to three magnitudes.
#' For visual meteor observations, where `psi` is typically around `4` to `6`
#' and limiting magnitudes lie near `6`, this leaves ample room.
#'
#' Beyond that the mean of the distribution flattens, and once it has reached
#' the limit it converges to, no data can distinguish a larger `psi` from a
#' smaller one. `predict(type = "psi")` then returns `Inf`. That is a result
#' rather than a failure: the magnitudes are geometric with
#' \eqn{r = 10^{0.4}}, and `Inf` says that `psi` is bounded below by the fit
#' rather than determined by it. It also arises when the mean magnitude of the
#' data lies above anything the model can produce, whether by chance in a small
#' sample or because the ideal distribution does not describe them. The linear
#' predictor stays finite throughout, so `type = "link"`, `summary` and `anova`
#' remain usable. Just short of that point the iteration may run out of steps;
#' such a fit is reported with a warning rather than returned as if it had
#' converged.
#'
#' An infinite `psi_range` bound is permitted here, unlike in
#' [dvmideal][vismeteor::vmideal], where an infinite `psi` is rejected.
#'
#' Predictions and their uncertainties are obtained on the scale of the linear
#' predictor. Since the link is the identity, the confidence limits are the
#' ordinary symmetric limits `fit +/- z * se.fit`.
#'
#' @return
#' * `vmideal_glm`: an object of class `vmideal_glm` inheriting from
#'   [stats::glm].
#' * `predict`: a numeric vector of predictions, or a list with the components
#'   `fit`, `se.fit`, `lwr` and `upr` if `se.fit = TRUE`.
#' @seealso [vismeteor::vmideal] [vismeteor::vmgeom_glm] [vismeteor::vmperception]
#'   [stats::glm] [stats::predict.glm]
#' @examples
#' # Simulate magnitudes whose location parameter depends on a covariate
#' set.seed(1)
#' lim_magn <- rep(c(5.8, 6.1, 6.4), each = 120)
#' x <- rep(seq(-1.5, 1.5, length.out = 12), each = 30)
#' psi <- 4.0 + 0.8 * x
#' magn <- mapply(\(lim_magn, psi) rvmideal(1L, lim_magn, psi), lim_magn, psi)
#' obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)
#'
#' # Fit the model
#' fit <- vmideal_glm(cbind(magn, lim_magn) ~ x, data = obs)
#' summary(fit)
#'
#' # Location parameter at selected covariate values
#' newdata <- data.frame(x = c(-1, 0, 1))
#' predict(fit, newdata) # psi
#'
#' # ... including confidence limits
#' predict(fit, newdata, se.fit = TRUE)
NULL

#' family object of the ideal distribution
#'
#' Builds the [stats::family] object that `vmideal_glm()` hands to
#' [stats::glm]. It is internal because the object carries the limiting
#' magnitudes of the fit it is used for and must therefore not be reused, and
#' because `vmideal_glm()` already exposes every argument it accepts.
#'
#' @noRd
.vmideal_family <- function(
    perception_fun = vmperception,
    psi_range = c(-Inf, Inf)
) {
    perception_fun <- match.fun(perception_fun)

    # Unlike the population index of the geometric model, `psi` is unbounded,
    # so an infinite bound is meaningful here. An infinite *lower* bound is
    # still only a limit, never a value the fit may take.
    if (2L != length(psi_range) || anyNA(psi_range) || Inf == psi_range[1L] || -Inf == psi_range[2L]) {
        stop("psi_range must be two values!")
    }

    if (psi_range[1L] >= psi_range[2L]) {
        stop("psi_range must be increasing!")
    }

    # Holds the limiting magnitudes of the current fit. It is filled by
    # `initialize`, which is evaluated inside the frame of stats::glm.fit.
    env <- new.env(parent = emptyenv())
    env$lim_magn <- NULL
    env$mu <- NULL
    env$psi <- NULL

    moments <- \(lim_magn, psi) .vmideal_moments(lim_magn, psi, perception_fun)

    # Inverting the mean numerically is by far the most expensive operation of
    # the family. During the iteration `mu` is always the result of the
    # preceding `linkinv()` call, whose location parameter is therefore known
    # exactly. Only that mapping is remembered, so a genuinely new `mu` still
    # falls back to the numerical inversion.
    psi_from_mu <- \(mu, lim_magn) {
        if (!is.null(env$mu) && length(env$mu) == length(mu) && identical(env$mu, mu)) {
            return(env$psi)
        }

        .vmideal_psi_from_mu(mu, lim_magn, perception_fun, psi_range)
    }

    structure(
        list(
            family = "vmideal",
            link = "identity",
            # The variance is fully determined by `psi` and the limiting
            # magnitude, so there is no free dispersion parameter to estimate.
            # Fixing it also makes the standard errors of a fit weighted by
            # magnitude counts agree with those of the replicated data.
            dispersion = 1.0,
            linkfun = \(mu) psi_from_mu(mu, env$lim_magn),
            linkinv = \(eta) {
                # During the fit `eta` always covers the observations of the
                # model frame. A different length means the family is being
                # used outside the fit, typically by
                # `stats::predict.glm(type = "response")`, which does not know
                # the limiting magnitudes of the new data and would silently
                # recycle them.
                if (!is.null(env$lim_magn) && length(eta) != length(env$lim_magn)) {
                    stop(paste(
                        "the limiting magnitudes of the new data are unknown;",
                        "fit the model with vmideal_glm() and use predict() instead"
                    ), call. = FALSE)
                }

                mu <- moments(env$lim_magn, eta)$mu
                env$mu <- mu
                env$psi <- eta

                mu
            },
            variance = \(mu) moments(env$lim_magn, psi_from_mu(mu, env$lim_magn))$var,
            mu.eta = \(eta) {
                # The mean is obtained from `pmideal()`, a monotone spline whose
                # relative accuracy in `psi` is about 1e-9 -- some seven orders
                # of magnitude above machine precision. A step of the size used
                # for the geometric model (1e-6) drowns in that noise and even
                # returns the wrong sign. This step keeps the truncation error
                # below one percent while staying clear of the noise floor.
                h <- 0.05
                derivative <- (moments(env$lim_magn, eta + h)$mu -
                    moments(env$lim_magn, eta - h)$mu) / (2.0 * h)

                # Once the mean has saturated the difference is pure rounding
                # noise, and its sign is arbitrary. The iteration divides by
                # this value, so leaving it there makes the step explode in a
                # random direction. The mean is strictly increasing in `psi`,
                # so the derivative is positive by construction, and flooring
                # it keeps the step bounded: the iteration then walks into the
                # saturated region instead of jumping across it, and settles
                # against the bound where `predict()` reports `psi` as
                # infinite.
                pmax(derivative, 0.05)
            },
            dev.resids = \(y, mu, wt) {
                psi <- psi_from_mu(mu, env$lim_magn)
                -2.0 * wt * .vmideal_log_dens(y, env$lim_magn, psi, perception_fun)
            },
            aic = \(y, n, mu, wt, dev) sum(dev),
            initialize = substitute(
                expression({
                    if (is.null(dim(y)) || 2L != NCOL(y)) {
                        stop(paste(
                            "for the 'vmideal' family, y must be a two column matrix",
                            "where col 1 is the meteor magnitude and col 2 the limiting magnitude"
                        ))
                    }
                    assign("lim_magn", as.numeric(y[, 2L]), envir = ENV)
                    y <- as.numeric(y[, 1L])
                    n <- rep.int(1, nobs)
                    mustart <- MUSTART(y, get("lim_magn", envir = ENV))
                }),
                list(ENV = env, MUSTART = \(y, lim_magn) {
                    range <- .vmideal_mu_range(lim_magn, perception_fun, psi_range)
                    interior <- .vmideal_moments(
                        lim_magn,
                        .vmideal_eta_range(lim_magn, psi_range)[2L] - 2.0,
                        perception_fun
                    )$mu

                    # Magnitudes fainter on average than the model can produce
                    # have their likelihood maximised where the mean has
                    # saturated. The iteration cannot travel there from the
                    # response, because a step towards an unattainable mean
                    # overshoots every bound, so it starts at that end instead.
                    if (mean(y) >= mean(range$upper)) {
                        return(interior)
                    }

                    pmin(pmax(y, range$lower), interior)
                })
            )[[2L]],
            validmu = \(mu) {
                if (!all(is.finite(mu))) {
                    return(FALSE)
                }

                # Only means that a location parameter within `psi_range` can
                # produce are attainable. The upper bound is the mean of the
                # geometric model, which the ideal distribution reaches as
                # `psi` grows without bound.
                range <- .vmideal_mu_range(env$lim_magn, perception_fun, psi_range)
                tol <- .Machine$double.eps^0.5 * (1.0 + abs(mu))

                all(mu >= range$lower - tol & mu <= range$upper + tol)
            },
            valideta = \(eta) {
                if (!all(is.finite(eta))) {
                    return(FALSE)
                }

                if (is.null(env$lim_magn)) {
                    return(TRUE)
                }

                # With an identity link the linear predictor *is* `psi`, so
                # this is the only place where `psi_range` can bind on the
                # coefficients. Without it the linear predictor runs off to
                # where the mean no longer identifies `psi` -- the derivative
                # is zero there, so the step of the iteration is unbounded --
                # and the fit reports an arbitrary value as converged.
                #
                # Magnitudes fainter than the model can produce push the
                # iteration against the upper bound rather than through it, so
                # a little room above it keeps the step halving of `glm.fit()`
                # working; `predict()` reports anything up there as infinite.
                range <- .vmideal_eta_range(env$lim_magn, psi_range)

                all(eta >= range[1L] & eta <= range[2L] + 2.0)
            }
        ),
        class = "family"
    )
}

#' @rdname vmideal_glm
#' @export
vmideal_glm <- function(formula, data, ...,
    perception_fun = vmperception,
    psi_range = c(-Inf, Inf)
) {
    cl <- match.call()

    # Build the stats::glm() call rather than forwarding `...`, so that
    # non-standard arguments such as `subset` and `weights` are evaluated in
    # the caller's frame.
    fit_call <- cl
    fit_call[[1L]] <- quote(stats::glm)
    fit_call$perception_fun <- NULL
    fit_call$psi_range <- NULL
    fit_call$family <- quote(family)

    env <- new.env(parent = parent.frame())
    assign(
        "family",
        .vmideal_family(perception_fun = perception_fun, psi_range = psi_range),
        envir = env
    )

    # Reaching the bound is a result, not a failure, and `predict()` reports it
    # as an infinite `psi`. Only the warnings describing the step halving that
    # gets there are mechanism rather than result and may be dropped. A failure
    # to converge is neither, and has to stay visible.
    quietly <- \(expr) {
        withCallingHandlers(
            expr,
            warning = \(w) {
                expected <- c("boundary value", "step size truncated")
                if (any(vapply(expected, grepl, logical(1L), conditionMessage(w), fixed = TRUE))) {
                    invokeRestart("muffleWarning")
                }

                # The iteration runs out of steps where the mean has flattened
                # so far that it barely moves with `psi`. Saying that is more
                # use than reporting the mechanism.
                if (grepl("did not converge", conditionMessage(w), fixed = TRUE)) {
                    warning(paste(
                        "the fit did not converge; psi is close to the point at which",
                        "the mean of the distribution no longer identifies it, so the",
                        "estimate is unreliable and its standard error understates that"
                    ), call. = FALSE)

                    invokeRestart("muffleWarning")
                }
            }
        )
    }

    fit <- tryCatch(
        quietly(eval(fit_call, env)),
        error = \(e) {
            # Magnitudes fainter than the model can produce have their
            # likelihood maximised where the mean has saturated. The iteration
            # cannot travel there from the default starting value, because a
            # step towards an unattainable mean overshoots every bound, so it
            # is restarted from that end instead.
            start <- .vmideal_bound_start(fit_call, env, perception_fun, psi_range)
            if (is.null(start)) {
                stop(e)
            }

            fit_call$start <- start
            quietly(eval(fit_call, env))
        }
    )
    fit$call <- cl
    class(fit) <- c("vmideal_glm", class(fit))

    fit
}

#' @rdname vmideal_glm
#' @export
predict.vmideal_glm <- function(
    object,
    newdata,
    type = c("psi", "link"),
    se.fit = FALSE,
    level = 0.95,
    ...
) {
    type <- match.arg(type)

    if (level <= 0.0 || level >= 1.0) {
        stop("level must be between 0.0 and 1.0!")
    }

    # Always predict on the link scale. `linkinv` would read the limiting
    # magnitudes of the fit instead of those of `newdata`. As the link is the
    # identity, that scale already is `psi`.
    pred <- if (missing(newdata) || is.null(newdata)) {
        stats::predict.glm(object, type = "link", se.fit = TRUE)
    } else {
        stats::predict.glm(object, newdata = newdata, type = "link", se.fit = TRUE)
    }

    to_psi <- if ("psi" == type) {
        lim_magn <- .vmideal_predict_lim_magn(
            object, if (missing(newdata)) NULL else newdata
        )
        \(eta) .vmideal_unbound(eta, lim_magn, object)
    } else {
        identity
    }

    if (!se.fit) {
        return(to_psi(pred$fit))
    }

    z <- stats::qnorm(1.0 - (1.0 - level) / 2.0)

    list(
        fit = to_psi(pred$fit),
        se.fit = pred$se.fit,
        lwr = to_psi(pred$fit - z * pred$se.fit),
        upr = to_psi(pred$fit + z * pred$se.fit)
    )
}

#' starting coefficients at the upper bound
#'
#' Returns `NULL` unless the mean magnitude of the data exceeds what the model
#' can produce, in which case the maximum of the likelihood lies at the bound
#' and the iteration has to be started from there. Only an intercept is placed
#' at the bound; any further coefficient starts at zero, which is what
#' `glm.fit()` would have used anyway.
#'
#' @noRd
.vmideal_bound_start <- function(fit_call, env, perception_fun, psi_range) {
    frame_call <- fit_call
    frame_call[[1L]] <- quote(stats::model.frame)
    frame_call$family <- NULL
    frame_call$start <- NULL
    frame_call$control <- NULL

    frame <- tryCatch(eval(frame_call, env), error = \(e) NULL)
    if (is.null(frame)) {
        return(NULL)
    }

    response <- stats::model.response(frame)
    if (is.null(dim(response)) || 2L != NCOL(response)) {
        return(NULL)
    }

    magn <- as.numeric(response[, 1L])
    lim_magn <- as.numeric(response[, 2L])
    upper <- .vmgeom_moments(lim_magn, .vmideal_r_limit, perception_fun)$mu

    if (mean(magn) < mean(upper)) {
        return(NULL)
    }

    p <- NCOL(stats::model.matrix(attr(frame, "terms"), frame))
    start <- rep(0.0, p)
    start[1L] <- .vmideal_eta_range(lim_magn, psi_range)[2L]

    start
}

#' limiting magnitudes a prediction refers to
#'
#' The limiting magnitude is part of the response, so `newdata` carries it only
#' if the caller supplied that column. Where it is missing, the widest limiting
#' magnitude of the fit is the conservative stand-in: it places the threshold
#' as high as any observation could justify.
#'
#' @noRd
.vmideal_predict_lim_magn <- function(object, newdata) {
    lim_magn <- object$model[[1L]][, 2L]

    if (is.null(newdata)) {
        return(lim_magn)
    }

    name <- all.vars(stats::formula(object)[[2L]])[2L]
    if (!is.null(newdata[[name]])) {
        return(as.numeric(newdata[[name]]))
    }

    max(lim_magn)
}

#' report an unidentified location parameter as infinite
#'
#' The mean of the ideal distribution rises towards the geometric limit as
#' `psi` grows, and reaches it to within rounding at
#' the upper end of `.vmideal_eta_range()`. A fit that ends there has not failed:
#' the magnitudes are geometric with the limiting population index, and every
#' larger `psi` describes them equally well. The estimate is unbounded, which
#' `Inf` states plainly, rather than a bound that would read like a number the
#' data supported.
#'
#' @noRd
.vmideal_unbound <- function(eta, lim_magn, object) {
    psi_range <- environment(object$family$linkinv)$psi_range

    if (is.null(lim_magn) || 0L == length(lim_magn)) {
        return(eta)
    }

    # Identifiability is governed by the distance between `psi` and the
    # limiting magnitude of the observation being predicted for, so the
    # threshold follows that row rather than the widest limiting magnitude in
    # the fit. Otherwise unrelated rows could decide whether `psi` counts as
    # identified here.
    lim_magn <- rep(lim_magn, length.out = length(eta))
    bound <- vapply(
        lim_magn, \(lim_magn) .vmideal_eta_range(lim_magn, psi_range)[2L],
        numeric(1L)
    )

    # The mean is already flat well before the bound itself, so the comparison
    # has to allow for the last stretch in which the iteration merely creeps
    # towards it. One magnitude of slack covers that without reaching into the
    # region where `psi` is still resolved.
    eta[eta >= bound - 1.0] <- Inf

    eta
}

#' population index of the geometric limit
#'
#' As `psi` grows without bound, the ideal distribution converges to the
#' geometric model with this population index.
#'
#' @noRd
.vmideal_r_limit <- 10^0.4

#' how far from the limiting magnitude the iteration may wander
#'
#' Both entries are distances from the limiting magnitude, but they bound the
#' iteration for different reasons.
#'
#' `below` is not about identifiability: a small `psi` is estimated perfectly
#' well, the mean simply follows it one to one. What grows is the width of the
#' table, which tracks the distance between `psi` and the limiting magnitude.
#' The bound therefore sits far outside anything observable -- a meteor stream
#' 40 magnitudes below the limiting magnitude would be invisible.
#'
#' `above` is about identifiability. Once `psi` exceeds the limiting magnitude
#' by this much, the mean equals its geometric limit to within rounding, its
#' derivative has dropped into the noise of the spline underlying it, and no
#' data can tell a larger `psi` from a smaller one. Carrying estimates past
#' this point produced arbitrary values instead of a signal that `psi` is
#' unidentified.
#'
#' @noRd
.vmideal_reach <- c(below = 40.0, above = 10.0)

#' bounds of the linear predictor implied by the limiting magnitudes
#'
#' Confining the iteration to that window is what keeps a diverging linear
#' predictor from being reported as a converged estimate.
#'
#' @noRd
.vmideal_eta_range <- function(lim_magn, psi_range) {
    c(
        max(psi_range[1L], min(lim_magn) - .vmideal_reach[["below"]]),
        min(psi_range[2L], max(lim_magn) + .vmideal_reach[["above"]])
    )
}

#' observations that have reached the geometric limit
#'
#' The ideal density is evaluated relative to `psi`, so a large `psi` shifts the
#' whole tabulated window into the far tail, where every term underflows to
#' zero. The normalization then vanishes and the moments would be `0/0`. Since
#' the distribution has long since become geometric at that point -- it is
#' indistinguishable from its limit in double precision from `psi` of about 17
#' onwards -- those observations are handed to the geometric model instead.
#'
#' @noRd
.vmideal_is_geometric <- function(lim_magn, psi, perception_fun) {
    # A missing value is not a limit and must not be silently answered with
    # one. Within the fit `na.action` has removed such rows long before.
    if (anyNA(psi) || anyNA(lim_magn)) {
        stop("NA's are not allowed!")
    }

    degenerate <- !is.finite(psi)

    # The normalization depends on `psi` only through its distance from the
    # limiting magnitude, and stays normal until that distance reaches about
    # 766. Screening arithmetically keeps the ordinary path free of the extra
    # tabulation; only the few rows near the threshold are examined in full.
    suspect <- !degenerate & psi - lim_magn > 700.0

    if (any(suspect)) {
        # A finite `psi` is only replaced once the tabulation has actually
        # collapsed, so that the ordinary path stays bit-for-bit unchanged. The
        # threshold is the smallest normal double: below it the terms of the
        # table lose their trailing digits before the normalization reaches
        # zero, which would otherwise leave a narrow band of degraded moments.
        tab <- .vmideal_table(lim_magn[suspect], psi[suspect], perception_fun)
        norm <- tab$norm[tab$row_of]

        degenerate[suspect] <- !is.finite(norm) | norm < .Machine$double.xmin
    }

    degenerate
}

#' mean and variance of the ideal distribution
#'
#' Vectorized over `lim_magn` and `psi`. The moments are evaluated once per
#' unique combination of both, which is essential for performance since limiting
#' magnitudes are discrete in practice.
#'
#' @noRd
.vmideal_moments <- function(lim_magn, psi, perception_fun) {
    n <- max(length(lim_magn), length(psi))
    lim_magn <- rep(lim_magn, length.out = n)
    psi <- rep(psi, length.out = n)

    # Beyond the reach of double precision the distribution is geometric, and
    # the tabulation below would divide by a vanished normalization.
    geometric <- .vmideal_is_geometric(lim_magn, psi, perception_fun)
    if (any(geometric)) {
        moments <- list(mu = rep(NA_real_, n), var = rep(NA_real_, n))

        limit <- .vmgeom_moments(lim_magn[geometric], .vmideal_r_limit, perception_fun)
        moments$mu[geometric] <- limit$mu
        moments$var[geometric] <- limit$var

        if (any(!geometric)) {
            ideal <- .vmideal_moments(lim_magn[!geometric], psi[!geometric], perception_fun)
            moments$mu[!geometric] <- ideal$mu
            moments$var[!geometric] <- ideal$var
        }

        return(moments)
    }

    tab <- .vmideal_table(lim_magn, psi, perception_fun)

    magn_of <- tab$magn_of
    visible <- tab$visible

    # Contribution of the magnitudes below the tabulated ones. There the
    # distribution is purely geometric in steps of one magnitude and the
    # perception probability is 1, so its moments are known in closed form.
    # Unlike the geometric model, this remainder sits at the bright end and can
    # carry most of the mass, hence it is never negligible.
    ratio <- tab$ratio
    step_mean <- ratio / (1.0 - ratio)
    step_square <- ratio * (1.0 + ratio) / (1.0 - ratio)^2
    beyond_first <- tab$magn_min - 1.0
    beyond_mean <- tab$beyond * (beyond_first - step_mean)
    beyond_square <- tab$beyond * (
        beyond_first^2 - 2.0 * beyond_first * step_mean + step_square
    )

    mean_magn <- (rowSums(visible * magn_of) + beyond_mean) / tab$norm
    square_magn <- (rowSums(visible * magn_of^2) + beyond_square) / tab$norm

    row_of <- tab$row_of

    list(
        mu = mean_magn[row_of],
        var = (square_magn - mean_magn^2)[row_of]
    )
}

#' log density of the ideal distribution
#'
#' Equivalent to `dvmideal(m, lim_magn, psi, log = TRUE)`, but built on the same
#' table as `.vmideal_moments()`. `dvmideal()` groups by the pair of limiting
#' magnitude and location parameter, which degenerates to one group per
#' observation as soon as `psi` varies by row.
#'
#' @noRd
.vmideal_log_dens <- function(m, lim_magn, psi, perception_fun) {
    n <- max(length(m), length(lim_magn), length(psi))
    m <- rep(m, length.out = n)
    lim_magn <- rep(lim_magn, length.out = n)
    psi <- rep(psi, length.out = n)

    # Beyond the reach of double precision the distribution is geometric, and
    # the tabulation below would normalize by a vanished total.
    geometric <- .vmideal_is_geometric(lim_magn, psi, perception_fun)
    if (any(geometric)) {
        logdens <- rep(NA_real_, n)

        logdens[geometric] <- .vmgeom_log_dens(
            m[geometric], lim_magn[geometric], .vmideal_r_limit, perception_fun
        )

        if (any(!geometric)) {
            logdens[!geometric] <- .vmideal_log_dens(
                m[!geometric], lim_magn[!geometric], psi[!geometric], perception_fun
            )
        }

        return(logdens)
    }

    tab <- .vmideal_table(lim_magn, psi, perception_fun)
    row_of <- tab$row_of

    logdens <- rep(-Inf, n)

    # Anything fainter than the rounded limiting magnitude cannot be observed.
    observable <- m <= tab$magn_max[row_of]

    tabulated <- observable & m >= tab$magn_min[row_of]
    if (any(tabulated)) {
        # `magn_of` counts upwards from the brightest tabulated magnitude, so
        # the column of a magnitude is its distance from that lower end.
        cell <- cbind(row_of[tabulated], m[tabulated] - tab$magn_min[row_of[tabulated]] + 1L)
        logdens[tabulated] <- log(tab$visible[cell])
    }

    # Below the tabulated magnitudes the perception probability is 1, leaving
    # the ideal density itself.
    beyond <- observable & m < tab$magn_min[row_of]
    if (any(beyond)) {
        logdens[beyond] <- .dmideal_int(m[beyond], psi[beyond], log = TRUE)
    }

    logdens[observable] <- logdens[observable] -
        log(tab$norm[row_of[observable]])

    logdens
}

#' tabulate the visible part of the ideal distribution
#'
#' The observable magnitudes of one observation form a finite table, running
#' from the limiting magnitude down to the brightness below which the ideal
#' density is purely geometric. This helper builds that table once per distinct
#' combination of limiting magnitude and location parameter. Everything the
#' moments and the density need is derived from it.
#'
#' @noRd
.vmideal_table <- function(lim_magn, psi, perception_fun) {
    n <- max(length(lim_magn), length(psi))
    lim_magn <- rep(lim_magn, length.out = n)
    psi <- rep(psi, length.out = n)

    magn_max <- .vmideal_upper_lm(lim_magn)

    # The table has to reach below `psi - 10`, where `.dmideal_int()` becomes
    # exactly geometric, so that the untabulated remainder has closed-form
    # moments. Fifteen magnitudes below the limiting magnitude are required in
    # any case, because that is where the perception probability reaches 1.
    magn_min <- pmin(magn_max - 15L, floor(psi) - 11L)

    # Observations sharing both the limiting magnitude and the location
    # parameter share their entire table, so each combination is tabulated only
    # once. `row_of` maps every observation back to its row. The key is built
    # numerically from the positions within the distinct values, because
    # pasting one string per observation dominates the runtime of the whole
    # fit. Unlike `paste()`, which rounds to 15 significant digits, this
    # compares bitwise and may therefore split values agreeing to that many
    # digits into separate groups. That is on the safe side: every group is
    # still evaluated exactly, only a redundant tabulation may result.
    lim_magn_key <- match(lim_magn, unique(lim_magn))
    key <- lim_magn_key + (match(psi, unique(psi)) - 1L) * max(lim_magn_key)
    unique_idx <- !duplicated(key)
    row_of <- match(key, key[unique_idx])
    lim_magn <- lim_magn[unique_idx]
    psi <- psi[unique_idx]
    magn_max <- magn_max[unique_idx]
    magn_min <- magn_min[unique_idx]

    # Each row holds the magnitudes of one combination, counted upwards from
    # its brightest tabulated magnitude. The rows are padded to a common width;
    # the padding is switched off through the perception probability below.
    width <- max(magn_max - magn_min) + 1L
    magn_of <- outer(magn_min, seq_len(width) - 1L, "+")

    # Magnitudes above the limiting magnitude are unobservable. They only exist
    # to keep the matrix rectangular and are removed by a perception
    # probability of 0.
    padding <- magn_of > magn_max

    dens <- .dmideal_int(as.numeric(magn_of), rep(psi, times = width))
    perception <- perception_fun(rep(lim_magn, times = width) - as.numeric(magn_of))
    visible <- matrix(dens * perception, nrow = length(psi))
    visible[padding] <- 0.0

    # Below the tabulated magnitudes the ideal density decays geometrically by
    # this factor per magnitude step, which is what gives the remainder its
    # closed form.
    ratio <- 10^-0.4
    beyond <- vismeteor::pmideal(magn_min - 0.5, psi, lower.tail = TRUE)

    list(
        row_of = row_of,
        magn_max = magn_max,
        magn_min = magn_min,
        magn_of = magn_of,
        ratio = ratio,
        visible = visible,
        beyond = beyond,
        norm = rowSums(visible) + beyond
    )
}

#' attainable range of the mean
#'
#' @noRd
.vmideal_mu_range <- function(lim_magn, perception_fun, psi_range) {
    # As `psi` grows without bound the ideal distribution converges to the
    # geometric model, which `.vmideal_moments()` returns for an infinite
    # `psi`, so the upper limit of the mean is exact. Its lower counterpart
    # diverges instead and has to be supplied here.
    upper <- .vmideal_moments(lim_magn, psi_range[2L], perception_fun)$mu

    lower <- if (is.finite(psi_range[1L])) {
        .vmideal_moments(lim_magn, psi_range[1L], perception_fun)$mu
    } else {
        rep(-Inf, length(lim_magn))
    }

    list(lower = lower, upper = upper)
}

#' location parameter from the mean
#'
#' Inverts `.vmideal_moments()` numerically. The mean is strictly increasing in
#' `psi`, so a unique solution exists within the bracketing interval.
#'
#' @noRd
.vmideal_psi_from_mu <- function(mu, lim_magn, perception_fun, psi_range) {
    n <- max(length(mu), length(lim_magn))
    mu <- rep(mu, length.out = n)
    lim_magn <- rep(lim_magn, length.out = n)

    # Each distinct pair of mean and limiting magnitude is inverted only once.
    # The key is numeric for the same reason as in `.vmideal_table()`.
    mu_key <- match(mu, unique(mu))
    key <- mu_key + (match(lim_magn, unique(lim_magn)) - 1L) * max(mu_key)
    idx <- !duplicated(key)
    psi <- mapply(\(mu, lim_magn) {
        f <- \(psi) .vmideal_moments(lim_magn, psi, perception_fun)$mu - mu

        if (!is.finite(mu)) {
            return(psi_range[1L])
        }

        # The search is confined to the window in which `psi` is both
        # identifiable and representable. Outside it the mean carries no
        # information about `psi`, so returning the bound is the honest answer.
        range <- .vmideal_eta_range(lim_magn, psi_range)

        if (f(range[2L]) <= 0.0) {
            return(range[2L])
        }
        if (f(range[1L]) >= 0.0) {
            return(range[1L])
        }

        stats::uniroot(f, range, tol = .Machine$double.eps^0.5)$root
    }, mu[idx], lim_magn[idx])

    as.numeric(psi[match(key, key[idx])])
}
