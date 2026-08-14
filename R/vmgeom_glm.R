#' @name vmgeom_glm
#' @aliases predict.vmgeom_glm
#' @title Generalized Linear Model for the Geometric Model of Visual Meteor Magnitudes
#' @description
#' Fits the [geometric model][vismeteor::vmgeom] of visual meteor magnitudes as
#' a generalized linear model, so that the population index `r` can be
#' estimated as a function of covariates, and returns those estimates on the
#' scale of `r`.
#' @param perception_fun function; perception probability function (optional).
#'     Default is [vismeteor::vmperception].
#' @param r_range numeric; the range of population indices considered when
#'     inverting the mean. Used to bracket the numerical inversion and to keep
#'     the moments of the distribution finite.
#' @param formula an object of class [stats::formula]; the response must be a
#'     two-column matrix `cbind(magn, lim_magn)`.
#' @param data a data.frame containing the variables of the model.
#' @param object an object returned by `vmgeom_glm`.
#' @param newdata optional data.frame in which to look for variables with which
#'     to predict. If omitted, the fitted linear predictors are used.
#' @param type character; the scale on which to return predictions.
#'     One of `"r"`, `"inv_r"`, `"log_r"` or `"link"`.
#' @param se.fit logical; if `TRUE`, standard errors and confidence limits are
#'     returned in addition to the predictions.
#' @param level numeric; the confidence level of the returned interval.
#' @param ... further arguments passed to [stats::glm] respectively
#'     [stats::predict.glm].
#' @details
#' In the [geometric model][vismeteor::vmgeom] of visual meteor magnitudes, the
#' population index `r` is usually estimated as a single global value.
#' `vmgeom_glm` allows `r` to depend on covariates instead, such as the solar
#' longitude or limiting magnitude.
#'
#' The response must be given as a two-column matrix, where the first column
#' contains the integer meteor magnitudes and the second column the
#' corresponding limiting magnitudes:
#'
#' ```
#' vmgeom_glm(cbind(magn, lim_magn) ~ x, data = magnitudes)
#' ```
#'
#' The limiting magnitude is a fixed property of each observation and is never
#' estimated. Passing it as part of the response ensures that it is subsetted
#' together with the remaining data whenever [stats::glm] drops rows, for
#' example through `subset` or `na.action`.
#'
#' The link function is
#' \deqn{
#'     \eta = \mathrm{logit}(1/r) = -\log(r - 1),
#' }
#' with the inverse \eqn{r = 1 + \exp(-\eta)}. The quantity \eqn{1/r} is the
#' factor by which the observable rate drops when the limiting magnitude is
#' lowered by one, and thus lies between `0` and `1`, which is why a logit link
#' is the natural choice. Since
#' \eqn{1 + \exp(-\eta) > 1} holds for every finite \eqn{\eta}, this link also
#' enforces the constraint \eqn{r > 1} required by [vismeteor::dvmgeom].
#'
#' The fit is an ordinary [stats::glm] object, so `summary`, `anova`,
#' `AIC`/`BIC` and [vismeteor::select_knots] apply unchanged. `predict`
#' additionally returns the population index itself rather than the linear
#' predictor.
#'
#' Visual magnitude observations are usually reported as counts per magnitude
#' class rather than as individual meteors. Such data can be passed in
#' aggregated form by supplying the counts as `weights`:
#'
#' ```
#' vmgeom_glm(cbind(magn, lim_magn) ~ x, data = magnitudes, weights = Freq)
#' ```
#'
#' This is equivalent to repeating every row `Freq` times, and yields identical
#' coefficients, standard errors, deviance and `AIC`. Fractional counts should
#' be rounded with [vismeteor::vmtable] beforehand.
#'
#' Internally the distribution is evaluated through \eqn{1/r} rather than `r`
#' itself, which keeps every quantity bounded. The
#' population index is nevertheless bracketed by `r_range`, because the mean of
#' the distribution diverges as \eqn{r \to 1}: a population index of `1` implies
#' arbitrarily many bright meteors and therefore has no finite mean. The default
#' lower bound of `1.01` lies far below any population index encountered in
#' practice, which is typically between `1.5` and `3.0`.
#'
#' Predictions and their uncertainties are obtained on the scale of the linear
#' predictor and transformed afterwards. The confidence limits are
#' the transformed limits of the linear predictor. Since every transformation
#' offered here is monotonic, these limits are exact rather than an
#' approximation, they are asymmetric about the estimate, and they cannot fall
#' below `1`. The reported `se.fit` is a first-order delta-method value,
#' \eqn{|g'(\eta)|\,\mathrm{se}(\eta)}, and is intended as a measure of scale;
#' use the limits themselves for inference. A second-order correction is
#' deliberately omitted, as it stays below one percent for
#' \eqn{\mathrm{se}(\eta) \le 0.1}, which covers the sample sizes occurring in
#' practice.
#'
#' @return
#' * `vmgeom_glm`: an object of class `vmgeom_glm` inheriting from
#'   [stats::glm].
#' * `predict`: a numeric vector of predictions, or a list with the components
#'   `fit`, `se.fit`, `lwr` and `upr` if `se.fit = TRUE`.
#' @seealso [vismeteor::vmgeom] [vismeteor::vmperception] [vismeteor::vmgeom_vst_lm]
#'   [stats::glm] [stats::predict.glm]
#' @examples
#' # Simulate magnitudes whose population index depends on a covariate
#' set.seed(1)
#' lim_magn <- rep(c(5.8, 6.1, 6.4), each = 120)
#' x <- rep(seq(-1.5, 1.5, length.out = 12), each = 30)
#' r <- 1 + exp(-(0.2 + 0.4 * x))
#' magn <- mapply(\(lim_magn, r) rvmgeom(1L, lim_magn, r), lim_magn, r)
#' obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)
#'
#' # Fit the model
#' fit <- vmgeom_glm(cbind(magn, lim_magn) ~ x, data = obs)
#' summary(fit)
#'
#' # Population index at selected covariate values
#' newdata <- data.frame(x = c(-1, 0, 1))
#' predict(fit, newdata) # r
#' predict(fit, newdata, type = "inv_r") # 1/r
#' predict(fit, newdata, type = "log_r") # log(r)
#'
#' # ... including confidence limits
#' predict(fit, newdata, se.fit = TRUE)
NULL

#' family object of the geometric model
#'
#' Builds the [stats::family] object that `vmgeom_glm()` hands to
#' [stats::glm]. It is internal because the object carries the limiting
#' magnitudes of the fit it is used for and must therefore not be reused, and
#' because `vmgeom_glm()` already exposes every argument it accepts.
#'
#' @noRd
.vmgeom_family <- function(
    perception_fun = vmperception,
    r_range = c(1.01, 20.0)
) {
    perception_fun <- match.fun(perception_fun)

    if (2L != length(r_range) || !all(is.finite(r_range)) || r_range[1L] <= 1.0) {
        stop("r_range must be two finite values greater than 1.0!")
    }

    if (r_range[1L] >= r_range[2L]) {
        stop("r_range must be increasing!")
    }

    # Holds the limiting magnitudes of the current fit. It is filled by
    # `initialize`, which is evaluated inside the frame of stats::glm.fit.
    env <- new.env(parent = emptyenv())
    env$lim_magn <- NULL
    env$mu <- NULL
    env$r <- NULL

    moments <- \(lim_magn, r) .vmgeom_moments(lim_magn, r, perception_fun)

    # Inverting the mean numerically is by far the most expensive operation of
    # the family. During the iteration `mu` is always the result of the
    # preceding `linkinv()` call, whose population index is therefore known
    # exactly. Only that mapping is remembered, so a genuinely new `mu` still
    # falls back to the numerical inversion.
    r_from_mu <- \(mu, lim_magn) {
        if (!is.null(env$mu) && length(env$mu) == length(mu) && identical(env$mu, mu)) {
            return(env$r)
        }

        .vmgeom_r_from_mu(mu, lim_magn, perception_fun, r_range)
    }

    structure(
        list(
            family = "vmgeom",
            link = "logit(1/r)",
            # The variance is fully determined by `r` and the limiting
            # magnitude, so there is no free dispersion parameter to estimate.
            # Fixing it also makes the standard errors of a fit weighted by
            # magnitude counts agree with those of the replicated data.
            dispersion = 1.0,
            linkfun = \(mu) -log(r_from_mu(mu, env$lim_magn) - 1.0),
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
                        "fit the model with vmgeom_glm() and use predict() instead"
                    ), call. = FALSE)
                }

                r <- 1.0 + exp(-eta)
                mu <- moments(env$lim_magn, r)$mu
                env$mu <- mu
                env$r <- r

                mu
            },
            variance = \(mu) moments(env$lim_magn, r_from_mu(mu, env$lim_magn))$var,
            mu.eta = \(eta) {
                h <- 1e-6
                (moments(env$lim_magn, 1.0 + exp(-(eta + h)))$mu -
                    moments(env$lim_magn, 1.0 + exp(-(eta - h)))$mu) / (2.0 * h)
            },
            dev.resids = \(y, mu, wt) {
                r <- r_from_mu(mu, env$lim_magn)
                -2.0 * wt * .vmgeom_log_dens(y, env$lim_magn, r, perception_fun)
            },
            aic = \(y, n, mu, wt, dev) sum(dev),
            initialize = substitute(
                expression({
                    if (is.null(dim(y)) || 2L != NCOL(y)) {
                        stop(paste(
                            "for the 'vmgeom' family, y must be a two column matrix",
                            "where col 1 is the meteor magnitude and col 2 the limiting magnitude"
                        ))
                    }
                    assign("lim_magn", as.numeric(y[, 2L]), envir = ENV)
                    y <- as.numeric(y[, 1L])
                    n <- rep.int(1, nobs)
                    mustart <- MUSTART(y, get("lim_magn", envir = ENV))
                }),
                list(ENV = env, MUSTART = \(y, lim_magn) {
                    range <- .vmgeom_mu_range(lim_magn, perception_fun, r_range)
                    pmin(pmax(y, range$lower), range$upper)
                })
            )[[2L]],
            validmu = \(mu) {
                if (!all(is.finite(mu))) {
                    return(FALSE)
                }

                # Only means that a population index within `r_range` can
                # produce are attainable. A mean sitting exactly on a bound is
                # reached whenever the iteration runs into `r_range`, so the
                # comparison has to tolerate rounding.
                range <- .vmgeom_mu_range(env$lim_magn, perception_fun, r_range)
                tol <- .Machine$double.eps^0.5 * (1.0 + abs(mu))

                all(mu >= range$lower - tol & mu <= range$upper + tol)
            },
            valideta = \(eta) all(is.finite(eta))
        ),
        class = "family"
    )
}

#' @rdname vmgeom_glm
#' @export
vmgeom_glm <- function(formula, data, ...,
    perception_fun = vmperception,
    r_range = c(1.01, 20.0)
) {
    cl <- match.call()

    # Build the stats::glm() call rather than forwarding `...`, so that
    # non-standard arguments such as `subset` and `weights` are evaluated in
    # the caller's frame.
    fit_call <- cl
    fit_call[[1L]] <- quote(stats::glm)
    fit_call$perception_fun <- NULL
    fit_call$r_range <- NULL
    fit_call$family <- quote(family)

    env <- new.env(parent = parent.frame())
    assign(
        "family",
        .vmgeom_family(perception_fun = perception_fun, r_range = r_range),
        envir = env
    )

    fit <- eval(fit_call, env)
    fit$call <- cl
    class(fit) <- c("vmgeom_glm", class(fit))

    fit
}

#' @rdname vmgeom_glm
#' @export
predict.vmgeom_glm <- function(
    object,
    newdata,
    type = c("r", "inv_r", "log_r", "link"),
    se.fit = FALSE,
    level = 0.95,
    ...
) {
    type <- match.arg(type)

    if (level <= 0.0 || level >= 1.0) {
        stop("level must be between 0.0 and 1.0!")
    }

    # Always predict on the link scale. `linkinv` would read the limiting
    # magnitudes of the fit instead of those of `newdata`.
    pred <- if (missing(newdata) || is.null(newdata)) {
        stats::predict.glm(object, type = "link", se.fit = TRUE)
    } else {
        stats::predict.glm(object, newdata = newdata, type = "link", se.fit = TRUE)
    }
    eta <- pred$fit

    trans <- .vmgeom_eta_transform(type)
    est <- trans$fun(eta)

    if (!se.fit) {
        return(est)
    }

    z <- stats::qnorm(1.0 - (1.0 - level) / 2.0)
    lower <- trans$fun(eta - z * pred$se.fit)
    upper <- trans$fun(eta + z * pred$se.fit)

    list(
        fit = est,
        se.fit = abs(trans$deriv(eta)) * pred$se.fit,
        lwr = pmin(lower, upper),
        upr = pmax(lower, upper)
    )
}

#' transformations of the linear predictor
#'
#' @noRd
.vmgeom_eta_transform <- function(type) {
    switch(type,
        "link" = list(
            fun = \(eta) eta,
            deriv = \(eta) rep(1.0, length(eta))
        ),
        "r" = list(
            fun = \(eta) 1.0 + exp(-eta),
            deriv = \(eta) -exp(-eta)
        ),
        "inv_r" = list(
            fun = \(eta) stats::plogis(eta),
            deriv = \(eta) exp(-eta) / (1.0 + exp(-eta))^2
        ),
        "log_r" = list(
            fun = \(eta) log1p(exp(-eta)),
            deriv = \(eta) -exp(-eta) / (1.0 + exp(-eta))
        )
    )
}

#' mean and variance of the geometric model
#'
#' Vectorized over `lim_magn` and `r`. The moments are evaluated once per unique
#' combination of both, which is essential for performance since limiting
#' magnitudes are discrete in practice.
#'
#' @noRd
.vmgeom_moments <- function(lim_magn, r, perception_fun) {
    tab <- .vmgeom_table(lim_magn, r, perception_fun)

    r_inv <- tab$r_inv
    rate_norm <- tab$rate_norm
    magn_of <- tab$magn_of

    # Contribution of the magnitudes beyond the tabulated ones, where the
    # distribution is purely geometric and its moments are known in closed form.
    beyond_first <- tab$magn_max + 1.0
    beyond_mean <- tab$beyond * (beyond_first + r_inv / rate_norm)
    beyond_square <- tab$beyond * (
        beyond_first^2 + (2.0 * beyond_first + 1.0) * r_inv / rate_norm +
            2.0 * r_inv^2 / rate_norm^2
    )

    mean_magn <- (rowSums(tab$visible * magn_of) + beyond_mean) / tab$norm
    square_magn <- (rowSums(tab$visible * magn_of^2) + beyond_square) / tab$norm

    # Magnitudes are counted downwards from the rounded limiting magnitude, so
    # the mean has to be mirrored back onto the magnitude scale. The variance is
    # unaffected by that reflection.
    row_of <- tab$row_of

    list(
        mu = tab$lim_magn_round - mean_magn[row_of],
        var = (square_magn - mean_magn^2)[row_of]
    )
}

#' log density of the geometric model
#'
#' Equivalent to `dvmgeom(m, lim_magn, r, log = TRUE)`, but built on the same
#' table as `.vmgeom_moments()`. `dvmgeom()` groups by the pair of offset and
#' population index, which degenerates to one group per observation as soon as
#' `r` varies by row.
#'
#' @noRd
.vmgeom_log_dens <- function(m, lim_magn, r, perception_fun) {
    n <- max(length(m), length(lim_magn), length(r))
    m <- rep(m, length.out = n)
    tab <- .vmgeom_table(
        rep(lim_magn, length.out = n),
        rep(r, length.out = n),
        perception_fun
    )
    row_of <- tab$row_of

    # Magnitudes are counted in whole steps below the rounded limiting
    # magnitude, matching `.vmgeom_std()`. Anything fainter than the limiting
    # magnitude cannot be observed at all.
    magn_of <- tab$lim_magn_round - m
    logdens <- rep(-Inf, length(magn_of))

    observable <- magn_of > -1L

    tabulated <- observable & magn_of < tab$magn_max
    if (any(tabulated)) {
        cell <- cbind(row_of[tabulated], magn_of[tabulated] + 1L)
        logdens[tabulated] <- log(tab$visible[cell])
    }

    # Beyond the tabulated magnitudes the perception probability is 1, leaving
    # the plain geometric distribution.
    beyond <- magn_of >= tab$magn_max
    if (any(beyond)) {
        logdens[beyond] <- stats::dgeom(
            magn_of[beyond], tab$rate_norm[row_of[beyond]],
            log = TRUE
        )
    }

    logdens[observable] <- logdens[observable] -
        log(tab$norm[row_of[observable]])

    logdens
}

#' tabulate the visible part of the geometric model
#'
#' The observable magnitudes of one observation form a finite table, running
#' from the limiting magnitude up to the brightness at which the perception
#' probability reaches `1`. This helper builds that table once per distinct
#' combination of offset and population index, and the perception probabilities
#' once per distinct offset. Everything the moments and the density need is
#' derived from it.
#'
#' @noRd
.vmgeom_table <- function(lim_magn, r, perception_fun) {
    n <- max(length(lim_magn), length(r))
    lim_magn <- rep(lim_magn, length.out = n)
    r <- rep(r, length.out = n)

    # The distribution depends on `lim_magn` only through its fractional part,
    # so the perception probabilities can be shared by all observations having
    # the same offset. Standardizing at `m = 0` yields the rounded limiting
    # magnitude together with that offset, and keeps the treatment of the
    # half-magnitude boundaries identical to the distribution functions.
    std <- .vmgeom_std(0.0, lim_magn)
    lim_magn_round <- std$m
    offset <- std$offset

    # Observations sharing both offset and population index share their entire
    # table, so each combination is tabulated only once. `row_of` maps every
    # observation back to its row. The key is built numerically from the
    # positions within the distinct values, because pasting one string per
    # observation dominates the runtime of the whole fit. Unlike `paste()`,
    # which rounds to 15 significant digits, this compares bitwise and may
    # therefore split values agreeing to that many digits into separate
    # groups. That is on the safe side: every group is still evaluated
    # exactly, only a redundant tabulation may result.
    offset_key <- match(offset, unique(offset))
    key <- offset_key + (match(r, unique(r)) - 1L) * max(offset_key)
    unique_idx <- !duplicated(key)
    row_of <- match(key, key[unique_idx])
    offset <- offset[unique_idx]
    r <- r[unique_idx]

    # Magnitudes are counted in whole steps below the rounded limiting
    # magnitude: step 0 is the faintest still observable magnitude, `magn_max`
    # the brightest one whose perception probability is below 1. Each row of
    # `magn_of` holds those steps for one combination.
    magn_max <- 15L
    magn_seq <- as.numeric(seq(0L, magn_max))
    magn_of <- matrix(
        magn_seq,
        nrow = length(r), ncol = magn_max + 1L, byrow = TRUE
    )

    # Lowering the limiting magnitude by one step reduces the observable rate
    # by the factor `r_inv`, so `r_inv^magn_of` is the share of meteors left
    # after that many steps. `rate_norm` turns those shares into probabilities,
    # being the reciprocal of their sum over all steps.
    r_inv <- 1.0 / r
    rate_norm <- 1.0 - r_inv

    # The perception probabilities depend on the offset alone. Evaluating them
    # once per distinct offset, rather than once per combination, is what keeps
    # the family affordable when `r` differs for every observation.
    unique_offset <- unique(offset)
    perception <- vapply(
        unique_offset,
        \(offset) perception_fun(magn_seq + offset),
        numeric(magn_max + 1L)
    )
    perception <- t(perception)[match(offset, unique_offset), , drop = FALSE]

    # Probability of observing each tabulated magnitude, up to the overall
    # normalization. Beyond `magn_max` the perception probability is 1, so the
    # remaining magnitudes contribute the plain geometric tail `beyond`.
    visible <- rate_norm * r_inv^magn_of * perception
    beyond <- r_inv^(magn_max + 1L)

    list(
        lim_magn_round = lim_magn_round,
        row_of = row_of,
        magn_max = magn_max,
        magn_of = magn_of,
        rate_norm = rate_norm,
        r_inv = r_inv,
        visible = visible,
        beyond = beyond,
        norm = rowSums(visible) + beyond
    )
}

#' attainable range of the mean
#'
#' @noRd
.vmgeom_mu_range <- function(lim_magn, perception_fun, r_range) {
    list(
        lower = .vmgeom_moments(lim_magn, r_range[1L], perception_fun)$mu,
        upper = .vmgeom_moments(lim_magn, r_range[2L], perception_fun)$mu
    )
}

#' population index from the mean
#'
#' Inverts `.vmgeom_moments()` numerically. The mean is strictly increasing in
#' `r`, so a unique solution exists within the bracketing interval.
#'
#' @noRd
.vmgeom_r_from_mu <- function(mu, lim_magn, perception_fun, r_range) {
    n <- max(length(mu), length(lim_magn))
    mu <- rep(mu, length.out = n)
    lim_magn <- rep(lim_magn, length.out = n)
    r_min <- r_range[1L]
    r_max <- r_range[2L]

    # Each distinct pair of mean and limiting magnitude is inverted only once.
    # The key is numeric for the same reason as in `.vmgeom_table()`.
    mu_key <- match(mu, unique(mu))
    key <- mu_key + (match(lim_magn, unique(lim_magn)) - 1L) * max(mu_key)
    idx <- !duplicated(key)
    r <- mapply(\(mu, lim_magn) {
        f <- \(r) .vmgeom_moments(lim_magn, r, perception_fun)$mu - mu

        # Clamp to the attainable range. Magnitudes at the perception limit
        # have no finite population index.
        if (!is.finite(mu) || f(r_min) >= 0.0) {
            return(r_min)
        }
        if (f(r_max) <= 0.0) {
            return(r_max)
        }

        stats::uniroot(f, c(r_min, r_max), tol = .Machine$double.eps^0.5)$root
    }, mu[idx], lim_magn[idx])

    as.numeric(r[match(key, key[idx])])
}
