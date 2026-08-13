#' @name vmideal_vst
#' @aliases vmideal_vst_from_magn
#' @aliases vmideal_vst_to_psi
#' @title Variance-stabilizing Transformation for the Ideal Distribution of Visual Meteor Magnitudes
#'
#' @description
#' Applies a variance-stabilizing transformation to meteor magnitudes
#' under the assumption of the ideal magnitude distribution.
#'
#' @param m integer; the meteor magnitude.
#' @param lm numeric; limiting magnitude.
#' @param tm numeric; transformed magnitude.
#' @param deriv_degree integer; the degree of the derivative at `tm` to return
#'   instead of `psi`. Must be `0`, `1` or `2`.
#'
#' @details
#' Many linear models require the variance of visual meteor magnitudes to be
#' homoscedastic. The function `vmideal_vst_from_magn` applies a transformation
#' that produces homoscedastic distributions of visual meteor magnitudes if the
#' underlying magnitudes follow the ideal magnitude distribution.
#' In this sense, the transformation acts as a normalization of meteor magnitudes
#' and yields a variance close to `1.0`.
#'
#' The ideal distribution of visual meteor magnitudes
#' depends on the [parameter psi][vismeteor::vmideal] and the limiting magnitude `lm`,
#' resulting in a two-parameter distribution. Without detection probabilities,
#' the magnitude distribution reduces to a pure [ideal magnitude distribution][vismeteor::mideal],
#' which depends only on the parameter `psi`. Since the
#' limiting magnitude `lm` is a fixed parameter and never estimated
#' statistically, the magnitudes can be transformed such that, for example,
#' the mean of the transformed magnitudes directly provides an estimate of `psi`
#' using the function `vmideal_vst_to_psi`.
#'
#' What the transformation resolves is the difference between `psi` and `lm`,
#' so its range follows the limiting magnitude: `psi` is recovered from about
#' \eqn{\texttt{lm} - 16.5} up to about \eqn{\texttt{lm} + 3}. At a limiting
#' magnitude of around `6.0` this amounts to \eqn{-10 \le \texttt{psi} \le 9}.
#' Within that range the estimate is accurate to a few hundredths of a
#' magnitude, degrading towards the upper end, where the mean of the
#' distribution begins to flatten. Outside it `vmideal_vst_to_psi` returns `NA`
#' or `Inf` rather than an unreliable value; see below.
#' The numerical form of the transformation is version-specific and may change
#' substantially in future releases. Do not rely on equality of transformed
#' values across package versions.
#'
#' @return
#' * `vmideal_vst_from_magn`: a numeric value, the transformed meteor magnitude.
#' * `vmideal_vst_to_psi`: a numeric value of the parameter `psi`, derived from
#'   the mean of `tm`.
#' The argument `deriv_degree` can be used to apply the delta method.
#'
#' `vmideal_vst_to_psi` returns `Inf` for a `tm` below `0.02`. A vanishing mean
#' of the transformed magnitudes is what an unbounded `psi` produces, and the
#' magnitudes are then geometric with \eqn{r = 10^{0.4}}, which
#' [dvmideal][vismeteor::vmideal] evaluates for an infinite `psi`. The
#' derivatives do not exist there and are `NA`. Above the range the
#' transformation covers, at a `tm` beyond `8.22`, the estimate is `NA` as
#' well: a `psi` that far below the limiting magnitude is not a limit but
#' simply out of range.
#'
#' @note
#' The internal approximations used here are derived from the perception
#' probabilities produced by [vismeteor::vmperception].
#' For details on the derivation, see the script `inst/derivation/vmideal_vst.R` in the
#' package's source code.
#'
#' @seealso [vismeteor::vmideal] [vismeteor::mideal] [vismeteor::vmperception]
#' @examples
#' N <- 100
#' psi <- 5.0
#' limmag <- 6.3
#'
#' # Simulate magnitudes
#' m <- rvmideal(N, limmag, psi)
#'
#' # Variance-stabilizing transformation
#' tm <- vmideal_vst_from_magn(m, limmag)
#' tm_mean <- mean(tm)
#' tm_var <- var(tm)
#'
#' # Estimator for psi from the transformed mean
#' psi_hat <- vmideal_vst_to_psi(tm_mean, limmag)
#'
#' # Derivative d(psi)/d(tm) at tm_mean (needed for the delta method)
#' dpsi_dtm <- vmideal_vst_to_psi(tm_mean, limmag, deriv_degree = 1L)
#'
#' # Variance of the sample mean of tm
#' var_tm.mean <- tm_var / N
#'
#' # Delta method: variance and standard error of psi_hat
#' var_psi.hat <- (dpsi_dtm^2) * var_tm.mean
#' se_psi.hat <- sqrt(var_psi.hat)
#'
#' # Results
#' print(psi_hat)
#' print(se_psi.hat)

#' @rdname vmideal_vst
#' @export
vmideal_vst_from_magn <- function(m, lm) {
    offset <- lm - round(lm)
    if (1L == length(lm)) {
        offset <- rep(offset, length(m))
    }

    arg_data <- data.frame(
        x = lm - m,
        offset = offset
    )

    data_f <- as.factor(arg_data$offset)
    data_s <- split(arg_data, data_f)
    y <- lapply(data_s, \(data) {
        x <- data$x
        offset <- data$offset[1]
        params <- .vmideal_vst_from_magn_params
        sy <- c(
            params$p1_fun(offset),
            params$p2_fun(offset),
            params$p3_fun(offset),
            params$p4_fun(offset),
            params$p5_fun(offset),
            params$p6_fun(offset)
        )
        sx <- params$sx
        f_spline <- stats::splinefun(sx, sy)
        y <- rep(NA_real_, length(x))

        idx <- x < sx[1]
        if (any(idx)) {
            y[idx] <- f_spline(sx[1], deriv = 1) * (x[idx] - sx[1]) + sy[1]
        }

        idx <- x > sx[length(sx)]
        if (any(idx)) {
            y[idx] <- f_spline(sx[length(sx)], deriv = 1) * (x[idx] - sx[length(sx)]) + sy[length(sy)]
        }

        idx <- x >= sx[1] & x <= sx[length(sx)]
        if (any(idx)) {
            y[idx] <- f_spline(x[idx])
        }

        y - params$p0_fun(offset)
    })

    unsplit(y, data_f)
}

#' @rdname vmideal_vst
#' @export
vmideal_vst_to_psi <- function(tm, lm, deriv_degree = 0L) {
    poly_coef0 <- c(
        -2.97081142349313, -2.72861507781908, -0.683119285080101, -0.197203754814882,
        -0.0707733355336266, -0.0267800890185785, -0.00481827712622998,
        -1.03426983830783e-05, 5.99984544744159e-05
    )
    names(poly_coef0) <- seq_along(poly_coef0) - 1 # exponents

    # A vanishing mean of the transformed magnitudes is what an unbounded `psi`
    # produces, so it is reported as such rather than as a failure of the
    # transformation. Below the threshold the polynomial has left the range it
    # was fitted on, while the mean is already so flat that no finite `psi` is
    # resolved: the ideal distribution equals its geometric limit from about
    # `lm + 10` onwards, and the mean stops distinguishing values well before
    # that.
    #
    # The threshold is the same one that bounds the polynomial below, so that
    # every `tm` too small for it is answered with the limit rather than with a
    # missing value.
    unbound <- !is.na(tm) & tm < 0.02

    # x min 0.02 (psi approx 9 at limiting magnitude of 6.0)
    # x max 8.22(psi approx -10 at limiting maginitde of 6.5)
    # The lower end is excluded here as well, so that a negative `tm` does not
    # warn about a `NaN` that the limit replaces further down.
    tm[tm < 0.02 | tm > 8.22] <- NA_real_
    y <- log(tm)

    if (deriv_degree > 0L) {
        poly_coef1 <- .f_polynomial_coef(poly_coef0, deriv_degree = 1L)
    }
    if (deriv_degree > 1L) {
        poly_coef2 <- .f_polynomial_coef(poly_coef1, deriv_degree = 1L)
    }
    if (deriv_degree > 2L) {
        stop(paste("deriv_degree", deriv_degree, "not implemented!"))
    }

    psi <- if (2L == deriv_degree) {
        (.f_polynomial(y, poly_coef2) - .f_polynomial(y, poly_coef1)) / tm^2
    } else if (1L == deriv_degree) {
        .f_polynomial(y, poly_coef1) / tm
    } else {
        lm + .f_polynomial(y, poly_coef0)
    }

    # `psi` is unbounded rather than merely large there, so its derivatives do
    # not exist as finite numbers. They serve the delta method, which needs a
    # finite `psi` to be applied to in the first place, and are reported as
    # missing instead of being continued.
    psi[unbound] <- if (0L == deriv_degree) Inf else NA_real_

    psi
}

#' @keywords internal
.vmideal_vst_from_magn_params <- (function() {
    param_df <- data.frame(
        "p1" = c(
            0.379578706193683, 0.380865978506213,
            0.383646581863156, 0.386594955357005, 0.384516871380943, 0.378119364136071,
            0.39783016804395, 0.373393379287508, 0.366393502192895, 0.345593857555874,
            0.367086720870886, 0.359195376100079, 0.347082537633701, 0.312420722892609,
            0.312232925285191
        ), "p2" = c(
            0.758955184180252, 0.757434042449928,
            0.756545241530695, 0.754597324053272, 0.746096530130275, 0.736273688899513,
            0.754472950700777, 0.729839979801878, 0.724090251762525, 0.706272163875503,
            0.732424006749156, 0.730856947686901, 0.722541149425633, 0.690367122987635,
            0.691926213129086
        ), "p3" = c(
            1.89943177726206, 1.89844069948433,
            1.89825430762927, 1.89725244849407, 1.88997593310439, 1.88072563899649,
            1.89899661183851, 1.87396958675163, 1.86749053248022, 1.84882849631795,
            1.87410789172307, 1.87172629242088, 1.86305301012132, 1.83068592572223,
            1.83212684285442
        ), "p4" = c(
            3.27014476551314, 3.26892787903128,
            3.26844009590108, 3.26702893292816, 3.25921320867758, 3.24969559483436,
            3.26789784387066, 3.24296948683872, 3.23673983810267, 3.21843295052407,
            3.24413734533724, 3.24221462073629, 3.23377186059092, 3.20154167652271,
            3.20307252532895
        ), "p5" = c(
            4.50424034207103, 4.50305721713432,
            4.50261471059499, 4.5012654394704, 4.49353145072154, 4.48405231361262,
            4.50225758259143, 4.47729738300568, 4.47099710887246, 4.45259464315401,
            4.47819774376931, 4.4761817895577, 4.46769920746914, 4.43544755919314,
            4.43696568503357
        ), "p6" = c(
            5.65036712145262, 5.64927790140628,
            5.64896064069461, 5.6477824229466, 5.64028041416121, 5.63092781192841,
            5.64918249953455, 5.62420382559607, 5.61781703615347, 5.5992741104315,
            5.62469988859743, 5.62248449312786, 5.61389842883615, 5.58158409117071,
            5.58306055049816
        ), "offset" = c(
            -0.5, -0.48, -0.45, -0.4, -0.3,
            -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.48, 0.5
        ), "intercept" = c(
            1.4811934820941,
            1.47927325731333, 1.47783119789222, 1.47503283992891, 1.46516880013637,
            1.45440969975367, 1.47200033556686, 1.44700869631271, 1.44106727057982,
            1.42347451350523, 1.45029075335968, 1.44986985073928, 1.44229503571524,
            1.41062845239295, 1.41258415137388
        )
    )

    list(
        sx = c(1.0, 2.0, 4.0, 6.0, 8.0, 10.0),
        p0_fun = stats::approxfun(param_df$offset, param_df[["intercept"]]),
        p1_fun = stats::approxfun(param_df$offset, param_df[["p1"]]),
        p2_fun = stats::approxfun(param_df$offset, param_df[["p2"]]),
        p3_fun = stats::approxfun(param_df$offset, param_df[["p3"]]),
        p4_fun = stats::approxfun(param_df$offset, param_df[["p4"]]),
        p5_fun = stats::approxfun(param_df$offset, param_df[["p5"]]),
        p6_fun = stats::approxfun(param_df$offset, param_df[["p6"]])
    )
})()
