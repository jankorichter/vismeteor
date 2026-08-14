#' magnitudes of the vignette excerpt
#'
#' The same subset `vignette("vmgeom")` uses, so that the values below can be
#' compared with the ones shown there.
#'
#' @noRd
vmgeom_vst_lm_data <- function() {
    observations <- vismeteor::PER_2015_magn$observations
    observations <- subset(
        observations,
        !is.na(observations$lim_magn) &
            observations$sl_start > 135.81 &
            observations$sl_end < 135.87,
        select = c("magn_id", "lim_magn")
    )

    magnitudes <- merge(
        observations,
        as.data.frame(vismeteor::PER_2015_magn$magnitudes),
        by = "magn_id"
    )
    magnitudes$magn <- as.integer(as.character(magnitudes$magn))

    subset(magnitudes, (magnitudes$lim_magn - magnitudes$magn) > -0.5)
}

test_that("vmgeom_vst_lm fits the transformed magnitudes", {
    magnitudes <- vmgeom_vst_lm_data()

    fit <- vismeteor::vmgeom_vst_lm(
        cbind(magn, lim_magn) ~ 1,
        data = magnitudes,
        weights = magnitudes$Freq
    )

    expect_s3_class(fit, "vmgeom_vst_lm")
    expect_s3_class(fit, "lm")

    # the fit is on the transformed scale
    tm <- vismeteor::vmgeom_vst_from_magn(magnitudes$magn, magnitudes$lim_magn)
    expected <- stats::lm(tm ~ 1, weights = magnitudes$Freq)
    expect_equal(unname(stats::coef(fit)), unname(stats::coef(expected)))

    # the response must carry the limiting magnitude
    expect_error(
        vismeteor::vmgeom_vst_lm(magn ~ 1, data = magnitudes),
        "two-column matrix"
    )

    # `subset` and `na.action` reach the rows before the transformation does
    fit <- vismeteor::vmgeom_vst_lm(
        cbind(magn, lim_magn) ~ 1,
        data = magnitudes,
        weights = magnitudes$Freq,
        subset = magnitudes$Freq > 0.0
    )
    expect_length(stats::residuals(fit), sum(magnitudes$Freq > 0.0))

    incomplete <- magnitudes
    incomplete$magn[1L:3L] <- NA_integer_
    fit <- vismeteor::vmgeom_vst_lm(
        cbind(magn, lim_magn) ~ 1,
        data = incomplete,
        weights = incomplete$Freq
    )
    expect_length(stats::residuals(fit), nrow(incomplete) - 3L)

    # the formula may be held in a variable, whose environment is the caller's
    formula <- cbind(magn, lim_magn) ~ 1
    fit <- vismeteor::vmgeom_vst_lm(
        formula,
        data = magnitudes,
        weights = magnitudes$Freq
    )
    expect_equal(unname(stats::coef(fit)), 2.482395, tolerance = 1e-06)
})

test_that("predict.vmgeom_vst_lm rescales to the number of meteors", {
    magnitudes <- vmgeom_vst_lm_data()
    n <- sum(magnitudes$Freq)

    fit <- vismeteor::vmgeom_vst_lm(
        cbind(magn, lim_magn) ~ 1,
        data = magnitudes,
        weights = magnitudes$Freq
    )
    newdata <- data.frame(row.names = "")

    pred <- stats::predict(fit, newdata, type = "tm", se.fit = TRUE)
    expect_equal(unname(pred$fit), 2.482395, tolerance = 1e-06)

    # the standard error refers to the meteors, not to the rows: `lm()` alone
    # would report 0.1171802 here
    expect_equal(unname(pred$se.fit), 0.08199124, tolerance = 1e-06)

    tm <- vismeteor::vmgeom_vst_from_magn(magnitudes$magn, magnitudes$lim_magn)
    tm_mean <- sum(magnitudes$Freq * tm) / n
    tm_var <- sum(magnitudes$Freq * (tm - tm_mean)^2) / (n - 1)
    expect_equal(unname(pred$se.fit), sqrt(tm_var / n))

    # without weights the scale is the one `lm()` reports itself
    unweighted <- vismeteor::vmgeom_vst_lm(
        cbind(magn, lim_magn) ~ 1,
        data = magnitudes
    )
    pred <- stats::predict(unweighted, newdata, type = "tm", se.fit = TRUE)
    expect_equal(
        unname(pred$se.fit),
        summary(unweighted)$sigma / sqrt(length(stats::residuals(unweighted)))
    )

    # `scale` and `df` are determined by the model
    expect_error(stats::predict(fit, newdata, scale = 1.0), "must not be given")
    expect_error(stats::predict(fit, newdata, df = 10L), "must not be given")
})

test_that("predict.vmgeom_vst_lm returns the population index", {
    magnitudes <- vmgeom_vst_lm_data()

    fit <- vismeteor::vmgeom_vst_lm(
        cbind(magn, lim_magn) ~ 1,
        data = magnitudes,
        weights = magnitudes$Freq
    )
    newdata <- data.frame(row.names = "")

    r <- stats::predict(fit, newdata)
    expect_type(r, "double")
    expect_length(r, 1L)
    expect_equal(unname(r), 2.262779, tolerance = 1e-06)

    # the scales are consistent with each other
    expect_equal(
        unname(stats::predict(fit, newdata, type = "inv_r")),
        unname(1.0 / r)
    )
    expect_equal(
        unname(stats::predict(fit, newdata, type = "log_r")),
        unname(log(r))
    )

    # generic dispatch reaches this method rather than stats::predict.lm
    expect_equal(unname(stats::predict(fit, newdata, type = "r")), unname(r))
})

test_that("predict.vmgeom_vst_lm applies the bias correction", {
    magnitudes <- vmgeom_vst_lm_data()
    n <- sum(magnitudes$Freq)

    fit <- vismeteor::vmgeom_vst_lm(
        cbind(magn, lim_magn) ~ 1,
        data = magnitudes,
        weights = magnitudes$Freq
    )
    newdata <- data.frame(row.names = "")

    pred <- stats::predict(
        fit, newdata,
        se.fit = TRUE, bias_correction = TRUE
    )

    # the values `vignette("vmgeom")` reports
    expect_equal(unname(pred$fit), 2.2589086734166, tolerance = 1e-09)
    expect_equal(unname(pred$se.fit^2), 0.0100592098794041, tolerance = 1e-09)

    # the correction is exactly the second-order term of the delta method
    tm <- stats::predict(fit, newdata, type = "tm", se.fit = TRUE)
    tm_mean_var <- unname(tm$se.fit)^2
    uncorrected <- stats::predict(fit, newdata)
    expect_equal(
        unname(uncorrected - pred$fit),
        0.5 * vismeteor::vmgeom_vst_to_r(
            unname(tm$fit),
            deriv_degree = 2L
        ) * tm_mean_var
    )

    # by default no correction is applied
    expect_equal(
        unname(stats::predict(fit, newdata, bias_correction = FALSE)),
        unname(uncorrected)
    )
})

test_that("predict.vmgeom_vst_lm orders the confidence limits", {
    magnitudes <- vmgeom_vst_lm_data()

    fit <- vismeteor::vmgeom_vst_lm(
        cbind(magn, lim_magn) ~ 1,
        data = magnitudes,
        weights = magnitudes$Freq
    )
    newdata <- data.frame(row.names = "")

    # the limits rest on the t-distribution of the linear model; the normal
    # approximation would give 2.321695 and 2.643095 here
    limits <- stats::predict(fit, newdata, type = "tm", interval = "confidence")
    expect_equal(unname(limits[1L, "lwr"]), 2.319644, tolerance = 1e-06)
    expect_equal(unname(limits[1L, "upr"]), 2.645146, tolerance = 1e-06)

    # `r` and `log_r` decrease with `tm` while `inv_r` increases, so the limits
    # have to be reordered on some scales but not on others
    for (type in c("r", "inv_r", "log_r", "tm")) {
        limits <- stats::predict(
            fit, newdata,
            type = type, interval = "confidence"
        )
        expect_true(limits[1L, "lwr"] < limits[1L, "fit"], label = type)
        expect_true(limits[1L, "fit"] < limits[1L, "upr"], label = type)
    }

    limits <- stats::predict(fit, newdata, interval = "confidence")
    expect_equal(unname(limits[1L, "lwr"]), 2.0782, tolerance = 1e-05)
    expect_equal(unname(limits[1L, "upr"]), 2.478014, tolerance = 1e-05)
})

test_that("predict.vmgeom_vst_lm accepts covariates", {
    set.seed(3L)

    lim_magn <- rep(c(5.8, 6.1, 6.4), each = 200L)
    x <- rep(seq(-1.5, 1.5, length.out = 20L), each = 30L)
    r <- 1.0 + exp(-(0.2 + 0.4 * x))
    magn <- mapply(
        \(lim_magn, r) vismeteor::rvmgeom(1L, lim_magn, r),
        lim_magn, r
    )
    obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)

    fit <- vismeteor::vmgeom_vst_lm(cbind(magn, lim_magn) ~ x, data = obs)

    newdata <- data.frame(x = c(-1.0, 0.0, 1.0))
    pred <- stats::predict(fit, newdata)
    expect_length(pred, 3L)

    # `r` falls as the covariate grows, as it was simulated
    expect_true(all(diff(pred) < 0.0))

    pred <- stats::predict(fit, newdata, se.fit = TRUE)
    expect_length(pred$fit, 3L)
    expect_length(pred$se.fit, 3L)
    expect_true(all(pred$se.fit > 0.0))
})
