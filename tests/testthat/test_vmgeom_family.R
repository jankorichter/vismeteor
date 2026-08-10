test_that("vmgeom_glm reproduces the maximum likelihood estimate", {
    set.seed(42L)
    r <- 2.3
    lim_magn <- c(5.6, 6.0, 6.3, 6.5)
    n <- 800L

    magn <- unlist(lapply(lim_magn, \(lim_magn) rvmgeom(n, lim_magn, r)))
    lim_magn <- rep(lim_magn, each = n)
    obs <- data.frame(magn = magn, lim_magn = lim_magn)

    fit <- vmgeom_glm(cbind(magn, lim_magn) ~ 1, data = obs)
    r_glm <- as.numeric(predict(fit, obs[1L, ], type = "r"))

    ll <- \(r) -sum(dvmgeom(magn, lim_magn, r, log = TRUE))
    r_mle <- stats::optim(
        2.0, ll,
        method = "Brent", lower = 1.2, upper = 4.0
    )$par

    expect_equal(r_glm, r_mle, tolerance = 1e-4)
})

test_that("vmgeom_glm recovers known coefficients", {
    set.seed(11L)
    b0 <- 0.0
    b1 <- 0.4

    x <- rep(seq(-1.5, 1.5, length.out = 8L), each = 150L)
    lim_magn <- rep(c(5.8, 6.1, 6.4), length.out = length(x))
    magn <- mapply(
        \(lim_magn, r) rvmgeom(1L, lim_magn, r),
        lim_magn, 1.0 + exp(-(b0 + b1 * x))
    )
    obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)

    fit <- vmgeom_glm(cbind(magn, lim_magn) ~ x, data = obs)
    est <- summary(fit)$coefficients

    # Both coefficients must lie within three standard errors of the truth.
    expect_lt(abs(est[1L, 1L] - b0), 3.0 * est[1L, 2L])
    expect_lt(abs(est[2L, 1L] - b1), 3.0 * est[2L, 2L])
    expect_true(fit$converged)
})

test_that("vmgeom_glm keeps the limiting magnitudes aligned", {
    set.seed(7L)
    r <- 2.2
    lim_magn <- rep(c(5.7, 6.2, 6.5), each = 200L)
    magn <- mapply(\(lim_magn) rvmgeom(1L, lim_magn, r), lim_magn)
    obs <- data.frame(magn = magn, lim_magn = lim_magn)

    fit_all <- vmgeom_glm(cbind(magn, lim_magn) ~ 1, data = obs)

    # Dropping rows through na.action must drop the limiting magnitudes too.
    obs_na <- obs
    obs_na$magn[c(5L, 100L, 400L)] <- NA_integer_
    fit_na <- vmgeom_glm(cbind(magn, lim_magn) ~ 1, data = obs_na)
    expect_equal(stats::nobs(fit_na), stats::nobs(fit_all) - 3L)

    # The same fit, once via subset and once on pre-subsetted data.
    idx <- seq_len(300L)
    fit_subset <- vmgeom_glm(cbind(magn, lim_magn) ~ 1, data = obs, subset = idx)
    fit_manual <- vmgeom_glm(cbind(magn, lim_magn) ~ 1, data = obs[idx, ])
    expect_equal(unname(coef(fit_subset)), unname(coef(fit_manual)))
})

test_that("vmgeom_glm converges for a small population index", {
    set.seed(3L)
    r <- 1.2
    lim_magn <- rep(c(5.8, 6.1, 6.4), each = 500L)
    magn <- mapply(\(lim_magn) rvmgeom(1L, lim_magn, r), lim_magn)
    obs <- data.frame(magn = magn, lim_magn = lim_magn)

    # A log link would allow r < 1 here, which the geometric model forbids.
    fit <- vmgeom_glm(cbind(magn, lim_magn) ~ 1, data = obs)
    r_est <- as.numeric(predict(fit, obs[1L, ], type = "r"))

    expect_true(fit$converged)
    expect_true(r_est > 1.0)
    expect_equal(r_est, r, tolerance = 0.1)
})

test_that("vmgeom_glm accepts weights as replication", {
    set.seed(5L)
    obs <- data.frame(
        magn = c(1L, 2L, 3L, 4L, 2L, 3L),
        lim_magn = c(6.1, 6.1, 6.1, 6.1, 5.8, 5.8),
        Freq = c(3L, 5L, 2L, 4L, 6L, 1L)
    )
    replicated <- obs[rep(seq_len(nrow(obs)), obs$Freq), ]

    fit_weighted <- vmgeom_glm(
        cbind(magn, lim_magn) ~ 1,
        data = obs, weights = Freq
    )
    fit_replicated <- vmgeom_glm(cbind(magn, lim_magn) ~ 1, data = replicated)

    expect_equal(unname(coef(fit_weighted)), unname(coef(fit_replicated)))
    expect_equal(deviance(fit_weighted), deviance(fit_replicated))
    expect_equal(AIC(fit_weighted), AIC(fit_replicated))

    # The dispersion is fixed, so aggregating the observations must not change
    # the standard errors even though `df.residual` differs.
    expect_equal(
        unname(summary(fit_weighted)$coefficients[, 2L]),
        unname(summary(fit_replicated)$coefficients[, 2L])
    )

    newdata <- obs[1L, ]
    expect_equal(
        predict(fit_weighted, newdata, se.fit = TRUE),
        predict(fit_replicated, newdata, se.fit = TRUE)
    )
})

test_that("predict returns one value per row of newdata", {
    set.seed(11L)
    x <- rep(seq(-1.5, 1.5, length.out = 8L), each = 60L)
    lim_magn <- rep(c(5.8, 6.1, 6.4), length.out = length(x))
    magn <- mapply(
        \(lim_magn, r) rvmgeom(1L, lim_magn, r),
        lim_magn, 1.0 + exp(-(0.2 + 0.4 * x))
    )
    obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)

    fit <- vmgeom_glm(cbind(magn, lim_magn) ~ x, data = obs)
    newdata <- data.frame(x = c(-1.0, 0.0, 1.0))

    # Every scale is derived from the linear predictor, so the result must
    # follow `newdata` rather than the observations of the fit.
    for (type in c("r", "inv_r", "log_r", "link")) {
        expect_length(predict(fit, newdata, type = type), nrow(newdata))
    }

    # The scales must be consistent with each other.
    eta <- predict(fit, newdata, type = "link")
    r <- predict(fit, newdata, type = "r")
    expect_equal(r, 1.0 + exp(-eta))
    expect_equal(predict(fit, newdata, type = "inv_r"), 1.0 / r)
    expect_equal(predict(fit, newdata, type = "log_r"), log(r))

    # Confidence limits bracket the estimate and respect r > 1.
    pred <- predict(fit, newdata, se.fit = TRUE)
    expect_true(all(pred$lwr < pred$fit))
    expect_true(all(pred$upr > pred$fit))
    expect_true(all(pred$lwr > 1.0))

    # The limits are the transformed limits of the linear predictor, not a
    # symmetric interval around the estimate.
    pred <- predict(fit, newdata, type = "r", se.fit = TRUE)
    eta <- predict(fit, newdata, type = "link")
    se_eta <- predict(fit, newdata, type = "link", se.fit = TRUE)$se.fit
    z <- stats::qnorm(0.975)
    expect_equal(pred$lwr, 1.0 + exp(-(eta + z * se_eta)))
    expect_equal(pred$upr, 1.0 + exp(-(eta - z * se_eta)))
    expect_false(isTRUE(all.equal(pred$upr - pred$fit, pred$fit - pred$lwr)))

    # se.fit is the first-order delta-method value.
    expect_equal(pred$se.fit, abs(-exp(-eta)) * se_eta)

    # The expected magnitude is not among the offered scales: it depends on
    # the limiting magnitudes, which `newdata` does not carry.
    expect_error(predict(fit, newdata, type = "mu"), "arg")

    # The family is internal, but a fit built with it directly must behave the
    # same. predict() then returns the linear predictor, from which the
    # population index follows. Only type = "response" is refused, because it
    # would evaluate the mean against the limiting magnitudes of the fit
    # instead of those belonging to the new data.
    plain <- stats::glm(
        cbind(magn, lim_magn) ~ x,
        family = .vmgeom_family(), data = obs
    )
    expect_equal(coef(plain), coef(fit))
    expect_equal(
        1.0 + exp(-stats::predict(plain, newdata)),
        predict(fit, newdata)
    )
    expect_error(
        stats::predict(plain, newdata, type = "response"),
        "vmgeom_glm"
    )
    expect_length(fitted(plain), nrow(obs))
})

test_that("vmgeom_glm supports the standard model generics", {
    set.seed(23L)
    x <- rep(seq(-1.0, 1.0, length.out = 6L), each = 80L)
    lim_magn <- rep(c(5.9, 6.3), length.out = length(x))
    magn <- mapply(
        \(lim_magn, r) rvmgeom(1L, lim_magn, r),
        lim_magn, 1.0 + exp(-(0.1 + 0.3 * x))
    )
    obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)

    fit <- vmgeom_glm(cbind(magn, lim_magn) ~ x, data = obs)

    expect_s3_class(fit, "vmgeom_glm")
    expect_s3_class(fit, "glm")
    expect_s3_class(summary(fit), "summary.glm")
    expect_true(is.finite(AIC(fit)))
    expect_true(is.finite(BIC(fit)))

    # update() must be able to re-evaluate the recorded call.
    refit <- update(fit, . ~ 1)
    expect_s3_class(refit, "vmgeom_glm")
    expect_length(coef(refit), 1L)

    # fitted() refers to the observations of the fit, whose limiting
    # magnitudes the family still holds, and therefore stays correct.
    expect_length(fitted(fit), nrow(obs))
    expect_true(all(is.finite(fitted(fit))))
})

test_that("vmgeom moments match a direct summation", {
    reference <- \(lim_magn, r) {
        m <- seq(floor(lim_magn + 0.5), -200.0, -1.0)
        d <- dvmgeom(m, lim_magn, r)
        mu <- sum(m * d)

        c(mu, sum((m - mu)^2 * d))
    }

    model <- expand.grid(
        lim_magn = c(4.5, 5.5, 5.52, 5.8, 6.0, 6.1, 6.3, 6.48, 6.5, 7.5),
        r = c(1.2, 1.8, 2.4, 3.2, 5.0)
    )

    for (i in seq_len(nrow(model))) {
        moments <- .vmgeom_moments(model$lim_magn[i], model$r[i], vmperception)
        expected <- reference(model$lim_magn[i], model$r[i])

        expect_equal(moments$mu, expected[1L], tolerance = 1e-8)
        expect_equal(moments$var, expected[2L], tolerance = 1e-8)
    }
})

test_that("vmgeom_glm reuses the population index without changing the fit", {
    set.seed(31L)
    x <- rep(seq(-1.2, 1.2, length.out = 8L), each = 60L)
    lim_magn <- rep(c(5.8, 6.2, 6.4), length.out = length(x))
    magn <- mapply(
        \(lim_magn, r) rvmgeom(1L, lim_magn, r),
        lim_magn, 1.0 + exp(-(0.1 + 0.35 * x))
    )
    obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)

    fit <- vmgeom_glm(cbind(magn, lim_magn) ~ x, data = obs)

    # `linkinv()` remembers which population index produced the current mean,
    # so that `variance()` need not invert it numerically. Discarding that
    # memory after every call forces the inversion and must not change the fit.
    family <- .vmgeom_family()
    env <- environment(family$linkinv)$env
    linkinv <- family$linkinv
    family$linkinv <- \(eta) {
        mu <- linkinv(eta)
        env$mu <- NULL
        env$r <- NULL

        mu
    }
    fit_uncached <- stats::glm(
        cbind(magn, lim_magn) ~ x,
        family = family, data = obs
    )

    expect_equal(coef(fit), coef(fit_uncached), tolerance = 1e-8)
})

test_that("vmgeom log density matches dvmgeom", {
    # Half-magnitude limits are the delicate cases: `round()` breaks ties to
    # the nearest even number, so 4.5 and 6.5 round down while 5.5 and 7.5
    # round up. The offset is then +/-0.5, and the brightest magnitude falls
    # outside the tabulated ones.
    model <- expand.grid(
        lim_magn = c(4.5, 5.5, 5.52, 5.8, 6.0, 6.1, 6.3, 6.48, 6.5, 7.5),
        r = c(1.2, 1.8, 2.4, 3.2, 5.0)
    )
    for (i in seq_len(nrow(model))) {
        magn <- seq(floor(model$lim_magn[i] + 0.5), -25.0)
        expect_equal(
            .vmgeom_log_dens(magn, model$lim_magn[i], model$r[i], vmperception),
            dvmgeom(magn, model$lim_magn[i], model$r[i], log = TRUE)
        )
    }

    # A population index that differs for every observation is the case the
    # family relies on, and the one dvmgeom() itself handles slowly.
    set.seed(3L)
    n <- 500L
    lim_magn <- rep(c(5.7, 6.0, 6.3, 6.5, 6.1), length.out = n)
    r <- 1.0 + exp(-rnorm(n, 0.2, 0.3))
    magn <- mapply(\(lim_magn, r) rvmgeom(1L, lim_magn, r), lim_magn, r)

    expect_equal(
        .vmgeom_log_dens(magn, lim_magn, r, vmperception),
        dvmgeom(magn, lim_magn, r, log = TRUE)
    )

    # Magnitudes outside the support are impossible, not merely unlikely.
    expect_equal(.vmgeom_log_dens(-Inf, 6.0, 2.4, vmperception), -Inf)
    expect_equal(.vmgeom_log_dens(8.0, 6.0, 2.4, vmperception), -Inf)
})

test_that("vmgeom_glm bounds the population index", {
    # The mean diverges as r -> 1, so the moments must stay finite at the
    # lower end of the bracketing range.
    moments <- .vmgeom_moments(6.0, 1.01, vmperception)
    expect_true(is.finite(moments$mu))
    expect_true(is.finite(moments$var))

    expect_error(.vmgeom_family(r_range = c(1.0, 20.0)), "greater than 1")
    expect_error(.vmgeom_family(r_range = c(3.0, 2.0)), "increasing")
    expect_error(.vmgeom_family(r_range = 2.0), "two finite values")

    # Degenerate data, where every meteor is as faint as the limiting
    # magnitude, imply a population index below the bracketing range and must
    # fail with a comprehensible message rather than an internal error.
    obs <- data.frame(magn = rep(6L, 50L), lim_magn = rep(6.0, 50L))
    expect_error(
        suppressWarnings(vmgeom_glm(cbind(magn, lim_magn) ~ 1, data = obs)),
        "coefficients"
    )
})

test_that("vmgeom_glm passes its arguments to the model family", {
    set.seed(9L)
    lim_magn <- rep(c(5.8, 6.1, 6.4), each = 200L)
    magn <- mapply(\(lim_magn) rvmgeom(1L, lim_magn, 2.4), lim_magn)
    obs <- data.frame(magn = magn, lim_magn = lim_magn)

    # `r_range` and `perception_fun` are the only arguments the model family
    # takes, and the wrapper is the only way to reach them.
    expect_equal(
        coef(vmgeom_glm(cbind(magn, lim_magn) ~ 1, data = obs, r_range = c(1.2, 5.0))),
        coef(stats::glm(
            cbind(magn, lim_magn) ~ 1,
            family = .vmgeom_family(r_range = c(1.2, 5.0)), data = obs
        ))
    )

    perception_fun <- \(dm) ifelse(dm > -0.5, 1.0 - exp(-0.005 * (dm + 0.5)^3), 0.0)
    expect_equal(
        coef(vmgeom_glm(
            cbind(magn, lim_magn) ~ 1,
            data = obs, perception_fun = perception_fun
        )),
        coef(stats::glm(
            cbind(magn, lim_magn) ~ 1,
            family = .vmgeom_family(perception_fun = perception_fun), data = obs
        ))
    )
})

test_that("vmgeom_glm rejects an invalid response", {
    obs <- data.frame(
        magn = c(1L, 2L, 3L, 4L),
        lim_magn = c(6.1, 6.1, 5.8, 5.8)
    )

    expect_error(
        vmgeom_glm(magn ~ 1, data = obs),
        "two column matrix"
    )
})
