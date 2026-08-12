test_that("vmideal_glm approximates the maximum likelihood estimate", {
    set.seed(42L)
    psi <- 5.0
    lim_magn <- c(5.6, 6.0, 6.3, 6.5)
    n <- 800L

    magn <- unlist(lapply(lim_magn, \(lim_magn) rvmideal(n, lim_magn, psi)))
    lim_magn <- rep(lim_magn, each = n)
    obs <- data.frame(magn = magn, lim_magn = lim_magn)

    fit <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs)
    psi_glm <- as.numeric(predict(fit, obs[1L, ]))

    ll <- \(psi) -sum(dvmideal(magn, lim_magn, psi, log = TRUE))
    psi_mle <- stats::optim(
        4.0, ll,
        method = "Brent", lower = 0.0, upper = 10.0
    )$par

    # The mean is not a sufficient statistic for psi, so the fit is a
    # quasi-likelihood estimate rather than the maximum likelihood one. Both
    # are consistent, and they agree to within a few hundredths.
    expect_equal(psi_glm, psi_mle, tolerance = 0.1)
    expect_false(isTRUE(all.equal(psi_glm, psi_mle, tolerance = 1e-4)))
})

test_that("vmideal_glm recovers known coefficients", {
    set.seed(11L)
    b0 <- 4.0
    b1 <- 0.8

    x <- rep(seq(-1.5, 1.5, length.out = 8L), each = 250L)
    lim_magn <- rep(c(5.8, 6.1, 6.4), length.out = length(x))
    magn <- mapply(
        \(lim_magn, psi) rvmideal(1L, lim_magn, psi),
        lim_magn, b0 + b1 * x
    )
    obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)

    fit <- vmideal_glm(cbind(magn, lim_magn) ~ x, data = obs)
    est <- summary(fit)$coefficients

    # Both coefficients must lie within three standard errors of the truth.
    expect_lt(abs(est[1L, 1L] - b0), 3.0 * est[1L, 2L])
    expect_lt(abs(est[2L, 1L] - b1), 3.0 * est[2L, 2L])
    expect_true(fit$converged)
})

test_that("vmideal_glm keeps the limiting magnitudes aligned", {
    set.seed(7L)
    psi <- 4.5
    lim_magn <- rep(c(5.7, 6.2, 6.5), each = 200L)
    magn <- mapply(\(lim_magn) rvmideal(1L, lim_magn, psi), lim_magn)
    obs <- data.frame(magn = magn, lim_magn = lim_magn)

    fit_all <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs)

    # Dropping rows through na.action must drop the limiting magnitudes too.
    obs_na <- obs
    obs_na$magn[c(5L, 100L, 400L)] <- NA_integer_
    fit_na <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs_na)
    expect_equal(stats::nobs(fit_na), stats::nobs(fit_all) - 3L)

    # The same fit, once via subset and once on pre-subsetted data.
    idx <- seq_len(300L)
    fit_subset <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs, subset = idx)
    fit_manual <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs[idx, ])
    expect_equal(unname(coef(fit_subset)), unname(coef(fit_manual)))
})

test_that("vmideal_glm accepts weights as replication", {
    obs <- data.frame(
        magn = c(1L, 2L, 3L, 4L, 2L, 3L),
        lim_magn = c(6.1, 6.1, 6.1, 6.1, 5.8, 5.8),
        Freq = c(3L, 5L, 2L, 4L, 6L, 1L)
    )
    replicated <- obs[rep(seq_len(nrow(obs)), obs$Freq), ]

    fit_weighted <- vmideal_glm(
        cbind(magn, lim_magn) ~ 1,
        data = obs, weights = Freq
    )
    fit_replicated <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = replicated)

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
    x <- rep(seq(-1.5, 1.5, length.out = 8L), each = 100L)
    lim_magn <- rep(c(5.8, 6.1, 6.4), length.out = length(x))
    magn <- mapply(
        \(lim_magn, psi) rvmideal(1L, lim_magn, psi),
        lim_magn, 4.0 + 0.8 * x
    )
    obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)

    fit <- vmideal_glm(cbind(magn, lim_magn) ~ x, data = obs)
    newdata <- data.frame(x = c(-1.0, 0.0, 1.0))

    # Every scale is derived from the linear predictor, so the result must
    # follow `newdata` rather than the observations of the fit.
    for (type in c("psi", "link")) {
        expect_length(predict(fit, newdata, type = type), nrow(newdata))
    }

    # The link is the identity, so both scales coincide.
    expect_equal(predict(fit, newdata, type = "psi"), predict(fit, newdata, type = "link"))

    # The confidence limits are symmetric about the estimate.
    pred <- predict(fit, newdata, se.fit = TRUE)
    expect_named(pred, c("fit", "se.fit", "lwr", "upr"))
    expect_true(all(pred$lwr < pred$fit))
    expect_true(all(pred$upr > pred$fit))

    z <- stats::qnorm(0.975)
    expect_equal(pred$lwr, pred$fit - z * pred$se.fit)
    expect_equal(pred$upr, pred$fit + z * pred$se.fit)

    # A narrower level yields narrower limits.
    expect_true(all(predict(fit, newdata, se.fit = TRUE, level = 0.8)$upr < pred$upr))
    expect_error(predict(fit, newdata, se.fit = TRUE, level = 1.5), "level")

    # The expected magnitude is not among the offered scales: it depends on
    # the limiting magnitudes, which `newdata` does not carry.
    expect_error(predict(fit, newdata, type = "mu"), "arg")

    # The family is internal, but a fit built with it directly must behave the
    # same. Only type = "response" is refused, because it would evaluate the
    # mean against the limiting magnitudes of the fit instead of those
    # belonging to the new data.
    plain <- stats::glm(
        cbind(magn, lim_magn) ~ x,
        family = .vmideal_family(), data = obs
    )
    expect_equal(coef(plain), coef(fit))
    expect_equal(stats::predict(plain, newdata), predict(fit, newdata))
    expect_error(
        stats::predict(plain, newdata, type = "response"),
        "vmideal_glm"
    )
    expect_length(fitted(plain), nrow(obs))
})

test_that("vmideal_glm supports the standard model generics", {
    set.seed(23L)
    x <- rep(seq(-1.0, 1.0, length.out = 6L), each = 100L)
    lim_magn <- rep(c(5.9, 6.3), length.out = length(x))
    magn <- mapply(
        \(lim_magn, psi) rvmideal(1L, lim_magn, psi),
        lim_magn, 4.5 + 0.6 * x
    )
    obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)

    fit <- vmideal_glm(cbind(magn, lim_magn) ~ x, data = obs)

    expect_s3_class(fit, "vmideal_glm")
    expect_s3_class(fit, "glm")
    expect_s3_class(summary(fit), "summary.glm")
    expect_true(is.finite(AIC(fit)))
    expect_true(is.finite(BIC(fit)))

    # update() must be able to re-evaluate the recorded call.
    refit <- update(fit, . ~ 1)
    expect_s3_class(refit, "vmideal_glm")
    expect_length(coef(refit), 1L)

    # fitted() refers to the observations of the fit, whose limiting
    # magnitudes the family still holds, and therefore stays correct.
    expect_length(fitted(fit), nrow(obs))
    expect_true(all(is.finite(fitted(fit))))
})

test_that("vmideal moments match a direct summation", {
    # The reference must not go through dvmideal(), which switches to the
    # geometric model once `lm + 10 < psi`. That shortcut differs from
    # .dmideal_int() by a small amount right at the switchover, and the
    # moments are built on the latter.
    reference <- \(lim_magn, psi) {
        m <- seq(.vmideal_upper_lm(lim_magn), -400.0, -1.0)
        d <- .dmideal_int(m, psi) * vmperception(lim_magn - m)
        d <- d / sum(d)
        mu <- sum(m * d)

        c(mu, sum((m - mu)^2 * d))
    }

    model <- expand.grid(
        lim_magn = c(4.5, 5.5, 5.52, 5.8, 6.0, 6.1, 6.3, 6.48, 6.5, 7.5),
        psi = c(-10.0, -5.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 16.0, 30.0)
    )

    for (i in seq_len(nrow(model))) {
        moments <- .vmideal_moments(model$lim_magn[i], model$psi[i], vmperception)
        expected <- reference(model$lim_magn[i], model$psi[i])

        expect_equal(moments$mu, expected[1L], tolerance = 1e-8)
        expect_equal(moments$var, expected[2L], tolerance = 1e-8)
    }

    # Evaluating a whole model frame at once must not change the result.
    set.seed(2L)
    lim_magn <- sample(c(5.5, 5.8, 6.0, 6.3, 6.5), 50L, replace = TRUE)
    psi <- runif(50L, -8.0, 9.0)
    moments <- .vmideal_moments(lim_magn, psi, vmperception)

    expect_equal(
        moments$mu,
        vapply(seq_along(psi), \(i) .vmideal_moments(lim_magn[i], psi[i], vmperception)$mu, numeric(1L))
    )
    expect_equal(
        moments$var,
        vapply(seq_along(psi), \(i) .vmideal_moments(lim_magn[i], psi[i], vmperception)$var, numeric(1L))
    )
})

test_that("the ideal distribution converges to the geometric model", {
    # As psi grows without bound the ideal distribution becomes geometric with
    # r = 10^0.4. That limit supplies the upper end of the attainable mean.
    for (lim_magn in c(5.5, 6.0, 6.5)) {
        ideal <- .vmideal_moments(lim_magn, 30.0, vmperception)
        geom <- .vmgeom_moments(lim_magn, 10^0.4, vmperception)

        expect_equal(ideal$mu, geom$mu, tolerance = 1e-8)
        expect_equal(ideal$var, geom$var, tolerance = 1e-8)
    }

    range <- .vmideal_mu_range(6.0, vmperception, c(-Inf, Inf))
    expect_equal(range$lower, -Inf)
    expect_equal(range$upper, .vmgeom_moments(6.0, 10^0.4, vmperception)$mu)
})

test_that("an unreachable psi falls back to the geometric model", {
    # The ideal density is evaluated relative to psi, so a large psi drives the
    # whole table into the far tail until it underflows. The moments would then
    # be 0/0, and a merely large but finite bound would poison `validmu` and
    # abort the fit deep inside glm.fit().
    for (lim_magn in c(5.5, 6.0, 6.5)) {
        expected <- .vmgeom_moments(lim_magn, 10^0.4, vmperception)

        for (psi in c(Inf, 1e3, 850.0, 800.0, 780.0)) {
            moments <- .vmideal_moments(lim_magn, psi, vmperception)

            expect_equal(moments$mu, expected$mu, tolerance = 1e-8)
            expect_equal(moments$var, expected$var, tolerance = 1e-8)
        }
    }

    # The log density has to follow the same limit.
    magn <- seq(6.0, -10.0)
    expect_equal(
        .vmideal_log_dens(magn, 6.0, Inf, vmperception),
        .vmgeom_log_dens(magn, 6.0, 10^0.4, vmperception)
    )

    # Rows that have reached the limit must not disturb their neighbours.
    lim_magn <- c(6.0, 6.0, 6.0)
    psi <- c(4.0, Inf, 5.0)
    moments <- .vmideal_moments(lim_magn, psi, vmperception)
    expect_equal(moments$mu[2L], .vmgeom_moments(6.0, 10^0.4, vmperception)$mu, tolerance = 1e-8)
    expect_equal(moments$mu[c(1L, 3L)], .vmideal_moments(c(6.0, 6.0), c(4.0, 5.0), vmperception)$mu)

    expect_equal(
        .vmideal_log_dens(c(3.0, 3.0, 3.0), lim_magn, psi, vmperception)[c(1L, 3L)],
        .vmideal_log_dens(c(3.0, 3.0), c(6.0, 6.0), c(4.0, 5.0), vmperception)
    )

    # An infinite upper bound is a legitimate range, and a large finite one
    # must fit exactly as an unbounded one does.
    set.seed(9L)
    lim_magn <- rep(c(5.8, 6.1, 6.4), each = 100L)
    magn <- mapply(\(lim_magn) rvmideal(1L, lim_magn, 4.5), lim_magn)
    obs <- data.frame(magn = magn, lim_magn = lim_magn)

    unbounded <- coef(vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs))
    for (upper in c(500.0, 900.0, Inf)) {
        expect_equal(
            coef(vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs, psi_range = c(-Inf, upper))),
            unbounded
        )
    }
})

test_that("vmideal_glm reuses the location parameter without changing the fit", {
    set.seed(31L)
    x <- rep(seq(-1.2, 1.2, length.out = 8L), each = 60L)
    lim_magn <- rep(c(5.8, 6.2, 6.4), length.out = length(x))
    magn <- mapply(
        \(lim_magn, psi) rvmideal(1L, lim_magn, psi),
        lim_magn, 4.5 + 0.5 * x
    )
    obs <- data.frame(x = x, lim_magn = lim_magn, magn = magn)

    fit <- vmideal_glm(cbind(magn, lim_magn) ~ x, data = obs)

    # `linkinv()` remembers which location parameter produced the current mean,
    # so that `variance()` need not invert it numerically. Discarding that
    # memory after every call forces the inversion and must not change the fit.
    family <- .vmideal_family()
    env <- environment(family$linkinv)$env
    linkinv <- family$linkinv
    family$linkinv <- \(eta) {
        mu <- linkinv(eta)
        env$mu <- NULL
        env$psi <- NULL

        mu
    }
    fit_uncached <- stats::glm(
        cbind(magn, lim_magn) ~ x,
        family = family, data = obs
    )

    expect_equal(coef(fit), coef(fit_uncached), tolerance = 1e-6)
})

test_that("vmideal log density matches dvmideal", {
    # Half-magnitude limits are the delicate cases: `round()` breaks ties to
    # the nearest even number, so 4.5 and 6.5 round down while 5.5 and 7.5
    # round up.
    model <- expand.grid(
        lim_magn = c(4.5, 5.5, 5.8, 6.0, 6.1, 6.3, 6.48, 6.5, 7.5),
        psi = c(-10.0, -5.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 30.0)
    )
    for (i in seq_len(nrow(model))) {
        magn <- seq(.vmideal_upper_lm(model$lim_magn[i]), -25.0)
        expect_equal(
            .vmideal_log_dens(magn, model$lim_magn[i], model$psi[i], vmperception),
            dvmideal(magn, model$lim_magn[i], model$psi[i], log = TRUE)
        )
    }

    # A location parameter that differs for every observation is the case the
    # family relies on, and the one dvmideal() itself handles slowly.
    set.seed(3L)
    n <- 500L
    lim_magn <- rep(c(5.7, 6.0, 6.3, 6.5, 6.1), length.out = n)
    psi <- rnorm(n, 5.0, 0.5)
    magn <- mapply(\(lim_magn, psi) rvmideal(1L, lim_magn, psi), lim_magn, psi)

    expect_equal(
        .vmideal_log_dens(magn, lim_magn, psi, vmperception),
        dvmideal(magn, lim_magn, psi, log = TRUE)
    )

    # Magnitudes outside the support are impossible, not merely unlikely.
    expect_equal(.vmideal_log_dens(-Inf, 6.0, 5.0, vmperception), -Inf)
    expect_equal(.vmideal_log_dens(8.0, 6.0, 5.0, vmperception), -Inf)
})

test_that("the mean is strictly increasing in psi", {
    # `mu.eta` differentiates the mean numerically. The mean comes from
    # pmideal(), a spline whose relative accuracy in psi is about 1e-9, so a
    # step of 1e-6 -- the one the geometric family uses -- drowns in that noise
    # and even returns the wrong sign. This is the regression guard for the
    # step size.
    h <- 0.05
    psi <- seq(-12.0, 9.0, 0.5)

    for (lim_magn in c(5.5, 6.0, 6.5)) {
        derivative <- vapply(psi, \(psi) {
            (.vmideal_moments(lim_magn, psi + h, vmperception)$mu -
                .vmideal_moments(lim_magn, psi - h, vmperception)$mu) / (2.0 * h)
        }, numeric(1L))

        expect_true(all(derivative > 0.0))
    }

    # The inversion of the mean must return the location parameter it started
    # from. It loses accuracy as the mean saturates, hence the widening
    # tolerance at the upper end.
    for (psi in c(-10.0, -5.0, 0.0, 3.0, 5.0, 7.0)) {
        mu <- .vmideal_moments(6.0, psi, vmperception)$mu
        expect_equal(
            .vmideal_psi_from_mu(mu, 6.0, vmperception, c(-Inf, Inf)),
            psi,
            tolerance = 1e-3
        )
    }
})

test_that("vmideal_glm reports an unidentified psi as infinite", {
    # The mean of the distribution rises towards the geometric limit as psi
    # grows. Magnitudes whose mean lies at or above that limit are not an
    # error: the model then says they are geometric, and every larger psi
    # describes them equally well. The estimate is unbounded above, and `Inf`
    # is what says so.
    magn <- c(rep(3L, 40L), rep(4L, 40L), rep(5L, 20L))
    obs <- data.frame(magn = magn, lim_magn = 5.5)

    # This sample mean is far above anything the model can produce.
    expect_gt(mean(magn), .vmgeom_moments(5.5, 10^0.4, vmperception)$mu)

    fit <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs)
    expect_true(fit$converged)
    expect_equal(as.numeric(predict(fit, obs[1L, ])), Inf)

    # The linear predictor itself stays finite; only its interpretation as a
    # location parameter is unbounded.
    expect_true(is.finite(coef(fit)))
    expect_true(is.finite(predict(fit, obs[1L, ], type = "link")))

    # Data drawn from a psi deep inside the flat region are equally
    # compatible with every larger psi, so they too report `Inf`.
    set.seed(5L)
    obs <- data.frame(magn = rvmideal(2000L, 5.5, 15.0), lim_magn = 5.5)

    fit <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs)
    expect_true(fit$converged)
    expect_equal(as.numeric(predict(fit, obs[1L, ])), Inf)

    # A sample mean a hair above the attainable bound is ordinary sampling
    # noise at a large psi, and must not abort the fit either.
    set.seed(202L)
    obs <- data.frame(magn = rvmideal(3000L, 5.5, 15.0), lim_magn = 5.5)
    for (psi_range in list(c(-20.0, 20.0), c(-Inf, 20.0), c(-Inf, Inf))) {
        fit <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs, psi_range = psi_range)

        expect_true(fit$converged)
        expect_true(is.finite(coef(fit)))
    }

    # Ordinary data must remain untouched by any of this.
    set.seed(1L)
    obs <- data.frame(magn = rvmideal(2000L, 6.0, 4.5), lim_magn = 6.0)
    expect_silent(fit <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs))
    expect_equal(as.numeric(predict(fit, obs[1L, ])), 4.5, tolerance = 0.2)
})

test_that("vmideal_glm reports a fit that did not converge", {
    # Near the saturation the mean barely moves with psi, and the iteration can
    # run out of steps. The estimate then looks ordinary -- finite, with a
    # standard error -- so the warning of glm.fit() must not be swallowed along
    # with the step-halving messages that merely describe the way there.
    set.seed(2L)
    obs <- data.frame(magn = rvmideal(1500L, 6.0, 12.0), lim_magn = 6.0)

    fit <- NULL
    warnings <- character(0L)
    withCallingHandlers(
        fit <- vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs),
        warning = \(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    # Either the fit converges, or it says so.
    expect_true(fit$converged || any(grepl("did not converge", warnings, fixed = TRUE)))

    # The mechanism of the step halving stays quiet, whatever the outcome.
    expect_false(any(grepl("step size truncated", warnings, fixed = TRUE)))
    expect_false(any(grepl("boundary value", warnings, fixed = TRUE)))
})

test_that("an unidentified psi does not depend on unrelated observations", {
    # Identifiability is governed by the distance between psi and the limiting
    # magnitude of the observation itself. Adding rows with a different
    # limiting magnitude must not change the verdict for the existing ones.
    set.seed(11L)
    obs <- data.frame(magn = rvmideal(2000L, 5.5, 15.0), lim_magn = 5.5)
    mixed <- rbind(obs, data.frame(magn = c(-2L, -3L, -2L), lim_magn = 7.5))

    alone <- suppressWarnings(vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs))
    together <- suppressWarnings(vmideal_glm(cbind(magn, lim_magn) ~ 1, data = mixed))

    expect_equal(as.numeric(predict(alone, obs[1L, ])), Inf)
    expect_equal(as.numeric(predict(together, obs[1L, ])), Inf)
})

test_that("vmideal_glm validates the range of psi", {
    expect_error(.vmideal_family(psi_range = c(3.0, 2.0)), "increasing")
    expect_error(.vmideal_family(psi_range = 2.0), "two values")
    expect_error(.vmideal_family(psi_range = c(NA_real_, 2.0)), "two values")

    # An infinite bound is a limit, never a value the fit may take, so it is
    # only meaningful on the side that opens the range.
    expect_error(.vmideal_family(psi_range = c(Inf, Inf)), "two values")
    expect_error(.vmideal_family(psi_range = c(-Inf, -Inf)), "two values")
    expect_silent(.vmideal_family(psi_range = c(-Inf, Inf)))
    expect_silent(.vmideal_family(psi_range = c(0.0, Inf)))

    # A missing location parameter is not a limit and must not be answered
    # with one.
    expect_error(.vmideal_is_geometric(6.0, NA_real_, vmperception), "NA")

    # A finite range must reach the family and restrict the fit.
    set.seed(9L)
    lim_magn <- rep(c(5.8, 6.1, 6.4), each = 200L)
    magn <- mapply(\(lim_magn) rvmideal(1L, lim_magn, 4.5), lim_magn)
    obs <- data.frame(magn = magn, lim_magn = lim_magn)

    expect_equal(
        coef(vmideal_glm(cbind(magn, lim_magn) ~ 1, data = obs, psi_range = c(-12.0, 9.0))),
        coef(stats::glm(
            cbind(magn, lim_magn) ~ 1,
            family = .vmideal_family(psi_range = c(-12.0, 9.0)), data = obs
        ))
    )
})

test_that("vmideal_glm passes its arguments to the model family", {
    set.seed(9L)
    lim_magn <- rep(c(5.8, 6.1, 6.4), each = 200L)
    magn <- mapply(\(lim_magn) rvmideal(1L, lim_magn, 4.5), lim_magn)
    obs <- data.frame(magn = magn, lim_magn = lim_magn)

    perception_fun <- \(dm) ifelse(dm > -0.5, 1.0 - exp(-0.005 * (dm + 0.5)^3), 0.0)
    expect_equal(
        coef(vmideal_glm(
            cbind(magn, lim_magn) ~ 1,
            data = obs, perception_fun = perception_fun
        )),
        coef(stats::glm(
            cbind(magn, lim_magn) ~ 1,
            family = .vmideal_family(perception_fun = perception_fun), data = obs
        ))
    )
})

test_that("vmideal_glm rejects an invalid response", {
    obs <- data.frame(
        magn = c(1L, 2L, 3L, 4L),
        lim_magn = c(6.1, 6.1, 5.8, 5.8)
    )

    expect_error(
        vmideal_glm(magn ~ 1, data = obs),
        "two column matrix"
    )
})
