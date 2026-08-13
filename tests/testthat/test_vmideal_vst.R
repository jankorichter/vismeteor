test_that("vmideal_vst", {
    psi <- seq(-3.5, 6.5, 0.5)
    limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)
    m <- seq(-200, 6, 1)

    model <- expand.grid(psi = psi, limmag = limmag)
    model <- do.call(
        rbind.data.frame,
        mapply(function(psi, limmag) {
            p <- dvmideal(m, limmag, psi)
            t <- vmideal_vst_from_magn(m, limmag)
            t.mean <- sum(p * t)
            t.var <- sum(p * (t - t.mean)^2)

            list(
                psi = psi,
                limmag = limmag,
                t.mean = t.mean,
                t.var = t.var
            )
        }, model$psi, model$limmag, SIMPLIFY = FALSE)
    )

    #
    # test vmideal_vst_from_magn
    #
    expect_true(all(abs(model$t.var - 1.0) < 0.0017))

    #
    # test vmideal_vst_to_psi
    #
    model$psi.est <- vmideal_vst_to_psi(model$t.mean, model$limmag)
    expect_true(mean(abs(model$psi - model$psi.est)) < 0.0072)
    expect_true(all(abs(model$psi - model$psi.est) < 0.03))

    # test with low psi values
    model0 <- with(model, {
        subset(model, psi < 5.0)
    })
    expect_true(mean(abs(model0$psi - model0$psi.est)) < 0.006)
    expect_true(all(abs(model0$psi - model0$psi.est) < 0.023))

    # test first derivative
    f <- function(x) {
        vmideal_vst_to_psi(x, 6.0, deriv_degree = 1L)
    }
    y <- vmideal_vst_to_psi(5.5, 6.0) - vmideal_vst_to_psi(4.5, 6.0)
    expect_true(abs(y - stats::integrate(f, 4.5, 5.5)$value) < 1e-10)

    # test second derivative
    f <- function(x) {
        vmideal_vst_to_psi(x, 6.0, deriv_degree = 2L)
    }
    y <- vmideal_vst_to_psi(5.5, 6.0, deriv_degree = 1L) - vmideal_vst_to_psi(4.5, 6.0, deriv_degree = 1L)
    expect_true(abs(y - stats::integrate(f, 4.5, 5.5)$value) < 1e-10)
})

test_that("vmideal_vst_to_psi covers the documented range of psi", {
    # What the transformation resolves is the distance between `psi` and the
    # limiting magnitude, so the range is stated relative to it. These bounds
    # are the ones the documentation gives; the estimate has to be usable
    # throughout and must not degrade silently at the edges.
    m <- seq(-200, 6, 1)
    round_trip <- function(psi, limmag) {
        p <- dvmideal(m, limmag, psi)
        tm <- sum(p * vmideal_vst_from_magn(m, limmag))
        as.vector(vmideal_vst_to_psi(tm, limmag))
    }

    for (limmag in c(5.5, 6.0, 6.5)) {
        for (delta in c(-16, -10, -5, 0, 1, 2, 3)) {
            psi <- limmag + delta
            expect_equal(round_trip(psi, limmag), psi, tolerance = 0.15)
        }
    }

    # Beyond the upper end the mean no longer identifies `psi`, which is
    # reported as the limit rather than as a value the data do not support.
    expect_equal(round_trip(6.0 + 6, 6.0), Inf)
})

test_that("vmideal_vst_to_psi reports an unbounded psi", {
    # A vanishing mean of the transformed magnitudes is what an unbounded `psi`
    # produces, and is reported as such rather than as a missing value.
    # The result carries the 1d-array shape the function has always returned.
    expect_equal(as.vector(vmideal_vst_to_psi(0.0005, 5.5)), Inf)
    expect_equal(as.vector(vmideal_vst_to_psi(0.0, 5.5)), Inf)

    # `psi` there is the geometric limit, which `dvmideal()` evaluates for an
    # infinite location parameter. This is what makes `Inf` the answer rather
    # than an approximation of one.
    m <- seq(5L, 1L)
    expect_identical(
        vismeteor::dvmideal(m, 5.5, as.vector(vmideal_vst_to_psi(0.0005, 5.5))),
        vismeteor::dvmgeom(m, 5.5, 10^0.4)
    )

    # Everything below the threshold is the limit, not a missing value. The
    # threshold that bounds the polynomial and the one that reports an
    # unbounded `psi` are the same, so no band between them is left as `NA`.
    expect_equal(as.vector(vmideal_vst_to_psi(0.001, 5.5)), Inf)
    expect_equal(as.vector(vmideal_vst_to_psi(0.005, 5.5)), Inf)
    expect_equal(as.vector(vmideal_vst_to_psi(0.016, 5.5)), Inf)
    expect_equal(as.vector(vmideal_vst_to_psi(0.0199, 5.5)), Inf)

    # The threshold is exclusive: at and above it the polynomial answers.
    expect_false(is.na(vmideal_vst_to_psi(0.02, 5.5)))
    expect_false(is.na(vmideal_vst_to_psi(0.05, 5.5)))

    # The upper threshold is untouched.
    expect_false(is.na(vmideal_vst_to_psi(8.22, 6.5)))
    expect_true(is.na(vmideal_vst_to_psi(8.23, 6.5)))

    # The derivatives serve the delta method, which needs a finite `psi`.
    for (deriv_degree in c(1L, 2L)) {
        expect_true(is.na(vmideal_vst_to_psi(0.0005, 5.5, deriv_degree = deriv_degree)))
        expect_false(is.na(vmideal_vst_to_psi(0.5, 5.5, deriv_degree = deriv_degree)))
    }

    # Vectorized, and a missing `tm` stays missing rather than turning into a
    # limit.
    psi <- vmideal_vst_to_psi(c(0.5, 0.0005, NA, 9.0), 5.5)
    expect_false(is.na(psi[1]))
    expect_equal(as.vector(psi[2]), Inf)
    expect_true(is.na(psi[3]))
    expect_true(is.na(psi[4]))

    # `psi` is unbounded whatever the limiting magnitude is.
    expect_equal(
        as.vector(vmideal_vst_to_psi(c(0.0005, 0.0005), c(5.5, 6.5))),
        c(Inf, Inf)
    )

    # The threshold is one-sided. A negative mean lies below it as well, and is
    # answered with the same limit rather than with a separate value.
    expect_equal(as.vector(vmideal_vst_to_psi(-0.001, 5.5)), Inf)
    expect_equal(as.vector(vmideal_vst_to_psi(-1.0, 5.5)), Inf)
})
