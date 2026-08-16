test_that("vmgeom_vst", {
    # the range both functions are calibrated on
    r <- seq(1.4, 4.0, 0.1)
    limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)
    m <- seq(-200, 6, 1)
    model <- expand.grid(r = r, limmag = limmag)
    model <- do.call(
        rbind.data.frame,
        mapply(function(r, limmag) {
            p <- dvmgeom(m, limmag, r)
            t <- vmgeom_vst_from_magn(m, limmag)
            t.mean <- sum(p * t)
            t.var <- sum(p * (t - t.mean)^2)

            # the rate-based statistic the transformation is built on
            den <- vmperception(limmag - m)
            g <- ifelse(den > 0.0, vmperception(limmag - m - 1.0) / den, 0.0)

            list(
                r = r,
                q = log(r),
                limmag = limmag,
                g.mean = sum(p * g),
                t.mean = t.mean,
                t.var = t.var
            )
        }, model$r, model$limmag, SIMPLIFY = FALSE)
    )

    #
    # test vmgeom_vst_from_magn
    #

    # the transformation rests on E[g] = 1/r, which holds exactly and for
    # every limiting magnitude; the calibration of c and d relies on it
    expect_equal(model$g.mean, 1.0 / model$r)

    # outside the support the transform is zero; it must stay finite, as
    # otherwise sum(p * t) above would be poisoned by 0 * NaN
    expect_equal(vmgeom_vst_from_magn(6, 5.5), 0.0)
    expect_false(any(is.nan(vmgeom_vst_from_magn(c(6, 7, 8), 5.5))))

    # a missing magnitude propagates as NA instead of reaching vmperception(),
    # which does not accept one
    expect_true(is.na(vmgeom_vst_from_magn(NA_real_, 5.5)))
    expect_equal(
        vmgeom_vst_from_magn(c(2, NA_real_), 6.0),
        c(vmgeom_vst_from_magn(2, 6.0), NA_real_)
    )

    expect_true(all(abs(model$t.var - 1.0) < 0.17))
    expect_true(all(abs(subset(model, r <= 3.5)$t.var - 1.0) < 0.16))

    # the perception probability increases monotonically, hence the ratio it is
    # built from lies between 0 and 1, which bounds the transform by the value
    # it takes at that ratio
    tm_max <- vmgeom_vst_from_magn(-100, 0)
    dm <- seq(-1.0, 200.0, 0.01)
    expect_true(all(vmgeom_vst_from_magn(0.0, dm) <= tm_max))
    expect_equal(max(vmgeom_vst_from_magn(0.0, dm)), tm_max)

    #
    # test vmgeom_vst_to_r
    #

    # `vmgeom_vst_to_r` is applied to the mean of the transformed magnitudes,
    # so that is what has to be tested. Testing it against
    # `a * mean(g)^b` instead would only restate how c and d are defined and
    # would hold for any b, hence catch nothing.
    # q = log(r) is the scale the deviation is even on; in r it grows with r
    # itself, which is why the bound below is the looser one of the two
    model$q_est <- vmgeom_vst_to_r(model$t.mean, log = TRUE)
    expect_true(all(abs(model$q - model$q_est) < 0.015))

    model$r_est <- vmgeom_vst_to_r(model$t.mean)
    expect_true(all(abs(model$r - model$r_est) < 0.05))
    expect_true(all(abs(subset(model, r <= 3.5)$r -
        subset(model, r <= 3.5)$r_est) < 0.03))

    # The estimator must be consistent: the error may not survive an
    # arbitrarily large sample. The delta-method correction of the vignette
    # scales with var/n and vanishes, so whatever is left at n -> Inf is a
    # property of the back-transformation itself.
    with(new.env(), {
        r <- 2.5
        limmag <- 6.0
        p <- dvmgeom(m, limmag, r)
        t <- vmgeom_vst_from_magn(m, limmag)
        t.mean <- sum(p * t)
        t.var <- sum(p * (t - t.mean)^2)

        for (n in c(1e3, 1e5, 1e7)) {
            r_hat <- vmgeom_vst_to_r(t.mean) -
                0.5 * vmgeom_vst_to_r(t.mean, deriv_degree = 2L) * t.var / n
            expect_true(abs(r_hat - r) < 0.07)
        }
    })

    # the upper bound of tm corresponds to the degenerate case r = 1, which
    # the back-transformation is constrained to reproduce, so the geometric
    # model's requirement r >= 1 holds over the whole domain
    expect_true(vmgeom_vst_to_r(tm_max) >= 1.0)
    expect_equal(vmgeom_vst_to_r(tm_max), 1.0, tolerance = 1e-05)

    # The back-transformation carries a term in `tm` beside the one in
    # log(tm). It must not cost monotonicity anywhere in (0, tm_max], since
    # `vmgeom_vst_to_r` promises an ordered answer for every value the model
    # can produce -- a curved form fitted without care turns around instead.
    with(new.env(), {
        tm_grid <- seq(tm_max / 1e04, tm_max, length.out = 20000L)
        r_grid <- vmgeom_vst_to_r(tm_grid)
        expect_true(all(diff(r_grid) < 0.0))
        expect_true(all(r_grid >= 1.0))
        # and r -> Inf at the faint end, not towards some finite value
        expect_true(vmgeom_vst_to_r(1e-08) > 1e06)
    })

    # unlike the former polynomial there is no calibration window:
    # the inverse stays finite and monotonic between its bounds. The last
    # probe is the upper bound itself, so that the range is covered whatever
    # the calibration makes of it.
    tm_probe <- c(0.5, 1.72, 3.53, tm_max)
    expect_true(all(is.finite(vmgeom_vst_to_r(tm_probe))))
    expect_true(all(diff(vmgeom_vst_to_r(tm_probe)) < 0.0))

    # tm = 0 means g = 0 and hence r = Inf; values outside [0, tm_max]
    # cannot occur under the model
    expect_identical(vmgeom_vst_to_r(0.0), Inf)
    expect_true(is.na(vmgeom_vst_to_r(-1.0)))
    expect_true(is.na(vmgeom_vst_to_r(tm_max + 0.1)))
    # both branches must reject impossible input alike, without a warning
    expect_silent(vmgeom_vst_to_r(-1.0, log = TRUE))

    # test first derivative
    f <- function(x) {
        vmgeom_vst_to_r(x, deriv_degree = 1L)
    }
    y <- vmgeom_vst_to_r(3.0) - vmgeom_vst_to_r(2.0)
    expect_true(abs(y - stats::integrate(f, 2.0, 3.0)$value) < 1e-10)

    # test second derivative
    f <- function(x) {
        vmgeom_vst_to_r(x, deriv_degree = 2L)
    }
    y <- vmgeom_vst_to_r(3.0, deriv_degree = 1L) - vmgeom_vst_to_r(2.0, deriv_degree = 1L)
    expect_true(abs(y - stats::integrate(f, 2.0, 3.0)$value) < 1e-10)

    # log ... both branches must agree with each other
    expect_equal(
        vmgeom_vst_to_r(model$t.mean, log = TRUE),
        base::log(vmgeom_vst_to_r(model$t.mean))
    )

    # test log first derivative
    f <- function(x) {
        vmgeom_vst_to_r(x, log = TRUE, deriv_degree = 1L)
    }
    y <- vmgeom_vst_to_r(3.0, log = TRUE) - vmgeom_vst_to_r(2.0, log = TRUE)
    expect_true(abs(y - stats::integrate(f, 2.0, 3.0)$value) < 1e-10)

    # test log second derivative
    f <- function(x) {
        vmgeom_vst_to_r(x, log = TRUE, deriv_degree = 2L)
    }
    y <- vmgeom_vst_to_r(3.0, log = TRUE, deriv_degree = 1L) -
        vmgeom_vst_to_r(2.0, log = TRUE, deriv_degree = 1L)
    expect_true(abs(y - stats::integrate(f, 2.0, 3.0)$value) < 1e-10)
})
