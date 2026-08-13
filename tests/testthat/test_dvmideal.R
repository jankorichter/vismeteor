test_that("dvmideal", {
    lm <- 6.5
    psi <- 4.0

    # from documentation
    dp_fun <- function(m, psi) {
        r <- 10^0.4
        1.5 * log(r) * sqrt(r^(3 * psi + 2 * m) / ((r^psi + r^m)^5))
    }

    for (current.psi in c(-10.0, psi, 5)) {
        m <- seq(6, current.psi - 25, -1)
        expected_p <- sapply(m, function(m) {
            stats::integrate(function(m) dp_fun(m, current.psi), m - 0.5, m + 0.5)$value *
                vismeteor::vmperception(lm - m)
        })
        expected_p <- expected_p / sum(expected_p)

        # density of meteor magnitudes
        p <- vismeteor::dvmideal(m, lm, current.psi)
        expect_type(p, "double")
        expect_length(p, length(m))
        expect_false(any(abs(log(p) - log(expected_p)) > 0.03), label = paste("psi =", current.psi))

        # log density of meteor magnitudes
        p <- vismeteor::dvmideal(m, lm, current.psi, log = TRUE)
        expect_type(p, "double")
        expect_length(p, length(m))
        expect_false(any(abs(p - log(expected_p)) > 0.03), label = paste("psi =", current.psi))

        # Is the density sum of meteor magnitudes equal to 1.0 ?
        p <- sum(vismeteor::dvmideal(m, lm, current.psi))
        expect_equal(round(p, 6), 1.0, label = paste("psi =", current.psi))
    }

    p <- vismeteor::dvmideal(7, lm, psi)
    expect_type(p, "double")
    expect_length(p, 1)
    expect_equal(p, 0.0)

    # log density of meteor magnitudes
    p <- vismeteor::dvmideal(7, lm, psi, log = TRUE)
    expect_type(p, "double")
    expect_length(p, 1)
    expect_equal(p, -Inf)

    p <- vismeteor::dvmideal(c(6, 6, 6), lm = c(5.4, 5.5, 5.6), psi)
    expect_type(p, "double")
    expect_length(p, 3)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_gt(p[3], 0.0)

    # density of infinite meteor magnitudes
    p <- vismeteor::dvmideal(c(-Inf, Inf), lm, psi)
    expect_type(p, "double")
    expect_length(p, 2)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)

    # log density of infinite meteor magnitudes
    p <- vismeteor::dvmideal(c(-Inf, Inf), lm, psi, log = TRUE)
    expect_type(p, "double")
    expect_length(p, 2)
    expect_equal(p[1], -Inf)
    expect_equal(p[2], -Inf)

    # test order of densities of meteor magnitudes
    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::dvmideal(rep(3, length(lm)), lm, psi)
    expect_equal(p[order(p)], p)

    # test maximum likelihood estimation (MLE) of psi
    lm <- 5.8
    m <- as.integer(seq(-20, 6))
    p <- vismeteor::dvmideal(m, lm, psi)
    ll <- function(psi) { # log likelihood function
        -sum(p * vismeteor::dvmideal(m, lm, psi, log = TRUE))
    }
    est <- optim(1, ll, method = "Brent", lower = 0, upper = 7, hessian = TRUE)
    expect_equal(round(est$par, 6), psi)

    # density of meteor magnitudes equals geometric distribution
    lm <- 6.3
    m <- as.integer(seq(-20, 6))

    p <- vismeteor::dvmideal(m, lm, 16.25, log = TRUE)
    expect_type(p, "double")
    expect_length(p, length(m))
    diff.p <- abs(p - dvmgeom(m, lm, 10^0.4, log = TRUE))
    expect_true(any(diff.p > 0.0))
    expect_false(any(diff.p > 1e-08))

    p <- vismeteor::dvmideal(m, lm, 16.35, log = TRUE)
    expect_type(p, "double")
    expect_length(p, length(m))
    expect_equal(p, dvmgeom(m, lm, 10^0.4, log = TRUE))
})

test_that(".dmideal_int is vectorized over psi", {
    # `psi` selects between the two asymptotic branches and the tabulated
    # middle, so a value paired with the wrong magnitude lands in the wrong
    # branch. dvmideal() splits to a scalar psi and never notices.
    m <- c(0, 0, 0, 0)
    psi <- c(-15.0, 5.0, 15.0, 5.0)

    expect_equal(
        .dmideal_int(m, psi),
        mapply(\(m, psi) .dmideal_int(m, psi), m, psi)
    )

    # Every branch, both scales, and no recycling warning.
    set.seed(17L)
    m <- sample(seq(-40L, 20L), 400L, replace = TRUE)
    psi <- runif(400L, -25.0, 25.0)

    for (log in c(FALSE, TRUE)) {
        expect_silent(d <- .dmideal_int(m, psi, log = log))
        expect_equal(d, mapply(\(m, psi) .dmideal_int(m, psi, log = log), m, psi))
    }

    # The shorter argument is recycled, whichever one it is.
    expect_equal(.dmideal_int(0.0, c(-15.0, 5.0, 15.0)), .dmideal_int(c(0.0, 0.0, 0.0), c(-15.0, 5.0, 15.0)))
    expect_equal(.dmideal_int(c(-2.0, 0.0, 2.0), 5.0), .dmideal_int(c(-2.0, 0.0, 2.0), c(5.0, 5.0, 5.0)))
})

test_that("dvmideal with an infinite psi", {
    m <- seq(6L, -4L)
    lm <- 6.5

    # `Inf` is the geometric limit itself, not an approximation of it, so it
    # has to agree with the geometric model exactly rather than closely.
    expect_identical(
        vismeteor::dvmideal(m, lm, Inf),
        vismeteor::dvmgeom(m, lm, 10^0.4)
    )
    expect_identical(
        vismeteor::dvmideal(m, lm, Inf, log = TRUE),
        vismeteor::dvmgeom(m, lm, 10^0.4, log = TRUE)
    )

    # A large finite `psi` reaches the same limit, so the two must not differ.
    expect_equal(vismeteor::dvmideal(m, lm, Inf), vismeteor::dvmideal(m, lm, 30.0))

    # The limit is taken per observation. An infinite `psi` in one row must
    # leave the remaining rows where they are.
    d <- vismeteor::dvmideal(c(3L, 3L, 3L), lm, c(4.0, Inf, 4.0))
    expect_equal(d[1], d[3])
    expect_false(isTRUE(all.equal(d[1], d[2])))

    # The geometric branch has to keep the perception probability function it
    # was given instead of falling back to the default.
    perception_fun <- function(x) vismeteor::vmperception(x)^0.9
    expect_equal(
        vismeteor::dvmideal(m, lm, Inf, perception_fun = perception_fun),
        vismeteor::dvmgeom(m, lm, 10^0.4, perception_fun = perception_fun)
    )

    # There is no limiting distribution at the opposite end.
    expect_error(vismeteor::dvmideal(m, lm, -Inf), "psi")
})

test_that("dvmideal falls back to the geometric model for a large finite psi", {
    m <- seq(6L, -4L)
    lm <- 6.5

    # The switch is `lm + 10 < psi`. Just below it the distribution is still
    # tabulated, just above it the geometric model answers, and both sides have
    # to agree -- otherwise the fallback would introduce a step.
    expect_equal(
        vismeteor::dvmideal(m, lm, lm + 10.0 - 0.1),
        vismeteor::dvmideal(m, lm, lm + 10.0 + 0.1)
    )
    expect_identical(
        vismeteor::dvmideal(m, lm, lm + 10.0 + 0.1),
        vismeteor::dvmgeom(m, lm, 10^0.4)
    )

    for (log in c(FALSE, TRUE)) {
        expect_equal(
            vismeteor::dvmideal(m, lm, 20.0, log = log),
            vismeteor::dvmgeom(m, lm, 10^0.4, log = log)
        )
    }
})
