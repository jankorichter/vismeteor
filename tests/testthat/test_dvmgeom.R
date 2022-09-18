test_that("dvmgeom", {
    r <- 1.8
    expected_p <- 0.113972

    # density of meteor magnitudes
    p <- dvmgeom(c(-0.5, 2.0), r)
    expect_type(p, 'double')
    expect_length(p, 2)
    expect_equal(p[1], 0.0)
    expect_equal(round(p[2], 6), expected_p)

    # log density of meteor magnitudes
    p <- dvmgeom(2.0, r, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 1)
    expect_equal(round(p, 4), round(log(expected_p), 4))

    # density of meteor magnitudes with limiting magnitude
    lm <- seq(5.5, 6.5)
    p <- dvmgeom(lm - 2.0, r, lm = lm)
    expect_type(p, 'double')
    expect_length(p, length(lm))
    expect_equal(round(p, 6), rep(expected_p, length(lm)))

    p <- dvmgeom(c(6, 6, 6), r, lm = c(5.4, 5.5, 5.6))
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_gt(p[3], 0.0)

    # density of infinite meteor magnitudes
    p <- dvmgeom(c(-Inf, Inf), r)
    expect_type(p, 'double')
    expect_length(p, 2)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)

    # log density of infinite meteor magnitudes
    p <- dvmgeom(c(-Inf, Inf), r, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 2)
    expect_equal(p[1], -Inf)
    expect_equal(p[2], -Inf)

    # test order of densities of meteor magnitudes
    lm <- seq(2.5, 6.5, 0.1)
    p <- dvmgeom(rep(3L, length(lm)), r, lm = lm)
    expect_equal(p[order(p)], p)

    # Is the density sum of meteor magnitudes equal to 1.0 ?
    for (lm in seq(5.4, 6.4, 0.1)) {
        p <- sum(dvmgeom(as.integer(seq(6, -100, -1)), r, lm = lm))
        expect_equal(round(p, 6), 1.0, label = paste('lm =', lm))
    }

    # test mean of meteor magnitudes
    q <- log(r)
    f <- function(m) {
        exp(-q * m) * vmperception(m)
    }
    norm <- stats::integrate(f, -0.5, Inf)$value
    f <- function(m) {
        m * exp(-q * m) * vmperception(m)/norm
    }
    expected_m.mean <- stats::integrate(f, -0.5, Inf)$value
    expect_equal(round(expected_m.mean, 2), 4.41)

    f <- function(m) {
        m * dvmgeom(m, r)
    }
    m.mean <- sum(f(as.integer(seq(0, 100, 1))))
    expect_equal(round(m.mean, 2), round(expected_m.mean, 2))

    f <- function(m) {
        m * dvmgeom(m, r, lm = 6.2)
    }
    m.mean <- sum(f(as.integer(seq(6, -100, -1))))
    expect_equal(round(6.2 - m.mean, 2), round(expected_m.mean, 2))

    # test maximum likelihood estimation (MLE) of q of meteor magnitudes
    lm <- 5.8
    m <- as.integer(seq(6, -60, -1))
    p <- dvmgeom(m, r, lm = lm)
    ll <- function(q) { # log likelihood function
        -sum(p * dvmgeom(m, exp(q), lm = lm, log = TRUE))
    }
    est <- optim(1, ll, method='Brent', lower=0.01, upper=4, hessian=TRUE)
    expect_equal(round(est$par, 6), round(base::log(r), 6))

    # test q value of meteor magnitudes
    lm <- 5.8
    m <- as.integer(seq(8, -60, -1))
    p <- dvmgeom(m, r, lm = lm)
    expect_equal(round(exp(sum(p * vismeteor::vmperception(lm - m, log = TRUE, deriv = TRUE))), 3), 1.799)
})
