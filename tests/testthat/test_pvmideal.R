test_that("pvmideal", {
    psi <- 4.0

    # probability of m >= 3
    expected_p <- round(sum(vismeteor::dvmideal(as.integer(seq(3, 6)), 6.3, psi)), 6)
    expect_equal(expected_p, 0.339967)

    # probability
    p <- vismeteor::pvmideal(c(3, 7, -100, -Inf), 6.3, psi, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(round(p[1], 6), 1 - expected_p)
    expect_equal(p[2], 1.0)
    expect_equal(p[3], 0.0)
    expect_equal(p[4], 0.0)

    p <- vismeteor::pvmideal(c(3, 7, -100, -Inf), 6.3, psi, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(round(p[1], 6), expected_p)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 1.0)
    expect_equal(p[4], 1.0)

    # log probability
    p <- vismeteor::pvmideal(c(3, 7, -Inf), 6.3, psi, lower.tail = TRUE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(round(p[1], 4), round(log(1 - expected_p), 4))
    expect_equal(p[2], 0.0)
    expect_equal(p[3], -Inf)

    p <- vismeteor::pvmideal(c(3, 7, -Inf), 6.3, psi, lower.tail = FALSE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(round(p[1], 4), round(log(expected_p), 4))
    expect_equal(p[2], -Inf)
    expect_equal(p[3], 0.0)

    # test order of probabilities
    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::pvmideal(rep(3L, length(lm)), lm, psi, lower.tail = TRUE)
    expect_equal(p[order(p, decreasing = TRUE)], p)

    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::pvmideal(rep(3L, length(lm)), lm, psi, lower.tail = FALSE)
    expect_equal(p[order(p, decreasing = FALSE)], p)

    # test lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- vismeteor::dvmideal(rep(6L, length(lm)), lm, psi)
    p2 <- vismeteor::pvmideal(rep(6L, length(lm)), lm, psi, lower.tail = FALSE)
    expect_equal(p1, p2)

    # test sum lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- vismeteor::dvmideal(rep(6L, length(lm)), lm, psi) +
        vismeteor::dvmideal(rep(5L, length(lm)), lm, psi)
    p2 <- vismeteor::pvmideal(rep(5L, length(lm)), lm, psi, lower.tail = FALSE)
    expect_equal(p1, p2)

    # density of meteor magnitudes equals geometric distribution
    lm <- 5.8
    m <- as.integer(seq(-20, 6))
    p <- vismeteor::pvmideal(m, lm, 30, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(pvmgeom(m, 10^(1/2.5), lm = lm, lower.tail = FALSE), p)

    p <- vismeteor::pvmideal(m, lm, 30, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(pvmgeom(m, 10^(1/2.5), lm = lm, lower.tail = TRUE), p)
})
