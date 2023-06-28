test_that("pvmgeom", {
    r <- 1.8
    perception.const <- function(m, log = FALSE) {
        rep(ifelse(log, 0.0, 1.0), length(m))
    }

    # from documentation
    dp.fun <- function(m, log = TRUE) {
        dgeom(m, 1 - 1/r) * vismeteor::vmperception(m)
    }
    m <- seq(0, 40)
    expected_p <- dp.fun(m)
    expected_p <- expected_p/sum(expected_p)
    expected_p <- cumsum(expected_p)

    # probability
    p <- vismeteor::pvmgeom(m, r, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, expected_p)

    p <- vismeteor::pvmgeom(c(-0.6, -0.5, 100, Inf), r, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 1.0)
    expect_equal(p[4], 1.0)

    p <- vismeteor::pvmgeom(c(2, -0.6, -0.5, 100, Inf), r, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(p[1], 1 - expected_p[3])
    expect_equal(p[2], 1.0)
    expect_equal(p[3], 1.0)
    expect_equal(p[4], 0.0)
    expect_equal(p[5], 0.0)

    # log probability
    p <- vismeteor::pvmgeom(m, r, log = TRUE, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(round(p, 4), round(log(expected_p), 4))

    p <- vismeteor::pvmgeom(c(-0.6, -0.5, 100, Inf), r, log = TRUE, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(p[1], -Inf)
    expect_equal(p[2], -Inf)
    expect_equal(p[3], 0.0)
    expect_equal(p[4], 0.0)

    p <- vismeteor::pvmgeom(c(2, -0.6, -0.5, 100, Inf), r, log = TRUE, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 4), round(log(1 - expected_p[3]), 4))
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 0.0)
    expect_lt(p[4], -40)
    expect_equal(p[5], -Inf)

    # probability with limiting magnitude
    p <- vismeteor::pvmgeom(6.0 - m, r, lm = 6.0, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, expected_p)

    p <- vismeteor::pvmgeom(c(6.6, 6.5, -100, -Inf), r, lm = 6.0, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 1.0)
    expect_equal(p[4], 1.0)

    p <- vismeteor::pvmgeom(c(4, 6.6, 6.5, -100, -Inf), r, lm = 6.0, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(p[1], 1.0 - expected_p[3])
    expect_equal(p[2], 1.0)
    expect_equal(p[3], 1.0)
    expect_equal(p[4], 0.0)
    expect_equal(p[5], 0.0)

    p <- vismeteor::pvmgeom(c(6, 6, 6), r, lm = c(5.4, 5.5, 5.6), lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_gt(p[3], 0.0)

    p <- vismeteor::pvmgeom(c(6, 6, 6), r, lm = c(5.4, 5.5, 5.6), lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], 1.0)
    expect_equal(p[2], 1.0)
    expect_lt(p[3], 1.0)

    # log probability with limiting magnitude
    p <- vismeteor::pvmgeom(c(4, 6.6, 6.5, -100, -Inf), r, lm = 6.0, log = TRUE, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 4), round(log(expected_p[3]), 4))
    expect_equal(p[2], -Inf)
    expect_equal(p[3], -Inf)
    expect_equal(p[4], 0.0)
    expect_equal(p[5], 0.0)

    p <- vismeteor::pvmgeom(c(4, 6.6, 6.5, -100, -Inf), r, lm = 6.0, log = TRUE, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 4), round(log(1.0 - expected_p[3]), 4))
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 0.0)
    expect_lt(p[4], -40)
    expect_equal(p[5], -Inf)

    # test probabilities of meteor magnitudes with different parameters
    p1 <- vismeteor::pvmgeom(3, 1.8, lm = 6.5)
    p2 <- vismeteor::pvmgeom(4, 2.0, lm = 6.3)
    p <- vismeteor::pvmgeom(c(3, 4), c(1.8, 2.0), lm = c(6.5, 6.3))
    expect_equal(p, c(p1, p2))

    # test order of probabilities
    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::pvmgeom(rep(3L, length(lm)), r, lm = lm, lower.tail = TRUE)
    expect_equal(p[order(p, decreasing = FALSE)], p)

    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::pvmgeom(rep(3L, length(lm)), r, lm = lm, lower.tail = FALSE)
    expect_equal(p[order(p, decreasing = TRUE)], p)

    # test lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- vismeteor::dvmgeom(rep(6L, length(lm)), r, lm = lm)
    p2 <- vismeteor::pvmgeom(rep(6L, length(lm)), r, lm = lm)
    expect_equal(p1, p2)

    # test sum lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- vismeteor::dvmgeom(rep(6L, length(lm)), r, lm = lm) +
        vismeteor::dvmgeom(rep(5L, length(lm)), r, lm = lm)
    p2 <- vismeteor::pvmgeom(rep(5L, length(lm)), r, lm = lm)
    expect_equal(p1, p2)

    # probability of meteor magnitudes equals geometric distribution
    m <- as.integer(seq(0, 30, 1))
    p <- vismeteor::pvmgeom(m, r, perception.fun = perception.const)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(stats::pgeom(m, 1 - 1/r), p)
    m <- as.integer(seq(0, 30, 1))

    p <- vismeteor::pvmgeom(m, r, lower.tail = FALSE, perception.fun = perception.const)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(stats::pgeom(m, 1 - 1/r, lower.tail = FALSE), p)
})
