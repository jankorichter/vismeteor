test_that("pvmgeom", {
    r <- 1.8
    perception.const <- function(m, log = FALSE) {
        rep(ifelse(log, 0.0, 1.0), length(m))
    }

    # probability of (limmag - m) = 2
    f <- function(m) {
        dvmgeom(m, r)
    }
    expected_p <- round(sum(dvmgeom(as.integer(seq(0, 2)), r)), 6)
    expect_equal(expected_p, 0.184265)

    # probability
    p <- pvmgeom(c(2, -0.6, -0.5, 100, Inf), r, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 6), expected_p)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 0.0)
    expect_equal(p[4], 1.0)
    expect_equal(p[5], 1.0)

    p <- pvmgeom(c(2, -0.6, -0.5, 100, Inf), r, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 6), 1 - expected_p)
    expect_equal(p[2], 1.0)
    expect_equal(p[3], 1.0)
    expect_equal(p[4], 0.0)
    expect_equal(p[5], 0.0)

    # log probability
    p <- pvmgeom(c(2, -0.6, -0.5, 100, Inf), r, log = TRUE, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 4), round(log(expected_p), 4))
    expect_equal(p[2], -Inf)
    expect_equal(p[3], -Inf)
    expect_equal(p[4], 0.0)
    expect_equal(p[5], 0.0)

    p <- pvmgeom(c(2, -0.6, -0.5, 100, Inf), r, log = TRUE, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 4), round(log(1 - expected_p), 4))
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 0.0)
    expect_lt(p[4], -40)
    expect_equal(p[5], -Inf)

    # probability with limiting magnitude
    p <- pvmgeom(c(4, 6.6, 6.5, -100, -Inf), r, lm = 6.0, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 6), expected_p)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 0.0)
    expect_equal(p[4], 1.0)
    expect_equal(p[5], 1.0)

    p <- pvmgeom(c(4, 6.6, 6.5, -100, -Inf), r, lm = 6.0, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 6), 1.0 - expected_p)
    expect_equal(p[2], 1.0)
    expect_equal(p[3], 1.0)
    expect_equal(p[4], 0.0)
    expect_equal(p[5], 0.0)

    p <- pvmgeom(c(6, 6, 6), r, lm = c(5.4, 5.5, 5.6), lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_gt(p[3], 0.0)

    p <- pvmgeom(c(6, 6, 6), r, lm = c(5.4, 5.5, 5.6), lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], 1.0)
    expect_equal(p[2], 1.0)
    expect_lt(p[3], 1.0)

    # log probability with limiting magnitude
    p <- pvmgeom(c(4, 6.6, 6.5, -100, -Inf), r, lm = 6.0, log = TRUE, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 4), round(log(expected_p), 4))
    expect_equal(p[2], -Inf)
    expect_equal(p[3], -Inf)
    expect_equal(p[4], 0.0)
    expect_equal(p[5], 0.0)

    p <- pvmgeom(c(4, 6.6, 6.5, -100, -Inf), r, lm = 6.0, log = TRUE, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 4), round(log(1.0 - expected_p), 4))
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 0.0)
    expect_lt(p[4], -40)
    expect_equal(p[5], -Inf)

    # test order of probabilities
    lm <- seq(2.5, 6.5, 0.1)
    p <- pvmgeom(rep(3L, length(lm)), r, lm = lm, lower.tail = TRUE)
    expect_equal(p[order(p, decreasing = FALSE)], p)

    lm <- seq(2.5, 6.5, 0.1)
    p <- pvmgeom(rep(3L, length(lm)), r, lm = lm, lower.tail = FALSE)
    expect_equal(p[order(p, decreasing = TRUE)], p)

    # test lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- dvmgeom(rep(6L, length(lm)), r, lm = lm)
    p2 <- pvmgeom(rep(6L, length(lm)), r, lm = lm)
    expect_equal(p1, p2)

    # test sum lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- dvmgeom(rep(6L, length(lm)), r, lm = lm) +
        dvmgeom(rep(5L, length(lm)), r, lm = lm)
    p2 <- pvmgeom(rep(5L, length(lm)), r, lm = lm)
    expect_equal(p1, p2)

    # probability of meteor magnitudes equals geometric distribution
    m <- as.integer(seq(0, 30, 1))
    p <- pvmgeom(m, r, perception.fun = perception.const)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(stats::pgeom(m, 1 - 1/r), p)
    m <- as.integer(seq(0, 30, 1))

    p <- pvmgeom(m, r, lower.tail = FALSE, perception.fun = perception.const)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(stats::pgeom(m, 1 - 1/r, lower.tail = FALSE), p)
})
