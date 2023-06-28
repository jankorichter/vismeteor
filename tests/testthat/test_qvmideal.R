test_that("qvmgeom", {
    psi <- 4.0

    lm <- 6.0
    p <- vismeteor::pvmideal(6, lm, psi, lower.tail = FALSE)
    q <- suppressWarnings(
        vismeteor::qvmideal(c(-1.0, 0.0, p - 1e-06, p, p + 1e-06, 1.0), lm, psi, lower.tail = FALSE)
    )
    expect_equal(q, c(NA, 6, 6, 6, 5, -Inf))

    p <- vismeteor::pvmideal(6, lm, psi, lower.tail = TRUE)
    q <- suppressWarnings(
        vismeteor::qvmideal(c(-1.0, 1.0, p + 1e-06, p, p - 1e-06, 0.0), lm, psi, lower.tail = TRUE)
    )
    expect_equal(q, c(NA, 6, 6, 6, 5, -Inf))

    lms <- seq(5.6, 6.5, 0.1)
    for (lm in lms) {
        q <- suppressWarnings(
            vismeteor::qvmideal(c(-1.0, 0.0, 1.0), lm, psi, lower.tail = TRUE)
        )
        expect_equal(q, c(NA, -Inf, 6), label = paste('lm =', lm))

        q <- suppressWarnings(
            vismeteor::qvmideal(c(-1.0, 0.0, 1.0), lm, psi, lower.tail = FALSE)
        )
        expect_equal(q, c(NA, 6, -Inf), label = paste('lm =', lm))
    }

    for (lm in lms) {
        p <- vismeteor::pvmideal(6, lm, psi, lower.tail = TRUE)
        q <- vismeteor::qvmideal(c(p + 1e-06, p -1e-06), lm, psi, lower.tail = TRUE)
        expect_equal(q, c(6, 5), label = paste('lm =', lm))

        p <- vismeteor::pvmideal(6, lm, psi, lower.tail = FALSE)
        q <- vismeteor::qvmideal(c(p - 1e-06, p, p + 1e-06), lm, psi, lower.tail = FALSE)
        expect_equal(q, c(6, 6, 5), label = paste('lm =', lm))
    }

    for (lm in lms) {
        p <- vismeteor::pvmideal(5, lm, psi, lower.tail = TRUE)
        q <- vismeteor::qvmideal(c(p + 1e-06, p - 1e-06), lm, psi, lower.tail = TRUE)
        expect_equal(q, c(5, 4), label = paste('lm =', lm))

        p <- vismeteor::pvmideal(5, lm, psi, lower.tail = FALSE)
        q <- vismeteor::qvmideal(c(p - 1e-06, p, p + 1e-06), lm, psi, lower.tail = FALSE)
        expect_equal(q, c(5, 5, 4), label = paste('lm =', lm))
    }

    # quantile of meteor magnitudes equals geometric distribution
    lm <- 6.0
    m <- as.integer(seq(-10, 6, 1))
    p <- round(vismeteor::pvmideal(m, lm, 30, lower.tail = TRUE), 6)
    m <- vismeteor::qvmideal(p, lm, 30, lower.tail = TRUE)
    expect_type(m, 'double')
    expect_length(m, length(m))
    expect_equal(qvmgeom(p, 10^(1/2.5), lm = lm, lower.tail = FALSE), m)

    lm <- 6.0
    m <- as.integer(seq(-10, 6, 1))
    p <- round(vismeteor::pvmideal(m, lm, 30, lower.tail = FALSE), 6)
    m <- vismeteor::qvmideal(p, lm, 30, lower.tail = FALSE)
    expect_type(m, 'double')
    expect_length(m, length(m))
    expect_equal(vismeteor::qvmgeom(p, 10^(1/2.5), lm = lm, lower.tail = TRUE), m)
})
