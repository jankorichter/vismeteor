test_that("qvmgeom", {
    r <- 1.8

    ####################################
    # Tests without limiting magnitude #
    ####################################

    p <- pvmgeom(0, r, lower.tail = TRUE)
    q <- suppressWarnings(
        qvmgeom(c(-1.0, 0.0, p - 1e-06, p, p + 1e-06, 1.0), r, lower.tail = TRUE)
    )
    expect_equal(q, c(NA, 0, 0, 0, 1, Inf))

    p <- pvmgeom(0, r, lower.tail = FALSE)
    q <- suppressWarnings(
        qvmgeom(c(-1.0, 1.0, p + 1e-06, p, p - 1e-06, 0.0), r, lower.tail = FALSE)
    )
    expect_equal(q, c(NA, 0, 0, 0, 1, Inf))

    p <- pvmgeom(1, r, lower.tail = TRUE)
    q <- suppressWarnings(
        qvmgeom(c(-1.0, 0.0, p - 1e-06, p, p + 1e-06, 1.0), r, lower.tail = TRUE)
    )
    expect_equal(q, c(NA, 0, 1, 1, 2, Inf))

    p <- pvmgeom(1, r, lower.tail = FALSE)
    q <- suppressWarnings(
        qvmgeom(c(-1.0, 1.0, p + 1e-06, p, p - 1e-06, 0.0), r, lower.tail = FALSE)
    )
    expect_equal(q, c(NA, 0, 1, 1, 2, Inf))

    m <- as.integer(seq(0, 19, 1))
    p <- pvmgeom(m, r, lower.tail = TRUE)
    expect_equal(qvmgeom(p, r, lower.tail = TRUE), m)

    m <- as.integer(seq(0, 19, 1))
    p <- pvmgeom(m, r, lower.tail = FALSE)
    expect_equal(qvmgeom(p, r, lower.tail = FALSE), m)

    #################################
    # Tests with limiting magnitude #
    #################################

    p <- pvmgeom(6, r, lm = 5.5, lower.tail = TRUE)
    q <- suppressWarnings(
        qvmgeom(c(-1.0, 0.0, p - 1e-06, p, p + 1e-06, 1.0), r, lm = 5.5, lower.tail = TRUE)
    )
    expect_equal(q, c(NA, 5, NA, 5, 5, -Inf))

    p <- pvmgeom(6, r, lm = 5.5, lower.tail = FALSE)
    q <- suppressWarnings(
        qvmgeom(c(-1.0, 1.0, p + 1e-06, p, p - 1e-06, 0.0), r, lm = 5.5, lower.tail = FALSE)
    )
    expect_equal(q, c(NA, 5, NA, 5, 5, -Inf))

    lms <- seq(5.6, 6.5, 0.1)

    for (lm in lms) {
        q <- suppressWarnings(
            qvmgeom(c(-1.0, 0.0, 1.0), r, lm = lm, lower.tail = TRUE)
        )
        expect_equal(q, c(NA, 6, -Inf), label = paste('lm =', lm))

        q <- suppressWarnings(
            qvmgeom(c(-1.0, 0.0, 1.0), r, lm = lm, lower.tail = FALSE)
        )
        expect_equal(q, c(NA, -Inf, 6), label = paste('lm =', lm))
    }

    for (lm in lms) {
        p <- pvmgeom(6, r, lm = lm, lower.tail = TRUE)
        q <- qvmgeom(c(p - 1e-06, p, p + 1e-06), r, lm = lm, lower.tail = TRUE)
        expect_equal(q, c(6, 6, 5), label = paste('lm =', lm))

        p <- pvmgeom(6, r, lm = lm, lower.tail = FALSE)
        q <- qvmgeom(c(p + 1e-06, p -1e-06), r, lm = lm, lower.tail = FALSE)
        expect_equal(q, c(6, 5), label = paste('lm =', lm))
    }

    for (lm in lms) {
        p <- pvmgeom(5, r, lm = lm, lower.tail = TRUE)
        q <- qvmgeom(c(p - 1e-06, p, p + 1e-06), r, lm = lm, lower.tail = TRUE)
        expect_equal(q, c(5, 5, 4), label = paste('lm =', lm))

        p <- pvmgeom(5, r, lm = lm, lower.tail = FALSE)
        q <- qvmgeom(c(p + 1e-06, p - 1e-06), r, lm = lm, lower.tail = FALSE)
        expect_equal(q, c(5, 4), label = paste('lm =', lm))
    }
})
