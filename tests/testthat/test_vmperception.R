test_that("vmperception", {
    # test approximation of `KOSREN90` perceptions published in
    # J. Rendtel, Handbook for meteor observers 2022 edition, p. 149, (c) International Meteor Organization
    data <- data.frame(
        m = seq(-0.6, 8.0, 0.2),
        p0 = c(
            0.00005, # virtual value
            0.00046, 0.0011, 0.0023, 0.0046, 0.0081, 0.0122, 0.0182,
            0.0257, 0.03444, 0.04365, 0.05495, 0.06918, 0.08511, 0.10351,
            0.13, 0.16, 0.2, 0.24, 0.29, 0.34674, 0.4, 0.46, 0.52, 0.57544,
            0.63096, 0.67764, 0.71, 0.74, 0.77, 0.79, 0.81, 0.83, 0.85,
            0.87, 0.89, 0.91, 0.93, 0.94, 0.96, 0.98, 1.0, 1.0, 1.0
        )
    )

    # regression
    # lm.data <- subset(data, data$p0 < 1.0)
    # lm.data$y <- log(lm.data$p0/(1 - lm.data$p0))
    # res <- stats::lm(y ~ poly(m, degree = 5, raw = TRUE), data = lm.data)
    # print(coef(res))

    # ignore lower magnitudes
    data <- subset(data, data$m >= 0)

    # check perceptions
    data$p <- vmperception(data$m)
    data$diff <- data$p - data$p0
    data$rdiff <- (data$p - data$p0)/data$p0
    expect_lt(mean(abs(data$rdiff)), 0.021)
    expect_true(all(abs(data$rdiff) < 0.096))

    # test first derivation
    f <- function(m) {
        vmperception(m, deriv = TRUE)
    }
    p <- vmperception(4.0, deriv = FALSE) -
        vmperception(1.0, deriv = FALSE)
    expect_equal(p, stats::integrate(f, 1.0, 4.0)$value)

    # test log first derivation
    f <- function(m) {
        vmperception(m, deriv = TRUE, log = TRUE)
    }
    p <- vmperception(4.0, deriv = FALSE, log = TRUE) -
        vmperception(1.0, deriv = FALSE, log = TRUE)
    expect_equal(p, stats::integrate(f, 1.0, 4.0)$value)
})
