test_that("vmperception", {
    # test approximation of `KOSREN90` perceptions published in
    # J. Rendtel, Handbook for meteor observers 2022 edition, p. 149, (c) International Meteor Organization
    data <- data.frame(
        m = seq(-0.4, 8.0, 0.2),
        p = c(
            0.00046, 0.0011, 0.0023, 0.0046, 0.0081, 0.0122, 0.0182,
            0.0257, 0.03444, 0.04365, 0.05495, 0.06918, 0.08511, 0.10351,
            0.13, 0.16, 0.2, 0.24, 0.29, 0.34674, 0.4, 0.46, 0.52, 0.57544,
            0.63096, 0.67764, 0.71, 0.74, 0.77, 0.79, 0.81, 0.83, 0.85,
            0.87, 0.89, 0.91, 0.93, 0.94, 0.96, 0.98, 1.0, 1.0, 1.0
        )
    )

    # ignore lower magnitudes
    data <- subset(data, data$m >= 0)

    # check perceptions
    p <- vmperception(data$m)
    expect_true(all(p < 1.0))
    expect_true(all(p > 0.0))
    rdiff <- (p - data$p)/data$p
    expect_lt(mean(abs(rdiff)), 0.084)
    expect_true(all(abs(rdiff) < 0.221))

    # test first derivation
    f <- function(m) {
        vmperception(m, deriv = TRUE)
    }
    p <- vmperception(4.0) - vmperception(1.0)

    expect_true(abs(p - stats::integrate(f, 1.0, 4.0)$value) < 1e-15)
})

# Approximation of perception probabilities
if(FALSE) {
    # library(ggplot2)

    data <- data.frame(
        m = seq(-0.4, 9.0, 0.2),
        p = c(
            0.00046, 0.0011, 0.0023, 0.0046, 0.0081, 0.0122, 0.0182,
            0.0257, 0.03444, 0.04365, 0.05495, 0.06918, 0.08511, 0.10351,
            0.13, 0.16, 0.2, 0.24, 0.29, 0.34674, 0.4, 0.46, 0.52, 0.57544,
            0.63096, 0.67764, 0.71, 0.74, 0.77, 0.79, 0.81, 0.83, 0.85,
            0.87, 0.89, 0.91, 0.93, 0.94, 0.96, 0.98, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0,1.0, 1.0
        )
    )

    # Model I: g(m) = exp(exp(a * m + b)/a)
    #
    # g'(m) = exp(exp(a * m + b)/a + a * m + b)
    # g'(m)/g(m) = exp(a * m + b)
    # log(g'(m)/g(m)) = a * m + b

    coef.model1 <- with(new.env(), {
        opt.data <- subset(data, data$m > 0.5 & data$p < 1.0)
        opt.fun <- function(params) {
            a <- params[1]
            b <- params[2]
            x <- opt.data$m
            y <- opt.data$p
            y0 <- exp(a * x + b)/a
            #print(data.frame(x = x, y = y, y0 = y0))

            sum((log(y)-y0)^2.0)
        }

        optim(c(1, 0), opt.fun)$par
    })
    print(paste(c('coef model I:', paste(coef.model1, collapse = ', '))))

    # print g of model I
    if (TRUE) {
        with(new.env(), {
            plot.data <- data.frame(
                m = data$m,
                p = data$p,
                g = exp(exp(coef.model1[1] * data$m + coef.model1[2])/coef.model1[1])
            )
            print(plot.data)

            p <- ggplot(plot.data) +
                theme_bw() +
                geom_point(aes(x = m, y = p), color='black') +
                geom_line(aes(x = m, y = g), color='blue') +
                # scale_y_continuous(
                #     trans = "logit",
                #     breaks = c(0.001,0.005,0.01,0.05,0.1,0.5, 0.9, 0.95, 0.98)
                # ) +
                scale_y_continuous(
                    trans = "log",
                    breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.7, 0.95)
                ) +
                scale_x_continuous(
                    name = "m",
                    breaks = seq(-0.5, 8.0, 0.5)
                )
            print(p)
        })
    }

    coef.model2 <- with(new.env(), {
        data <- subset(data, data$p < 1.0)
        g.model1 <- exp(exp(coef.model1[1] * data$m + coef.model1[2])/coef.model1[1])
        #g.model1.weights <- rep(0.0001, nrow(data))
        g.model1.weights <- rep(3, nrow(data))

        lm.data <- data.frame(
            m = c(data$m, data$m),
            p = c(data$p, g.model1),
            w = c(rep(1, nrow(data)), g.model1.weights)
        )
        lm.data$y <- log(lm.data$p/(1 - lm.data$p)) # logit
        res <- stats::lm(y ~ poly(m, degree = 5, raw = TRUE), data = lm.data, weights = lm.data$w)
        coef(res)
    })
    print(paste(c('coef model II:', paste(coef.model2, collapse = ', '))))

    # round
    coef.model2 <- c(-6.24, 3.33, -1.03, 0.241, -0.03, 0.0015)
    print(paste(c('coef model II (rounded):', paste(coef.model2, collapse = ', '))))

    deriv.polynomial <- function(poly.coef, degree) {
        if (0 == degree)
            return(poly.coef)

        if (1 == length(poly.coef))
            return(0)

        poly.coef <- poly.coef[-1]
        deriv.polynomial(poly.coef * seq(along = poly.coef), degree - 1)
    }

    f.polynomial <- function(m, poly.coef) {
        margin.table(poly.coef * t(outer(m, seq(along=poly.coef) - 1, "^")), 2)
    }

    # print g of model II
    with(new.env(), {
        g.model1 <- exp(exp(coef.model1[1] * data$m + coef.model1[2])/coef.model1[1])
        g.model2 <- 1/(1 + exp(-f.polynomial(data$m, coef.model2)))
        plot.data <- data.frame(
            m = data$m,
            p = c(data$p),
            # g.model1 = g.model1,
            g.model2 = g.model2
        )
        print(plot.data)

        p <- ggplot(plot.data) +
            theme_bw() +
            geom_point(aes(x = m, y = p), color='black') +
            # geom_point(aes(x = m, y = g.model1), color='blue') +
            geom_line(aes(x = m, y = g.model2), color='red') +
            # scale_y_continuous(
            #     trans = "logit",
            #     breaks = c(0.001,0.005,0.01,0.05,0.1,0.5, 0.9, 0.95, 0.98)
            # ) +
            scale_y_continuous(
                trans = "log",
                breaks = c(0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0)
            ) +
            scale_x_continuous(
                name = "m",
                breaks = seq(-0.5, 9.0, 0.5)
            )
        print(p)
    })

    # print q of model II
    with(new.env(), {
        m <- seq(-0.4, 9.0, 0.2)
        g0 <- 1/(1 + exp(-f.polynomial(m, coef.model2)))
        inner0 <- f.polynomial(m, coef.model2)
        inner1 <- f.polynomial(m, deriv.polynomial(coef.model2, 1))
        exp.inner0 <- exp(-inner0)
        g1 <- inner1 * exp.inner0/(exp.inner0 + 1)^2
        q <- g1/g0

        plot.data <- data.frame(
            m = m,
            q = q
        )

        p <- ggplot(plot.data) +
            theme_bw() +
            geom_line(aes(x = m, y = q), color='blue') +
            scale_y_continuous(
                trans = "log",
                #limits = c(0.001, 10),
                breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 4)
            ) +
            scale_x_continuous(
                name = "m",
                breaks = seq(-0.6, 9.0, 0.5)
            )
        print(p)
    })
}