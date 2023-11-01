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
    expect_lt(mean(abs(rdiff)), 0.0672)
    expect_true(all(abs(rdiff) < 0.2076))

    # test first derivation
    f <- function(m) {
        vmperception(m, deriv = TRUE)
    }
    p <- vmperception(4.0) - vmperception(1.0)
    expect_true(abs(p - stats::integrate(f, 1.0, 4.0)$value) < 1e-14)
})

# Approximation of perception probabilities
if(FALSE) {
    #library(ggplot2)

    deriv.polynomial <- function(poly.coef, degree = 1L) {
        if (0L == degree)
            return(poly.coef)

        if (1 == length(poly.coef))
            return(0)

        exponents <- as.numeric(names(poly.coef))
        poly.coef <- poly.coef * exponents
        if (0L %in% exponents) {
            intercept.idx <- 0L == exponents
            poly.coef <- poly.coef[! intercept.idx]
            exponents <- exponents[! intercept.idx]
        }
        exponents <- exponents - 1L
        names(poly.coef) <- exponents

        deriv.polynomial(poly.coef, degree - 1L)
    }

    f.polynomial <- function(m, poly.coef) {
        exponents <- as.numeric(names(poly.coef))
        margin.table(poly.coef * t(outer(m, exponents, "^")), 2)
    }

    data <- data.frame(
        m = seq(-0.4, 9.0, 0.2),
        p = c(
            0.00046, 0.0011, 0.0023, 0.0046, 0.0081, 0.0122, 0.0182,
            0.0257, 0.03444, 0.04365, 0.05495, 0.06918, 0.08511, 0.10351,
            0.13, 0.16, 0.2, 0.24, 0.29, 0.34674, 0.4, 0.46, 0.52, 0.57544,
            0.63096, 0.67764, 0.71, 0.74, 0.77, 0.79, 0.81, 0.83, 0.85,
            0.87, 0.89, 0.91, 0.93, 0.94, 0.96, 0.98, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0
        )
    )

    # Model 0: g(m) = exp(exp(a * m + b)/a)
    #
    # g'(m) = exp(exp(a * m + b)/a + a * m + b)
    # g'(m)/g(m) = exp(a * m + b)
    # log(g'(m)/g(m)) = a * m + b

    coef.model0 <- with(new.env(), {
        opt.data <- subset(data, data$m > 0.5 & data$p < 1.0)
        opt.fun <- function(params) {
            a <- params[1]
            b <- params[2]
            x <- opt.data$m
            y <- opt.data$p
            y0 <- exp(a * x + b)/a

            sum((log(y)-y0)^2.0)
        }

        optim(c(1, 0), opt.fun)$par
    })
    print(paste(c('coef model 0:', paste(coef.model0, collapse = ', '))))

    # print g of model 0
    if (TRUE) {
        with(new.env(), {
            plot.data <- data.frame(
                m = data$m,
                p = data$p,
                g = exp(exp(coef.model0[1] * data$m + coef.model0[2])/coef.model0[1])
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

    coef.model <- with(new.env(), {
        data <- subset(data, data$p < 1.0)

        g.model0 <- exp(exp(coef.model0[1] * data$m + coef.model0[2])/coef.model0[1])
        g.model0.weights <- rep(0.01, nrow(data))

        data <- data.frame(
            x = c(data$m + 0.5, data$m + 0.5),
            p = c(data$p, g.model0),
            w = c(rep(1, nrow(data)), g.model0.weights)
        )

        coef.model <- with(new.env(), {
            p.fun <- function(x, params) {
                names(params) <- seq(along = params) # exponents
                inner0 <- f.polynomial(x, params)
                1.0 - exp(-inner0)
            }

            ll <- function(params) {
                #sum(data$w * ((data$p - p.fun(data$x, params)))^2)
                sum(data$w * ((log(data$p) - log(p.fun(data$x, params))))^2)
            }

            optim(c(0.004, 0.0012, 0.0035, 0.0007), ll)$par
        })

        names(coef.model) <- seq(along = coef.model) # exponents
        coef.model
    })
    print(paste(c('coef model:', paste(coef.model, collapse = ', '))))

    # round
    coef.model <- c(0.004, 0.0012, 0.0035, 0.0007)
    names(coef.model) <- seq(along = coef.model) # exponents
    print(paste(c('coef model (rounded):', paste(coef.model, collapse = ', '))))

    # print g
    if (TRUE) {
        with(new.env(), {
            inner0 <- f.polynomial(data$m + 0.5, coef.model)
            g.model <- 1.0 - exp(-inner0)

            plot.data <- data.frame(
                m = data$m,
                p = data$p,
                g.model = g.model
            )
            print(plot.data)

            p <- ggplot(plot.data) +
                theme_bw() +
                geom_point(aes(x = m, y = p), color='black') +
                geom_line(aes(x = m, y = g.model), color='red') +
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
    }

    # print g'
    if (TRUE) {
        with(new.env(), {
            m <- c(-0.499, seq(-0.45, 8.5, 0.05))
            # d/dx(1 - exp(-f(x))) = e^(-f(x)) f'(x)
            inner0 <- f.polynomial(m + 0.5, coef.model)
            inner1 <- f.polynomial(m + 0.5, deriv.polynomial(coef.model, degree = 1L))
            g.model <- exp(-inner0) * inner1

            plot.data <- data.frame(
                m = m,
                g.model = g.model
            )

            p <- ggplot(plot.data) +
                theme_bw() +
                geom_line(aes(x = m, y = g.model), color='red') +
                scale_y_continuous(
                    trans = "log",
                    breaks = c(0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2)
                ) +
                scale_x_continuous(
                    name = "m",
                    breaks = seq(-0.5, 9.0, 0.5)
                )
            print(p)
        })
    }

    # print q
    if (TRUE) {
        with(new.env(), {
            m <- seq(-0.4, 8.4, 0.05)
            inner0 <- f.polynomial(m + 0.5, coef.model)
            g0 <- 1.0 - exp(-inner0)
            inner1 <- f.polynomial(m + 0.5, deriv.polynomial(coef.model, degree = 1L))
            g1 <- exp(-inner0) * inner1
            q <- g1/g0

            plot.data <- data.frame(
                m = m,
                q = q
            )

            p <- ggplot(plot.data) +
                theme_bw() +
                geom_line(aes(x = m, y = q), color='blue') +
                scale_x_continuous(
                    name = "m",
                    limits = c(-0.5, 8.5),
                    breaks = seq(-0.5, 8.5, 0.5)
                ) +
                scale_y_continuous(
                    trans = "log",
                    # labels = scales::comma,
                    limits = c(0.001, 20),
                    breaks = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
                )
            print(p)
        })
    }

    if (TRUE) {
        with(new.env(), {
            r <- seq(1.2, 3.6, 0.1)
            limmag <- seq(5.6, 6.5, 0.2)
            m <- seq(-100, 6, 1)

            g.fun <- function(m) {
                idx <- m >= -0.400001 & m <= 9.0
                p <- rep(1.0, length(m))
                p[m < -0.400001] <- 0.0
                if (any(idx)) {
                    x <- m[idx]
                    p[idx] <- data$p[match(as.integer(round(10*x)), as.integer(round(10*data$m)))]
                }

                p
            }

            result <- do.call(
                rbind.data.frame,
                sapply(r, function(r) {
                    obs <- do.call(
                        rbind.data.frame,
                        sapply(limmag, function(limmag) {
                            if (TRUE) {
                                p <- dgeom(6 - m, 1.0 - 1.0/r) * g.fun(limmag - m)
                                p <- p/sum(p)
                            } else {
                                p <- dvmgeom(m, limmag, r)
                            }

                            list(
                                limmag = rep(limmag, length(m)),
                                m = m,
                                p = p,
                                q = p * vmperception(limmag - m, TRUE)/vmperception(limmag - m)
                            )
                        }, simplify = FALSE)
                    )

                    r.approx <- with(obs,{
                        # log likelihood function
                        llr <- function(r) {
                            -sum(p * dvmgeom(m, limmag, r, log=TRUE))
                        }

                        # maximum likelihood estimation (MLE) of r
                        optim(2, llr, method='Brent', lower=1.1, upper=5, hessian=FALSE)
                    })

                    list(
                        r = r,
                        r.est = r.approx$par,
                        q.exp = exp(sum(obs$q)/length(limmag))
                    )
                }, simplify = FALSE)
            )

            print(result)
        })
    }
}
