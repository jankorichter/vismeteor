##
## Derivation of the variance-stabilizing transformation for the geometric
## model of visual meteor magnitudes (vmgeom). Developer notes -- not part of
## the user documentation.
##
## Method
##   g = f(limmag - m - 1) / f(limmag - m)     E[g] = 1/r, exact, free of limmag
##   t = a * g^b                               chosen so that Var[t] ~ 1
##   log(r) = c + d * log(E[t])                back-transformation
##
## g is the statistic of the rate-based estimator ("By Rate" in
## vignette("vmgeom")). (a, b) stabilise its variance, (c, d) invert the result.
##
## (c, d) must be FITTED, not obtained by rearranging t = a * g^b: the estimator
## averages transformed magnitudes, and E[g^b] != E[g]^b for b != 1. The
## algebraic inverse would stay wrong at any sample size.
##
## Workflow
## 1. Run top-to-bottom. The `if (FALSE)` block re-optimises (a, b) and takes
##    several minutes; the values in use are hard-coded.
## 2. Copy (a, b) and the printed (c, d) into `.vmgeom_vst_params` in
##    `R/vmgeom_vst.R`.
## 3. Re-run devtools::load_all() and tests/testthat/test_vmgeom_vst.R.
##
## Why no per-limmag parameter table (unlike the pre-3.1 approximation)
##   limmag enters through the perception probabilities alone. Fitting (a, b)
##   per limmag improves the variance but destroys the common scale across
##   limmag, and the back-transformation then degrades badly.
##
## On the calibrated range 1.4 <= r <= 4
## - Downwards nothing bounds it. Both error terms below shrink as r -> 1;
##   only Var[t] deteriorates there.
## - Upwards both terms grow, and the SPREAD across limmag at fixed r
##   dominates: E[g] is limmag-free exactly, E[g^b] is not, and one pair (c, d)
##   must serve every limmag. Spread in log(r) is 0.006 at r = 2, 0.028 at
##   r = 3.3, 0.096 at r = 6.4, 0.20 at r = 12; the part common to all limmag
##   stays under 0.012 only out to r ~ 4.3 and reaches 0.14 at r = 12.
## - Var[t] is the two-sided limit and what the range actually describes:
##   |Var[t] - 1| is 0.19 at r = 4, 0.36 at r = 6.3, and below 1.4 it decays
##   faster still (0.36 at r = 1.2), since nearly all meteors then fall into
##   the brightest classes.
## - Monotonicity never breaks: E[t] is strictly decreasing in r out to r = 30
##   for every limmag. This is why vmgeom_vst_to_r() needs no calibration
##   window and returns NA only for values the model cannot produce.
## - Re-optimising (a, b) does not help: the optimum on [1.4, 4] is
##   (4.502, 0.681), i.e. where we already are. On [1.2, 5] it is (4.577,
##   0.692), trading the lower end for very little at the upper one.
##
## Caveat when checking small r: the geometric weights decay like r^m, so a
## magnitude grid ending at m = -200 truncates enough mass near r = 1 to fake
## a bend that is not there. Pass a wider `m` to data.t.fun().
##
# Approximation of the Variance-stabilizing transformation
library(vismeteor)
library(ggplot2)

r <- seq(1.4, 4.0, 0.1)
m <- seq(-200, 6, 1)
limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)

# dput(params)
# NOTE: Keep these values in sync with `.vmgeom_vst_params` inside
# `R/vmgeom_vst.R`. After recalibration, copy the rounded values there so the
# package picks them up without sourcing this script.
params <- c(4.52, 0.68)

# t = a * g^b for vectorised m / limmag; params = c(a, b).
# g lies in [0, 1], so a is also the upper bound of t.
myVmgeomVstFromMagn <- function(m, limmag, params) {
    a <- params[1]
    b <- params[2]

    # 0 outside the support: the denominator vanishes there and 0 * NaN would
    # poison the sum(p * t) below
    den <- vismeteor::vmperception(limmag - m)
    g <- ifelse(den > 0, vismeteor::vmperception(limmag - m - 1) / den, 0)

    a * g^b
}

# Moments of m, g and t under vmgeom over a grid of r and limmag.
# Returns r, q, limmag, m.mean, m.var, g.mean, g.var, t.mean, t.var.
# g.var vs t.var is the diagnostic: g.var varies strongly with r, t.var must not.
data.t.fun <- function(params, r, limmag, m = seq(-200, 6, 1)) {
    model <- expand.grid(r = r, limmag = limmag)
    model <- do.call(
        rbind.data.frame,
        mapply(function(r, limmag) {
            p <- dvmgeom(m, limmag, r)

            m.mean <- sum(p * (limmag - m))
            m.var <- sum(p * (limmag - m - m.mean)^2)

            # the rate-based estimator: E[g] = 1/r exactly
            den <- vismeteor::vmperception(limmag - m)
            g <- ifelse(den > 0, vismeteor::vmperception(limmag - m - 1) / den, 0)
            g.mean <- sum(p * g)
            g.var <- sum(p * (g - g.mean)^2)

            t <- myVmgeomVstFromMagn(m, limmag, params)
            t.mean <- sum(p * t)
            t.var <- sum(p * (t - t.mean)^2)

            list(
                r = r,
                q = log(r),
                limmag = limmag,
                m.mean = m.mean,
                m.var = m.var,
                g.mean = g.mean,
                g.var = g.var,
                t.mean = t.mean,
                t.var = t.var
            )
        }, model$r, model$limmag , SIMPLIFY = FALSE)
    )
}
data.t <- data.t.fun(params, r = r, limmag = limmag)

if (FALSE) {
    # Estimate (a, b) by minimizing the deviation of the t-variance
    # from a target value (here 1.0)
    params <- with(new.env(), {
        opt.res <- optim(params, fn = function(params) {
            df <- data.t.fun(params, r = r, limmag = limmag)
            mean((df$t.var - 1.0)^2)
        })
        print(opt.res[['par']])
        round(opt.res[['par']], 2)
    })
    dput(params)
}

# Back-transformation (c, d): map E[t] onto E[g], fitted, not rearranged
# (see header). E[g] = 1/r is exact, so log(E[g]) is the response.
lm.res.t <- lm(log(g.mean) ~ log(t.mean), data = data.t)
print(summary(lm.res.t)$r.squared)

# log(r) = -log(g) = -(k0 + k1 * log(t)) = c + d * log(t)
lm.res.t.coef <- c(c = -coef(lm.res.t)[[1]], d = -coef(lm.res.t)[[2]])
# Paste into `.vmgeom_vst_params` in `R/vmgeom_vst.R` (entries c and d).
print(round(lm.res.t.coef, 6))

# max error on the r scale; the remainder is the limmag spread and does not
# go away with higher polynomial degrees
data.t$r.pred <- exp(lm.res.t.coef[['c']] + lm.res.t.coef[['d']] * log(data.t$t.mean))
print(max(abs(data.t$r - data.t$r.pred)))


if (TRUE) {
    with(new.env(), {
        m <- seq(7.0, -4.5, -0.1)
        limmag <- seq(5.6, 6.5, 0.1)

        p <- ggplot(data.t) + theme_bw()

        for (l in limmag) {
            # shape of the transform t versus x = limmag - m
            plot.data <- expand.grid(limmag=limmag, m=m)
            plot.data$dm <- plot.data$limmag - plot.data$m
            plot.data <- plot.data[order(plot.data$dm),]
            plot.data$x <-  myVmgeomVstFromMagn(m = plot.data$m, limmag = plot.data$limmag, params = params)
            plot.data <- subset(plot.data, plot.data$x > 0)
            p <- p + geom_line(aes(x = dm, y = x), color='blue', data = subset(plot.data, l == plot.data$limmag))
        }

        p <- p +
            scale_x_continuous(
                name = "m",
                breaks = seq(0, 11, 1.0)
            ) +
            scale_y_continuous(
                name = "t",
                breaks = seq(0, 5, 0.5)
            )
        print(p)
    })
}

if (TRUE) {
    # Plot 2: how far does the back-transformation reach?
    #
    # Plotting log(E[g]) rather than log(r) is the point: E[g] = 1/r is exact
    # and limmag-free, so every deviation seen here is introduced by t = a * g^b
    # alone. The fit uses 1.4 <= r <= 4 only, so the rest is real extrapolation.
    with(new.env(), {
        # wider `m` than the default: see the caveat on small r in the header
        r.wide <- exp(seq(log(1.1), log(10.0), length.out = 90))
        data.wide <- data.t.fun(
            params,
            r = r.wide, limmag = limmag,
            m = seq(-600, 8, 1)
        )

        # sanity check: E[g] must reproduce 1/r to numerical accuracy
        stopifnot(max(abs(data.wide$g.mean - 1 / data.wide$r)) < 1e-04)

        data.cal <- subset(data.wide, data.wide$r >= 1.4 & data.wide$r <= 4.0)
        cf <- coef(lm(log(g.mean) ~ log(t.mean), data = data.cal))
        data.wide$err <- (cf[[1]] + cf[[2]] * log(data.wide$t.mean)) -
            log(data.wide$g.mean)

        # one line per limmag, dashed the fit; straight over the whole range,
        # so linearity is not what bounds the calibrated window
        p <- ggplot(data.wide, aes(x = log(t.mean), y = log(g.mean), group = limmag)) +
            theme_bw() +
            annotate(
                'rect',
                xmin = min(log(data.cal$t.mean)), xmax = max(log(data.cal$t.mean)),
                ymin = -Inf, ymax = Inf, alpha = 0.08
            ) +
            geom_line(color = 'blue', alpha = 0.5) +
            geom_abline(
                intercept = cf[[1]], slope = cf[[2]],
                linetype = 'dashed', color = 'red'
            ) +
            scale_x_continuous(name = "log(E[t])") +
            scale_y_continuous(name = "log(E[g]) = -log(r)") +
            ggtitle("log(E[g]) vs log(E[t]); shaded: calibrated range 1.4 <= r <= 4")
        print(p)
    })
}

