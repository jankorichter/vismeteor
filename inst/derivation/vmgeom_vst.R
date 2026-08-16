##
## Derivation of the variance-stabilizing transformation for the geometric
## model of visual meteor magnitudes (vmgeom). Developer notes -- not part of
## the user documentation.
##
## Method
##   g = f(limmag - m - 1) / f(limmag - m)     E[g] = 1/r, exact, free of limmag
##   t = a * ((g + s)^lambda - s^lambda) / lambda    chosen so that Var[t] ~ 1
##   log(r) = c + d * log(E[t]) + e * E[t]     back-transformation,
##                                             constrained to r(tm_max) = 1
##
## g is the statistic of the rate-based estimator ("By Rate" in
## vignette("vmgeom")). (a, lambda, s) stabilise its variance, (c, d, e) invert
## the result. (c, d, e) must be FITTED, not obtained by rearranging: the
## estimator averages transformed magnitudes, and E[h(g)] != h(E[g]).
##
## Workflow
## 1. Run top-to-bottom. The `if (FALSE)` block re-optimises (a, lambda, s);
##    set it to TRUE to determine the parameters anew.
## 2. Copy (a, lambda, s) and the printed (c, d, e) into `.vmgeom_vst_params`
##    in `R/vmgeom_vst.R`.
## 3. Re-run devtools::load_all() and tests/testthat/test_vmgeom_vst.R.
##
## Fix and round (a, lambda, s) BEFORE fitting (c, d, e), so that the two
## halves belong to each other.
##
## Choice of the back-transformation (max rel. error in r over the grid)
##
##     c + d log t                    2.09%   r(tm_max) 1.040
##     c + d log t + e (log t)^2      0.95%   turns at r ~ 65, then r -> 0
##     c + d log t + e t              0.95%   r(tm_max) 0.98
##     c + d log t + e t, constrained 1.21%   r(tm_max) 1.000     <- in use
##
##   The quadratic in log(t) is excluded although it fits best: it reverses
##   inside the attainable range (0, tm_max]. The term in `t` keeps the
##   curvature without that -- it vanishes as t -> 0, and d/t + e < 0 for
##   d < 0, e < 0 makes monotonicity structural (hence the stopifnot below).
##   The constraint is exact, not cosmetic: g = 1 means r = 1.
##   A cubic term adds nothing (0.98%); the rest is the limmag spread.
##
## Why the shifted Box-Cox family and not a plain power t = a * g^b
##   Re-optimising the plain power gives (4.518, 0.677), i.e. its old values,
##   with mean((Var[t]-1)^2) = 0.00528 and max|Var[t]-1| = 0.190, against
##   0.00411 and 0.160 here. The plain power is the case s -> 0, lambda = b.
##
## Why not MASS::boxcox()
##   1. It needs y > 0, but P[g = 0] > 0 (1.8% at r = 2 / limmag 6.0, 20.6% at
##      r = 4 / limmag 6.4). Dropping the zeros breaks E[g] = 1/r: the
##      conditional mean is 0.3148 instead of 0.2500 at r = 4, limmag 6.4.
##   2. It maximises a normal likelihood for ONE sample, not variance constancy
##      across r. lambda drifts accordingly (limmag 6.0, zeros removed):
##      r = 1.4 -> 1.38, r = 2 -> 0.56, r = 3 -> 0.24, r = 4 -> -0.02.
##
## On re-running the optimisation
##   The optimum is unique but the valley is flat -- lambda and s are partly
##   interchangeable, so different starting points reach the same score with
##   parameters differing in the third decimal. Check the SCORE, not the
##   parameters. The starting point is therefore fixed and neutral rather than
##   the result itself.
##
## Why no per-limmag parameter table (unlike the pre-3.1 approximation)
##   Fitting per limmag improves the variance but destroys the common scale
##   across limmag, and the back-transformation then degrades badly.
##
## On the calibrated range 1.4 <= r <= 4
## - Downwards nothing bounds it; only Var[t] deteriorates as r -> 1.
## - Upwards the spread across limmag dominates: E[g] is limmag-free, E[t] is
##   not, and one set (c, d, e) must serve every limmag.
## - |Var[t] - 1| reaches 0.16 at the edges of the range.
## - E[t] is strictly decreasing in r out to r = 30 for every limmag, which is
##   why vmgeom_vst_to_r() needs no calibration window.
##
## Caveat when checking small r: the geometric weights decay like r^m, so a
## grid ending at m = -200 truncates enough mass near r = 1 to fake a bend that
## is not there. Pass a wider `m` to data.t.fun().
##
# Approximation of the Variance-stabilizing transformation
library(vismeteor)
library(ggplot2)

r <- seq(1.4, 4.0, 0.1)
m <- seq(-200, 6, 1)
limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)

# c(a, lambda, s); keep in sync with `.vmgeom_vst_params` in `R/vmgeom_vst.R`
params <- c(3.0086, 0.3714, 0.1731)

# starting point for the optimisation below; deliberately not `params`
# (see header), these are the old plain-power values plus a small shift
params.start <- c(4.52, 0.68, 0.05)

# t = a * ((g + s)^lambda - s^lambda) / lambda, vectorised over m / limmag.
# The shift lets g = 0 through, subtracting s^lambda keeps t(0) = 0.
myVmgeomVstFromMagn <- function(m, limmag, params) {
    a <- params[1]
    lambda <- params[2]
    s <- params[3]

    # 0 outside the support: the denominator vanishes there and 0 * NaN would
    # poison the sum(p * t) below
    den <- vismeteor::vmperception(limmag - m)
    g <- ifelse(den > 0, vismeteor::vmperception(limmag - m - 1) / den, 0)

    a * ((g + s)^lambda - s^lambda) / lambda
}

# Moments of m, g and t under vmgeom over a grid of r and limmag.
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
    # Estimate (a, lambda, s) by minimizing the deviation of Var[t] from 1.
    #
    # p and g depend only on the grid, so they are precomputed here rather
    # than inside the objective -- seconds instead of minutes.
    #
    # L-BFGS-B with real bounds, not Nelder-Mead: unbounded, the optimiser
    # walks into negative lambda, where t.mean turns negative and log() fails.
    params <- with(new.env(), {
        model <- expand.grid(r = r, limmag = limmag)
        p.grid <- mapply(
            function(r, limmag) dvmgeom(m, limmag, r),
            model$r, model$limmag, SIMPLIFY = FALSE
        )
        g.grid <- lapply(model$limmag, function(limmag) {
            den <- vismeteor::vmperception(limmag - m)
            ifelse(den > 0, vismeteor::vmperception(limmag - m - 1) / den, 0)
        })

        t.var.fun <- function(params) {
            vapply(seq_along(p.grid), function(i) {
                t <- params[1] *
                    ((g.grid[[i]] + params[3])^params[2] - params[3]^params[2]) /
                    params[2]
                t.mean <- sum(p.grid[[i]] * t)
                sum(p.grid[[i]] * (t - t.mean)^2)
            }, numeric(1))
        }

        opt.res <- optim(
            params.start,
            fn = function(params) mean((t.var.fun(params) - 1.0)^2),
            method = "L-BFGS-B",
            lower = c(0.1, 0.01, 1e-04),
            upper = c(20.0, 3.0, 1.0),
            control = list(factr = 1e07, maxit = 1000L)
        )
        # a silent non-convergence would be indistinguishable from a result
        stopifnot(0L == opt.res[['convergence']])

        # the score, not the parameters, is what a re-run has to reproduce
        print(opt.res[['par']])
        print(opt.res[['value']])
        round(opt.res[['par']], 4)
    })
    dput(params)

    # (c, d, e) below must be refitted from these rounded values
    data.t <- data.t.fun(params, r = r, limmag = limmag)
}

# upper bound of t, reached at g = 1; `a` is not that bound itself, so
# `R/vmgeom_vst.R` keeps it as a separate entry `tm_max`
tm.max <- params[1] * ((1.0 + params[3])^params[2] - params[3]^params[2]) / params[2]
print(round(tm.max, 6))

# Back-transformation, subject to r(tm_max) = 1, i.e.
# c + d * log(tm_max) + e * tm_max = 0. Eliminating c turns the constraint
# into a fit through the origin.
lm.res.t <- lm(
    log(r) ~ I(log(t.mean) - log(tm.max)) + I(t.mean - tm.max) - 1,
    data = data.t
)
lm.res.t.coef <- c(
    c = -coef(lm.res.t)[[1]] * log(tm.max) - coef(lm.res.t)[[2]] * tm.max,
    d = coef(lm.res.t)[[1]],
    e = coef(lm.res.t)[[2]]
)
# Paste into `.vmgeom_vst_params` in `R/vmgeom_vst.R` (entries c, d and e).
print(round(lm.res.t.coef, 6))

# both negative <=> r strictly decreasing on all of (0, tm_max] (see header)
stopifnot(lm.res.t.coef[['d']] < 0.0, lm.res.t.coef[['e']] < 0.0)

# max absolute and relative error on the r scale
data.t$r.pred <- exp(
    lm.res.t.coef[['c']] +
        lm.res.t.coef[['d']] * log(data.t$t.mean) +
        lm.res.t.coef[['e']] * data.t$t.mean
)
print(max(abs(data.t$r - data.t$r.pred)))
print(max(abs(data.t$r - data.t$r.pred) / data.t$r))


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
    # Plotted against log(E[g]), not log(r): E[g] = 1/r is exact and
    # limmag-free, so every deviation seen here comes from the transformation.
    # The fit uses 1.4 <= r <= 4, the rest is extrapolation.
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

        # a curve in these coordinates, hence a line layer, not geom_abline
        data.wide$fit <- -(
            lm.res.t.coef[['c']] +
                lm.res.t.coef[['d']] * log(data.wide$t.mean) +
                lm.res.t.coef[['e']] * data.wide$t.mean
        )
        data.wide$err <- data.wide$fit - log(data.wide$g.mean)

        # one line per limmag, dashed the fit
        p <- ggplot(data.wide, aes(x = log(t.mean), y = log(g.mean), group = limmag)) +
            theme_bw() +
            annotate(
                'rect',
                xmin = min(log(data.cal$t.mean)), xmax = max(log(data.cal$t.mean)),
                ymin = -Inf, ymax = Inf, alpha = 0.08
            ) +
            geom_line(color = 'blue', alpha = 0.5) +
            geom_line(
                aes(y = fit),
                linetype = 'dashed', color = 'red'
            ) +
            scale_x_continuous(name = "log(E[t])") +
            scale_y_continuous(name = "log(E[g]) = -log(r)") +
            ggtitle("log(E[g]) vs log(E[t]); shaded: calibrated range 1.4 <= r <= 4")
        print(p)
    })
}
