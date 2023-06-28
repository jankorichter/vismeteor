#' @name vmgeom
#' @aliases dvmgeom
#' @aliases pvmgeom
#' @aliases qvmgeom
#' @aliases rvmgeom
#' @title Visual magnitude distribution of geometric distributed meteor magnitudes
#' @description
#' Density, distribution function, quantile function and random generation for the
#' visual magnitude distribution of geometric distributed meteor magnitudes.
#' @param r numeric; the population index. It is the only parameter of the distribution.
#' @param m numeric; difference between the limiting magnitude and the meteor magnitude.
#' @param lm numeric; limiting magnitude (optional).
#' @param p numeric; probability.
#' @param n numeric; count of meteor magnitudes.
#' @param log logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE` (default) probabilities are
#'     \eqn{P[M \le m]}, otherwise, \eqn{P[M > m]}.
#' @param perception.fun function; perception probability function (optional).
#'     Default is [vismeteor::vmperception].
#' @details
#' In visual meteor observation, it is common to estimate meteor magnitudes in integer values.
#' Hence, this distribution is discrete and has the density
#' \deqn{
#'     {\displaystyle P[M = m] \sim f(m) \, \mathrm r^{-m}} \,\mathrm{,}
#' }
#' where \eqn{m \ge -0.5} is the difference between the limiting magnitude and the meteor magnitude and
#' \eqn{f(m)} is the perception probability.
#' This distribution is thus a convolution of the
#' [perception probabilities][vismeteor::vmperception] with the
#' actual [geometric distribution][stats::Geometric] of the meteor magnitudes.
#' Therefore, the parameter `p` of the geometric distribution is `p = 1 - 1/r`.
#'
#' The parameter `lm` indicate what the parameter `m` refers to.
#'
#' If `lm` is passed, then `m` must be an integer meteor magnitude.
#' The length of the vector `lm` must then be equal to the length of the vector `m`
#' or `lm` is a scalar value.
#' In case of `rvmgeom`, the length of the vector `lm` must be `n` or `lm` is a scalar value.
#' Otherwise `m` is the difference between the limiting magnitude and the meteor magnitude.
#'
#' If the perception probabilities function `perception.fun` is given,
#' it must have the signature `function(m, log = FALSE)` and must return the perception probabilities of
#' the difference `m` between the limiting magnitude and the meteor magnitude.
#' If `m >= 15.0`, the `perception.fun` function should return the perception probability of `1.0`.
#' If `log = TRUE` is given, the logarithm value of the perception probabilities
#' must be returned. `perception.fun` is resolved using [match.fun].
#' @return
#' `dvmgeom` gives the density, `pvmgeom` gives the distribution function,
#' `qvmgeom` gives the quantile function, and `rvmgeom` generates random deviates.
#'
#' The length of the result is determined by `n` for `rvmgeom`, and is the maximum
#' of the lengths of the numerical vector arguments for the other functions.
#'
#' Since the distribution is discrete, `qvmgeom` and `rvmgeom` always return integer values.
#' `qvmgeom` can return `NaN` value with a warning.
#' @seealso [vismeteor::vmperception]
#'   [stats::Geometric]
#' @examples
#' N <- 100
#' r <- 2.0
#' limmag <- 6.5
#' (m <- seq(6, -4))
#'
#' # discrete density of `N` meteor magnitudes
#' (freq <- round(N * dvmgeom(m, r, lm=limmag)))
#'
#' # log likelihood function
#' lld <- function(r) {
#'     -sum(freq * dvmgeom(m, r, lm=limmag, log=TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of r
#' est <- optim(2, lld, method='Brent', lower=1.1, upper=4, hessian=TRUE)
#'
#' # estimations
#' est$par # mean of r
#'
#' # generate random meteor magnitudes
#' m <- rvmgeom(N, r, lm=limmag)
#'
#' # log likelihood function
#' llr <- function(r) {
#'     -sum(dvmgeom(m, r, lm=limmag, log=TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of r
#' est <- optim(2, llr, method='Brent', lower=1.1, upper=4, hessian=TRUE)
#'
#' # estimations
#' est$par # mean of r
#' sqrt(1/est$hessian[1][1]) # standard deviation of r

#' @rdname vmgeom
#' @export
dvmgeom <- function(m, r, lm = NULL, log = FALSE, perception.fun = NULL) {
    if (1 == length(r)) {
        r <- rep(r, length(m))
    }

    if (anyNA(m) | anyNA(r) | anyNA(lm) ) {
        stop("NA's are not allowed!")
    }

    if (any(1.0 > r)) {
        stop(paste0('r must be greater than 1.0 instead of "', r, '"!'))
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    m.max <- 15L
    p.geom <- 1.0 - 1.0/r

    if (! is.null(lm)) {
        m = lm - m
    }

    m.round <- round(m)
    offset <- rep(0.0, length(m))
    idx <- is.infinite(m)
    offset[! idx] <- m[! idx] - m.round[! idx]
    m <- m.round

    idx <- -0.5 == offset
    m[idx] <- m[idx] - 1L
    offset[idx] <- offset[idx] + 1.0

    if (1 == length(offset)) {
        offset <- rep(offset, length(m))
    }

    f.norm <- function(m, offset, p.geom) {
        m <- as.integer(seq(0, m.max))
        sum(stats::dgeom(m, p.geom) * perception.fun(m + offset, log = FALSE)) +
            stats::pgeom(m.max, p.geom, lower.tail = FALSE)
    }

    f.density <- function(m, offset, p.geom) {
        idx <- m <= m.max
        d <- stats::dgeom(m, p.geom, log = TRUE)
        d[idx] <- d[idx] + perception.fun(m[idx] + offset, log = TRUE)

        d - base::log(f.norm(m, offset, p.geom))
    }

    arg.data <- data.frame(
        m = m,
        offset = offset,
        p.geom = p.geom
    )

    data.f <- as.factor(paste0(offset, '/', p.geom))
    data.s <- split(arg.data, data.f)
    d <- lapply(data.s, function(data) {
        m <- data$m
        offset <- data$offset[1]
        p.geom <- data$p.geom[1]
        d <- rep(-Inf, length(m))

        idx <- m > -1
        if(any(idx)) {
            d[idx] <- f.density(m[idx], offset, p.geom)
        }

        if (! log) {
            d[idx] <- exp(d[idx])
            d[!idx] <- 0.0
        }

        d
    })

    unsplit(d, data.f)
}

#' @rdname vmgeom
#' @export
pvmgeom <- function(m, r, lm = NULL, lower.tail = TRUE, log = FALSE, perception.fun = NULL) {
    if (1 == length(r)) {
        r <- rep(r, length(m))
    }

    if (anyNA(m) | anyNA(r) | anyNA(lm) ) {
        stop("NA's are not allowed!")
    }

    if (any(1.0 > r)) {
        stop(paste0('r must be greater than 1.0 instead of "', r, '"!'))
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    m.max <- 15L
    p.geom <- 1.0 - 1.0/r

    if (! is.null(lm)) {
        m = lm - m
    }

    m.round <- round(m)
    offset <- rep(0.0, length(m))
    idx <- is.infinite(m)
    offset[! idx] <- m[! idx] - m.round[! idx]
    m <- m.round

    idx <- -0.5 == offset
    m[idx] <- m[idx] - 1L
    offset[idx] <- offset[idx] + 1.0

    if (1 == length(offset)) {
        offset <- rep(offset, length(m))
    }

    f.norm <- function(offset, p.geom){
        m <- as.integer(seq(0, m.max))
        sum(stats::dgeom(m, p.geom) * perception.fun(m + offset, log = FALSE)) +
            stats::pgeom(m.max, p.geom, lower.tail = FALSE)
    }

    f.density <- function(m, offset, p.geom) {
        stats::dgeom(m, p.geom) * perception.fun(m + offset, log = FALSE)
    }

    f.sum <- Vectorize(function(m, offset, p.geom) {
        m <- as.integer(seq(0, m))
        sum(f.density(m, offset, p.geom))
    })

    f.prob <- function(m, offset, p.geom) {
        norm <- f.norm(offset, p.geom)
        p <- rep(0.0, length(m))

        if (lower.tail) {
            idx <- m <= m.max
            if(any(idx)) {
                p[idx] <- f.sum(m[idx], offset, p.geom)/norm
            }

            idx <- m > m.max
            if(any(idx)) {
                p[idx] <- 1.0 - stats::pgeom(m[idx], p.geom, lower.tail = FALSE)/norm
            }
        } else {
            idx <- m <= m.max
            if(any(idx)) {
                p[idx] <- 1 - f.sum(m[idx], offset, p.geom)/norm
            }

            idx <- m > m.max
            if(any(idx)) {
                p[idx] <- stats::pgeom(m[idx], p.geom, lower.tail = FALSE)/norm
            }
        }

        p
    }

    arg.data <- data.frame(
        m = m,
        offset = offset,
        p.geom = p.geom
    )

    data.f <- as.factor(paste0(offset, '/', p.geom))
    data.s <- split(arg.data, data.f)
    p <- lapply(data.s, function(data) {
        m <- data$m
        offset <- data$offset[1]
        p.geom <- data$p.geom[1]

        if (lower.tail) {
            p <- rep(0.0, length(m))
        } else {
            p <- rep(1.0, length(m))
        }

        idx <- m > -1
        if(any(idx)) {
            p[idx] <- f.prob(m[idx], offset, p.geom)
        }

        if (log) {
            p[idx] <- base::log(p[idx])
            if (lower.tail) {
                p[!idx] <- -Inf
            } else {
                p[!idx] <- 0.0
            }
        }

        p
    })

    unsplit(p, data.f)
}

#' @rdname vmgeom
#' @export
qvmgeom <- function(p, r, lm = NULL, lower.tail = TRUE, perception.fun = NULL) {
    if (1 == length(r)) {
        r <- rep(r, length(p))
    }

    if (anyNA(p) | anyNA(r) | anyNA(lm) ) {
        stop("NA's are not allowed!")
    }

    if (any(1.0 > r)) {
        stop(paste0('r must be greater than 1.0 instead of "', r, '"!'))
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    m.max <- 15L
    p.geom <- 1.0 - 1.0/r

    offset <- 0.0
    if (! is.null(lm)) {
        lm.round <- round(lm)
        offset <- rep(0.0, length(lm))
        idx <- is.infinite(lm)
        offset[! idx] <- lm[! idx] - lm.round[! idx]
        lm <- lm.round
    }

    idx <- -0.5 == offset
    lm[idx] <- lm[idx] - 1L
    offset[idx] <- offset[idx] + 1.0

    if (1 == length(offset)) {
        offset <- rep(offset, length(p))
    }

    arg.data <- data.frame(
        p = p,
        offset = offset,
        p.geom = p.geom,
        r = r
    )

    data.f <- as.factor(paste0(offset, '/', p.geom))
    data.s <- split(arg.data, data.f)
    m <- lapply(data.s, function(data) {
        p <- data$p
        offset <- data$offset[1]
        p.geom <- data$p.geom[1]
        r <- data$r[1]
        m <- rep(NA, length(p))

        if(lower.tail) {
            m[0.0 == p] <- 0
            p.max <- vismeteor::pvmgeom(m.max + offset, r, perception.fun = perception.fun)
            idx <- p>p.max & p<=1.0
            if (any(idx)) {
                m[idx] <- m.max + 1 + stats::qgeom((p[idx] - p.max)/(1.0 - p.max), p.geom)
            }

            idx <- p>0.0 & p<=p.max
            if (any(idx)) {
                m0 <- seq(1, -m.max, -1)
                p0 <- vismeteor::pvmgeom(m0, r, lm = offset, perception.fun = perception.fun)
                p.idx <- findInterval(p[idx], p0, left.open = TRUE) + 1
                m[idx] <- -m0[p.idx]
                m[m<0] <- NA
            }
        } else {
            m[1.0 == p] <- 0
            p.max <- 1.0 - vismeteor::pvmgeom(m.max + offset, r, perception.fun = perception.fun)

            idx <- p>=0.0 & p<p.max
            if (any(idx)) {
                m[idx] <- m.max + 1 + stats::qgeom(p[idx]/p.max, p.geom, lower.tail = FALSE)
            }

            idx <- p>=p.max & p<1.0
            if (any(idx)) {
                m0 <- seq(-m.max, 0, 1)
                p0 <- c(vismeteor::pvmgeom(m0, r, lm = offset, lower.tail = FALSE, perception.fun = perception.fun), 1.0)
                m0 <- c(m0, 1)
                p.idx <- findInterval(p[idx], p0, left.open = FALSE)
                m[idx] <- -m0[p.idx]
            }
        }

        m
    })

    m <- unsplit(m, data.f)

    if (! is.null(lm)) {
        m <- lm - m
    }

    if (anyNA(m)) {
        warning('NaNs produced')
    }

    as.numeric(m)
}

#' @rdname vmgeom
#' @export
rvmgeom <- function(n, r, lm = NULL, perception.fun = NULL) {
    if (1.0 > r) {
        stop(paste0('r must be greater than 1.0 instead of "', r, '"!'))
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    p <- stats::runif(n)
    m <- rep(NA, n)
    idx <- p < 0.5

    if (any(idx)) {
        m[idx] <- vismeteor::qvmgeom(p[idx], r, lm = lm, lower.tail = TRUE, perception.fun = perception.fun)
    }

    if (any(!idx)) {
        m[!idx] <- vismeteor::qvmgeom(1.0 - p[!idx], r, lm = lm, lower.tail = FALSE, perception.fun = perception.fun)
    }

    m
}
