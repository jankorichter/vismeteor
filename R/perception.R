#' @title Perception Probabilities of Visual Meteor Magnitudes
#' @description
#' Provides the perception probability of visual meteor magnitudes and its first derivative.
#' @param m numerical; difference between the limiting magnitude and the meteor magnitude.
#' @param deriv boolean; when set to true, it returns the first derivative of the the perception probability.
#' @details
#' The perception probabilities of _Koschack R., Rendtel J., 1990b_
#' are estimated with the formula
#' \deqn{
#'     p(m) = \begin{cases}
#'         1.0 - \exp\left(-z(m + 0.5)\right)\  & \text{ if } m > -0.5,\\
#'         0.0 \  & \text{ otherwise,}
#'     \end{cases}
#' }
#' where
#' \deqn{
#' z(x) = 0.003 \, x + 0.0056 \, x^2 + 0.0014 \, x^4
#' }
#' and `m` is the difference between the limiting magnitude and the meteor magnitude.
#' @return This function returns the visual perception probabilities.
#' If `deriv` is set to `TRUE`, it will return the first derivative of these probabilities.
#' @references Koschack R., Rendtel J., 1990b _Determination of spatial number density and mass index from visual meteor observations (II)._ WGN 18, 119–140.
#' @examples
#' # Perception probability of visually estimated meteor of magnitude 6.0
#' # with a limiting magnitude of 5.6.
#' vmperception(5.6 - 6.0)
#' @export
vmperception <- function(m, deriv = FALSE) {
    poly.coef <- c(0.0, 0.003, 0.0056, 0, 0.0014)
    names(poly.coef) <- seq(along = poly.coef) - 1 # exponents

    m <- m + 0.5
    p <- rep(0.0, length(m))
    idx <- m > .Machine$double.eps
    if (any(idx)) {
        if (deriv) {
            inner0 <- f.polynomial(m[idx], poly.coef)
            poly.coef1 <- f.polynomial.coef(poly.coef, deriv.degree = 1L)
            inner1 <- f.polynomial(m[idx], poly.coef1)
            p[idx] <- exp(-inner0) * inner1
        } else {
            inner0 <- f.polynomial(m[idx], poly.coef)
            p[idx] <- 1.0 - exp(-inner0)
        }
    }

    p
}

#' @title Laplace-Transformed Perception Probabilities of Visual Meteor Magnitudes
#' @description
#' Provides the Laplace-transformed perception probability of visual meteor magnitudes
#' and its first derivative.
#' @param s numerical; Real (non-complex) parameter for the Laplace transformation.
#' @param deriv boolean; when set to true, it returns the first derivative of the transformation.
#' @details
#' The Laplace-transformed [perception probabilities][vismeteor::vmperception] `P(s)`, given as
#' \deqn{
#' P(s) = s \, \mathcal{B} \left\{p\right\}(s)
#' = \int_{-\infty}^{\infty} \, p(m) \, s \, \mathrm e^{-s \, m} \,\mathrm{d}m \,,
#' }
#' are approximately
#' \deqn{
#'     P(s) = \begin{cases}
#'         \exp\left(-4.11 \, s + 1.32 \, s^2 - 0.15 \, s^3\right)\ & \text{ if } s >= 0.0,\\
#'         \text{undefined} \  & \text{ otherwise.}
#'     \end{cases}
#' }
#' Here, `m` is the difference between the limiting magnitude and the meteor magnitude,
#' and `p(m)` denotes the perception probabilities as a function of `m`.
#' The \eqn{\mathcal{B}} recalls here the "two-sided Laplace transform" or
#' "bilateral Laplace transform". The term \eqn{s \, \mathcal{B}} is used to generate
#' the Laplace transform of the first derivative of the perception probabilities.
#' This ensures that \eqn{s \, \mathcal{B}} always lies between `0.0` and `1.0`.
#' On the other hand, \eqn{s^{-1} \, P(s)} yields the two-sided Laplace transform
#' of the perception probabilities.
#' @return This function returns the Laplace-transformed of the first derivative of
#' the perception probabilities. If `deriv` is set to `TRUE`, it will return the
#' first derivative of these Laplace-transformed values.
#' @seealso [vismeteor::vmperception]
#' @examples
#' vmperception.l(c(0, 0.5, Inf))
#' @export
vmperception.l <- function(s, deriv = FALSE) {
    poly.coef <- c(0.0, -4.11, 1.32, -0.15)
    names(poly.coef) <- seq(along = poly.coef) - 1 # exponents

    L <- rep(NA, length(s))
    idx <- s >= -.Machine$double.eps & s != Inf
    if (any(idx)) {
        if (deriv) {
            # exp(f(s)) * f'(s)
            f0 <- f.polynomial(s[idx], poly.coef)
            poly.coef1 <- f.polynomial.coef(poly.coef, deriv.degree = 1L)
            f1 <- f.polynomial(s[idx], poly.coef1)
            L[idx] <- exp(f0) * f1
        } else {
            # exp(f(s))
            f0 <- f.polynomial(s[idx], poly.coef)
            L[idx] <- exp(f0)
        }
    }
    L[Inf == s] <- 0.0

    L
}

f.polynomial <- function(m, poly.coef) {
    exponents <- as.numeric(names(poly.coef))
    margin.table(poly.coef * t(outer(m, exponents, "^")), 2)
}

f.polynomial.coef <- function(poly.coef, deriv.degree = 1L) {
    if (0L == deriv.degree)
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

    f.polynomial.coef(poly.coef, deriv.degree - 1L)
}
