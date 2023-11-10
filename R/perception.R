#' @title Perception probability of visual meteor magnitudes
#' @description
#' Returns the perception probability of visual meteor magnitudes and its first derivation.
#' @param m numerical; difference between the limiting magnitude and the meteor magnitude.
#' @param deriv boolean; if true, the first derivation of the perception probability is returned.
#' @details
#' The perception probabilities according to _Koschack R., Rendtel J., 1990b_
#' are approximated with the formula
#' \deqn{
#'     p(m) = \begin{cases}
#'         1.0 - \exp\left(-z(m + 0.5)\right)\  & \text{ if } m > -0.5,\\
#'         0.0 \  & \text{ otherwise,}
#'     \end{cases}
#' }
#' with
#' \deqn{
#' z(x) = 0.00028 \, x + 0.0081 \, x^2 + 0.0012 \, x^4
#' }
#' where `m` is the difference between the limiting magnitude and the meteor magnitude.
#' @return This function returns the visual perception probabilities.
#' If `deriv` is `TRUE`, it will return the first derivative
#' of the visual perception probabilities.
#' @references Koschack R., Rendtel J., 1990b _Determination of spatial number density and mass index from visual meteor observations (II)._ WGN 18, 119–140.
#' @examples
#' # Perception probability of visually estimated meteor of magnitude 6.0
#' # with a limiting magnitude of 5.6.
#' vmperception(5.6 - 6.0)
#' @export
vmperception <- function(m, deriv = FALSE) {
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

    m <- m + 0.5
    poly.coef <- c(0.00028, 0.0081, 0, 0.0012)
    names(poly.coef) <- seq(along = poly.coef) # exponents

    p <- rep(0.0, length(m))
    if (deriv) {
        idx <- m > .Machine$double.eps
        if (any(idx)) {
            inner0 <- f.polynomial(m[idx], poly.coef)
            inner1 <- f.polynomial(m[idx], deriv.polynomial(poly.coef, degree = 1L))
            p[idx] <- exp(-inner0) * inner1
        }
    } else {
        idx <- m > .Machine$double.eps
        if (any(idx)) {
            inner0 <- f.polynomial(m[idx], poly.coef)
            p[idx] <- 1.0 - exp(-inner0)
        }
    }

    p
}
