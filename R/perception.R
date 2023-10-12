#' @title Perception probability of visual meteor magnitudes
#' @description
#' Returns the perception probability of visual meteor magnitudes and its derivations.
#' @param m numerical; difference between the limiting magnitude and the meteor magnitude.
#' @param deriv.degree integer; degree of derivation (0, 1, 2). Default is `0`.
#' @details
#' The perception probabilities according to _Koschack R., Rendtel J., 1990b_
#' are approximated with the formula
#' \deqn{
#' p(m) = (1 + \exp\left(-z(m)\right))^{-1}
#' }
#' with
#' \deqn{
#' z(m) = -6.01 + 4.05 \, m - 2.15 \, m^2 + 0.72 \, m^3 - 0.1087 \, m^4 + 0.00595 \, m^5
#' }
#' where `m` is the difference between the limiting magnitude and the meteor magnitude.
#' @return This function returns the visual perception probabilities.
#' If `deriv.degree` is set, it will return the corresponding derivative
#' of the visual perception probabilities. Currently, only the first and second derivatives
#' are supported.
#' @references Koschack R., Rendtel J., 1990b _Determination of spatial number density and mass index from visual meteor observations (II)._ WGN 18, 119–140.
#' @examples
#' # Perception probability of visually estimated meteor of magnitude 6.0
#' # with a limiting magnitude of 5.6.
#' vmperception(5.6 - 6.0)
#' @export
vmperception <- function(m, deriv.degree = 0) {
    deriv.polynomial <- function(poly.coef, degree) {
        if (0 == degree)
            return(poly.coef)

        if (1 == length(poly.coef))
            return(0)

        poly.coef <- poly.coef[-1]
        deriv.polynomial(poly.coef * seq(along = poly.coef), degree - 1)
    }

    f.inner <- function(m, poly.coef) {
        margin.table(poly.coef * t(outer(m, seq(along=poly.coef) - 1, "^")), 2)
    }

    poly.coef <- c(-6.1, 4.05, -2.15, 0.72, -0.1087, 0.00595)

    if (0 == deriv.degree) {
        1/(1 + exp(-f.inner(m, poly.coef)))
    } else if (1 == deriv.degree) {
        inner0 <- f.inner(m, poly.coef)
        inner1 <- f.inner(m, deriv.polynomial(poly.coef, 1))
        exp.inner0 <- exp(-inner0)
        inner1 * exp.inner0/(exp.inner0 + 1)^2
    } else if (2 == deriv.degree) {
        inner0 <- f.inner(m, poly.coef)
        inner1 <- f.inner(m, deriv.polynomial(poly.coef, 1))
        inner2 <- f.inner(m, deriv.polynomial(poly.coef, 2))
        exp.inner0 <- exp(-inner0)
        (exp.inner0 * ((exp.inner0 - 1) * inner1^2 + (exp.inner0 + 1) * inner2))/(exp.inner0 + 1)^3
    } else {
        stop(paste("Invalid degree of derivation:", deriv.degree))
    }
}
