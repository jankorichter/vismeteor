#' @title Perception probability of visual meteor magnitudes
#' @description
#' Returns the perception probability of visual meteor magnitudes and its first derivation.
#' @param m numerical; difference between the limiting magnitude and the meteor magnitude.
#' @param deriv boolean; if true, the first derivation of the perception probability is returned.
#' @details
#' The perception probabilities according to _Koschack R., Rendtel J., 1990b_
#' are approximated with the formula
#' \deqn{
#' p(m) = (1 + \exp\left(-z(m)\right))^{-1}
#' }
#' with
#' \deqn{
#' z(m) = -6.24 + 3.33 \, m - 1.03 \, m^2 + 0.241 \, m^3 - 0.03 \, m^4 + 0.0015 \, m^5
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

    poly.coef <- c(-6.24, 3.33, -1.03, 0.241, -0.03, 0.0015)

    if (deriv) {
        inner0 <- f.inner(m, poly.coef)
        inner1 <- f.inner(m, deriv.polynomial(poly.coef, 1))
        exp.inner0 <- exp(-inner0)
        inner1 * exp.inner0/(exp.inner0 + 1)^2
    } else {
        1/(1 + exp(-f.inner(m, poly.coef)))
    }
}
