#' @title Perception probability of visual meteor magnitudes
#' @description
#' Returns the perception probability of visual meteor magnitudes and its first derivation.
#' @param m numerical; difference between the limiting magnitude and the meteor magnitude.
#' @param log logical; if `TRUE`, the logarithmic value is returned.
#' @param deriv logical; if `TRUE`, its first derivation is returned.
#' @details
#' The perception probabilities according to _Koschack R., Rendtel J., 1990b_
#' are approximated with the formula
#' \deqn{
#' p(m) = (1 + \exp\left(-z(m)\right))^{-1}
#' }
#' with
#' \deqn{
#' z(m) = -6.01 + 3.54 \, m - 1.655 \, m^2 + 0.5468 \, m^3 - 0.0832 \, m^4 + 0.0046 \, m^5
#' }
#' where `m` is the difference between the limiting magnitude and the meteor magnitude.
#' @return This function returns the visual perception probabilities.
#' If `deriv = TRUE`, its first derivation is returned.
#' If `log = TRUE`, the logarithmic value of the perception probabilities is returned.
#' If both, `deriv = TRUE` and `log = TRUE`, the first derivation of its logarithmic value is returned.
#' @references Koschack R., Rendtel J., 1990b _Determination of spatial number density and mass index from visual meteor observations (II)._ WGN 18, 119–140.
#' @examples
#' # Perception probability of visually estimated meteor of magnitude 6.0
#' # with a limiting magnitude of 5.6.
#' vmperception(5.6 - 6.0)
#' @export
vmperception <- function(m, log = FALSE, deriv = FALSE) {
    p.inner <- -6.1 + 4.05*m - 2.15*m^2 + 0.72*m^3 - 0.1087*m^4 + 0.00595*m^5

    if (deriv) {
        p.inner.derived <- 4.05 - 4.3*m + 2.16*m^2 - 0.4348*m^3 + 0.02975*m^4

        if (log) {
            p.inner.derived/(1 + exp(p.inner))
        } else {
            p.inner.derived * exp(-p.inner)/(exp(-p.inner) + 1)^2
        }
    } else {
        p <- 1/(1 + exp(-p.inner))
        if(log) {
            base::log(p)
        } else {
            p
        }
    }
}
