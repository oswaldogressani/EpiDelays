#' Distribution function of the Skew-Normal Distribution
#'
#' @param q A vector of quantiles.
#' @param par1 Location parameter of the skew-normal distribution.
#' @param par2 Scale parameter of the skew-normal distribution.
#' @param par3 Slant parameter of the skew-normal distribution.
#'
#' @details The \code{pskewnorm} routine uses Owen's (1956)
#' \eqn{T(h,a)} function to evaluate the distribution function.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references
#' Owen, D. B. (1956). Tables for Computing Bivariate Normal Probabilities.
#' \emph{The Annals of Mathematical Statistics}, \strong{27}(4), 1075–1090.
#'
#' Azzalini, A., and Capitanio, A. (2014). The Skew-Normal and Related Families.
#' \emph{Cambridge University Press}.
#'
#' @examples
#' pskewnorm(q = c(0.1,1.3))
#'
#' @export

pskewnorm <- function(q, par1 = 0, par2 = 1, par3 = 0) {
  z <- (q - par1) / par2
  o <- stats::pnorm(z) - 2 * sapply(z, kerTOw, a = par3)
  return(o)
}
