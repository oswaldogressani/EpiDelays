#' Distribution function of the Skew-Normal Distribution
#'
#' @param x A vector of quantiles.
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
#' pskewnorm(x = c(0.1,1.3))
#'
#' @export

pskewnorm <- function(x, par1 = 0, par2 = 1, par3 = 0) {
  z <- (x - par1) / par2
  o <- stats::pnorm(z) - 2 * sapply(z, kerTOw, a = par3)
  # Owen's T can return values a hair outside [0, 1] in extreme parameter
  # regions (large |slant| or deep tails), and floating-point noise can leave
  # the returned vector slightly non-monotone in x. A CDF must satisfy both,
  # so we clamp into [0, 1] and lift any rounding-noise dips with cummax over
  # the order of x. primarycensored::check_pdist enforces both properties on
  # any pdist passed to dprimarycensored() and would otherwise abort during
  # an optim step that visits the saturation region.
  o <- pmin(1, pmax(0, o))
  ord <- order(x)
  o[ord] <- cummax(o[ord])
  o
}
