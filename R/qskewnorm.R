#' Quantile function of the Skew-Normal Distribution
#'
#' @param p A vector of probabilities.
#' @param par1 Location parameter of the skew-normal distribution.
#' @param par2 Scale parameter of the skew-normal distribution.
#' @param par3 Slant parameter of the skew-normal distribution.
#'
#' @details The \code{qskewnorm} routine use Owen's (1956)
#' \eqn{T(h,a)} function and a \code{uniroot} procedure for quantile
#' computation.
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
#' qskewnorm(p = c(0.1,0.9))
#'
#' @export

qskewnorm <- function(p, par1 = 0, par2 = 1, par3 = 0) {
  if (any(p <= 0) | any(p >= 1))
    stop("Values in p must be in (0,1)")
  d3 <- par3 / sqrt(1 + par3^2)
  mean <- par1 + par2 * d3 * sqrt(2 / pi)
  sd <- sqrt(par2^2 * (1 - (2 / pi) * d3^2))
  xpl <- mean - sd * sqrt((1 - p) / p)
  xpu <- mean + sd * sqrt(p / (1 - p))
  g <- function(q, prob)
    pskewnorm(q, par1 = par1, par2 = par2, par3 = par3) - prob
  xma <- cbind(p, xpl, xpu)
  frow <- function(r)
    stats::uniroot(f = g, lower = r[2], upper = r[3], p = r[1],
                   extendInt = "upX")$root
  o <- apply(xma, 1, frow)
  return(o)
}
