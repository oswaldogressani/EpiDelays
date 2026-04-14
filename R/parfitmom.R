#' Method of moments for single or doubly interval-censored data
#'
#' @description
#' This routine admits as data input a data frame \code{x} with either two
#' columns (named \code{xl} and \code{xr}) representing the left
#' and right bound, respectively, of the delay variable, or four columns
#' (named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}) representing the left
#' and right bound, respectively, of the primary and secondary events of the
#' delay variable. The naming convention of the columns is strict and different
#' namings are not allowed. When data frame \code{x} has two columns, the left
#' bound should be strictly smaller than the right bound, i.e. \code{xl < xr}
#' must be satisfied for all rows in \code{x}. When data frame \code{x} has four
#' columns, \code{x1l < x1r} and \code{x2l < x2r} must hold for all rows in
#' \code{x}. Moreover, \code{NA} values are not allowed in data frame \code{x}.
#'
#' @details
#' When a data frame with four columns is specified (i.e. doubly
#' interval-censored data), the routine reduces the data to a single interval
#' with endpoints \code{xl} and \code{xr}. The empirical observations for
#' moment matching are computed according to a midpoint imputation rule on the
#' intervals.
#'
#' @param x A data frame with either two columns named \code{xl} and \code{xr},
#' or four columns named \code{x1l} and \code{x1r}, \code{x2l}, \code{x2r}. See
#' description for constraints imposed on the columns.
#' @param family A character string specifying the name of the parametric
#' family. Can be one of the following: \code{"gaussian"}, \code{"gamma"},
#' \code{"lognormal"}, \code{"weibull"}, or \code{"skewnorm"}.
#' @param incheck Should internal checks be implemented to check that data
#' frame \code{x} is in correct format and that \code{family} has been
#' correctly specified by the user? Default is TRUE.
#' @param L Lower truncation point. Accepted for signature compatibility with
#' \code{parfitml()}; currently ignored. Under non-trivial truncation the
#' moment estimators below are biased (sample mean and variance drift toward
#' the centre of \code{[L, D]}), but the result is only used as an optimiser
#' starting value, so the bias is tolerated.
#' @param D Upper truncation point. See \code{L}.
#' @param dprimary Ignored. Accepted for signature consistency with
#' \code{\link{parfitml}} so the same argument set can be threaded through
#' both routines. Moment matching assumes uniform primary onset.
#' @param dprimary_args Ignored. Accepted for signature consistency with
#' \code{\link{parfitml}}.
#'
#' @return A list containing information on the chosen parametric family and
#' method of moments (MoM) estimates of the model parameters. The value
#' \code{mompoint_ub} gives the MoM estimates of the transformed parameters in
#' an unbounded parameter space (useful for starting point in maximum likelihood
#' estimation).
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export

parfitmom <- function(x, family, incheck = TRUE, L = 0, D = Inf,
                      dprimary = stats::dunif, dprimary_args = list()) {
  # L and D are accepted so parfitml() can forward them, but the moment
  # estimators below do not correct for truncation. The resulting
  # `mompoint_ub` is a known-biased optimiser start under non-trivial
  # truncation; parfitml()'s optim call still converges in practice.
  if (!is.numeric(L) || length(L) != 1L || is.na(L) || L < 0) {
    stop("L must be a non-negative scalar.")
  }
  if (!is.numeric(D) || length(D) != 1L || is.na(D) || L >= D) {
    stop("L must be less than D.")
  }
  if(!is.logical(incheck)) {
    stop("incheck must be either TRUE or FALSE")
  } else if (isTRUE(incheck)) {
    # Input checks
    dfck <- kerdata_check(x = x) # data frame check
    if (dfck$result == "fail") {
      stop(dfck$message)
    }
    famck <- kerfamily_check(x = family) # family check
    if (famck$result == "fail") {
      stop(famck$message)
    }
    domck <- kerdomain_check(x = x, family = family) # domain check
    if (domck$result == "fail") {
      stop(domck$message)
    }
  }
  fset <- kerfamilies()
  famdesc <- fset[[match(family, sapply(fset, "[[", "fname"))]]
  n <- nrow(x)
  nc <- ncol(x)
  if (nc == 2) {
    xl <- x$xl
    xr <- x$xr
  } else if (nc == 4) {
    xl <- x$x2l - x$x1r
    xr <- x$x2r - x$x1r
  }
  y <- 0.5 * (xl + xr)
  m1 <- mean(y)                   # First sample moment
  m2 <- (1 / n) * sum((y - m1)^2) # Second sample central moment
  m3 <- (1 / n) * sum((y - m1)^3) # Third sample central moment
  if(family == "gaussian") {
    par1approx <- m1
    par2approx <- sqrt(m2)
    mompoint_ub <- c(par1approx, log(par2approx))
    lpout <- list(par1approx = par1approx, par2approx = par2approx,
                  mompoint_ub = mompoint_ub)
  } else if (family == "skewnorm") {
    rsn <- 2 * m3 / (4 - pi)
    # Use a signed cube root so negative sample skew does not produce NaN
    # seeds. The original `rsn^(2/3)` returns NaN for rsn < 0 in R, which
    # propagates into v0 and crashes optim before the first likelihood
    # evaluation.
    cbrt <- function(z) sign(z) * abs(z)^(1 / 3)
    rsn_13 <- cbrt(rsn)
    rsn_23 <- rsn_13^2
    par2approx <- sqrt(m2 + rsn_23)
    dapprox <- rsn_13 / (par2approx * sqrt(2 / pi))
    # Clamp |dapprox| < 1 so the slant remains finite even in small-sample
    # cases where the moment identity pushes dapprox past the feasible
    # boundary. The MoM estimate is only an optimiser seed.
    dapprox_clamped <- max(min(dapprox, 0.99), -0.99)
    par3approx <- dapprox_clamped / sqrt(1 - dapprox_clamped^2)
    par1approx <- m1 - par2approx * dapprox_clamped * sqrt(2 / pi)
    if (!all(is.finite(c(par1approx, par2approx, par3approx))) ||
        par2approx <= 0) {
      par1approx <- m1
      par2approx <- sqrt(max(m2, 1e-6))
      par3approx <- 0
    }
    mompoint_ub <- c(par1approx, log(par2approx), par3approx)
    lpout <- list(par1approx = par1approx, par2approx = par2approx,
                  par3approx = par3approx, mompoint_ub = mompoint_ub)
  } else if (family == "gamma") {
    par1approx <- (m1^2 / m2)
    par2approx <- m1 / m2
    mompoint_ub <- log(c(par1approx, par2approx))
    lpout <- list(par1approx = par1approx, par2approx = par2approx,
                  mompoint_ub = mompoint_ub)
  } else if (family == "lognormal") {
    par1approx <- 2 * log(m1) - 0.5 * log(m1^2 + m2)
    par2approx <- sqrt(log(1 + m2 / (m1^2)))
    mompoint_ub <- c(par1approx, log(par2approx))
    lpout <- list(par1approx = par1approx, par2approx = par2approx,
                  mompoint_ub = mompoint_ub)
  } else if (family == "weibull") {
    mr <- m2 / m1^2
    f <- function(par1) {
      (gamma(1 + 2 / par1) / (gamma(1 + 1 / par1)^2)) - 1 - mr
    }
    lb <- 0
    flb <- f(lb)
    while (is.na(flb) || is.infinite(flb)) {
      lb <- lb + 0.01
      flb <- f(lb)
    }
    ub <- lb + 1
    par1approx <- stats::uniroot(f, lower = lb, upper = ub,
                                 extendInt = "downX")$root
    par2approx <- m1 / gamma(1 + 1 / par1approx)
    mompoint_ub <- log(c(par1approx, par2approx))
    lpout <- list(par1approx = par1approx, par2approx = par2approx,
                  mompoint_ub = mompoint_ub)
  }
  o <- c(famdesc, lpout)
  return(o)
}
