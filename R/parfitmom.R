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
#' or four columns named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}. See
#' description for constraints imposed on the columns.
#' @param family A character string specifying the name of the parametric
#' family. Can be one of the following: \code{"gaussian"}, \code{"gamma"},
#' \code{"lognormal"}, \code{"weibull"}, or \code{"skewnorm"}.
#' @param incheck Should internal checks be implemented to check that data
#' frame \code{x} is in correct format and that \code{family} has been
#' correctly specified by the user? Default is TRUE.
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

parfitmom <- function(x, family, incheck = TRUE){
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
    par2approx <- sqrt(m2 + rsn^(2 / 3))
    dapprox <- (rsn^(1 / 3)) / (par2approx * sqrt(2 / pi))
    if(abs(dapprox) >= 1) {
      par1approx <- m1
      par2approx <- sqrt(m2)
      par3approx <- 0
    } else{
      par3approx <- dapprox / sqrt(1 - dapprox^2)
      par1approx <- m1 - par2approx * dapprox * sqrt(2 / pi)
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









