#' Compute initial values for maximum likelihood estimation via moment matching
#'
#' @keywords internal
#'
#' @param x A data frame with n rows and two columns indicating the left and
#'  right bound, respectively, of the interval-censored observations. \code{NA}
#'  values are not allowed.
#'
#' @param family The name of the parametric family used in modeling the data.
#'  Must be one of: "gamma", "lognormal", "weibull" or "gaussian".
#'
#' @return A vector of parameter values.
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export

midmom <- function(x, family = "gamma") {
  tl <- x[, 1]           # left bound of censored observation
  tr <- x[, 2]           # right bound of censored observation
  tm <- (tl + tr) * 0.5  # midpoints of censoring interval
  Etm <- mean(tm)        # mean of midpoints
  Vtm <- stats::var(tm)         # variance of midpoints
  if (family == "gamma") {
    uparinit <- log(c(Etm^2 / Vtm, Etm / Vtm))
  } else if (family == "lognormal") {
    locmom <-  2 * log(Etm) - 0.5 * log(Etm^2 + Vtm)
    uparinit <- c(locmom, log(sqrt(2 * (log(Etm) - locmom))))
  } else if (family == "weibull") {
    fw <- function(shape) {
      val <- Etm^2 * (gamma(1 + 2 / shape) / (gamma(1 + 1 / shape)^2) - 1) - Vtm
      return(val)
    }
    lb <- 1e-5
    fwlb <- fw(lb)
    while (is.na(fwlb) || is.infinite(fwlb)) {
      lb <- lb + 0.1
      fwlb <- fw(lb)
    }
    ub <- lb + 0.1
    fwub <- fw(ub)
    while (fwub >= 0) {
      ub <- ub + 0.1
      fwub <- fw(ub)
    }
    shapemom <- stats::uniroot(fw, lower = lb, upper = ub)$root
    uparinit <- log(c(shapemom, Etm / gamma(1 + 1 / shapemom)))
  } else if (family == "gaussian") {
    uparinit <- c(Etm, log(sqrt(Vtm)))
  }
  return(uparinit)
}
