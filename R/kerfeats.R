#' Computes features for a given parametric family
#'
#' @description
#' Computes the mean, variance, standard deviation and different quantiles
#' given a parametric family and parameter values.
#'
#' @param par A vector of parameter values.
#' @param family A characther string specifying the parametric family.
#' @param p A vector of probabilities.
#'
#' @returns
#' A numeric vector of features.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerfeats <- function(par, family, p) {
  if (family == "gaussian") {
    mean <- par[1]
    var <- par[2]^2
    sd <- par[2]
    qp <- stats::qnorm(p = p, mean = par[1], sd = par[2])
  } else if (family == "skewnorm") {
    d3 <- par[3] / sqrt(1 + par[3]^2)
    mean <- par[1] + par[2] * d3 * sqrt(2 / pi)
    var <- par[2]^2 * (1 - (2 / pi) * d3^2)
    sd <- sqrt(var)
    qp <- qskewnorm(p = p, par1 = par[1], par2 = par[2], par3 = par[3])
  } else if (family == "gamma") {
    mean <- par[1] / par[2]
    var <- par[1] / (par[2]^2)
    sd <- sqrt(var)
    qp <- stats::qgamma(p = p, shape = par[1], rate = par[2])
  } else if (family == "lognormal") {
    mean <- exp(par[1] + 0.5 * par[2]^2)
    var <- exp(2 * par[1] + par[2]^2) * (exp(par[2]^2) - 1)
    sd <- sqrt(var)
    qp <- stats::qlnorm(p = p, meanlog = par[1], sdlog = par[2])
  } else if (family == "weibull") {
    mean <- par[2] * gamma(1 + 1 / par[1])
    var <- par[2]^2 * (gamma(1 + 2 / par[1]) - (gamma(1 + 1 / par[1])^2))
    sd <- sqrt(var)
    qp <- stats::qweibull(p = p, shape = par[1], scale = par[2])
  }
  o <- c(mean, var, sd, qp)
  return(o)
}






