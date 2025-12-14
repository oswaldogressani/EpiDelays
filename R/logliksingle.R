#' Log likelihood functions and model features for single interval-censored data
#'
#' @description
#' Provides the log likelihood function, features, and parameter space for a
#' given parametric family.
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
#' @return A list containing the log likelihood function and other model
#'  features.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} and
#'  Dongxuan Chen.
#'
#' @export

logliksingle <- function(x, family) {
  qfeat <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)
  sstats <- data.frame(mean = numeric(0), var = numeric(0), sd = numeric(0),
                       q1 = numeric(0), q5 = numeric(0), q25 = numeric(0),
                       median = numeric(0), q75 = numeric(0), q95 = numeric(0),
                       q99 = numeric(0))
  if (family == "gamma") {
    ll <- function(upar, x) {
      # Loglik upar: unbounded parameter; x: data
      tl <- x[, 1]
      tr <- x[, 2]
      par  <- exp(upar)
      Ftr  <- stats::pgamma(q = tr, shape = par[1], rate = par[2])
      Ftl  <- stats::pgamma(q = tl, shape = par[1], rate = par[2])
      llv  <- sum(log(Ftr - Ftl))
      return(llv)
    }
    originsc <- function(upar) exp(upar)
    features <- function(par) {
      meang <- par[1] / par[2]
      varg <- par[1] / (par[2]^2)
      sdg <- sqrt(varg)
      qg <- stats::qgamma(p = qfeat, shape = par[1], rate = par[2])
      sstats[1, ] <- c(meang, varg, sdg, qg)
      return(sstats)
    }
    parspace <- c("par[1]: shape > 0", "par[2]: rate > 0")
  } else if (family == "lognormal") {
    ll <- function(upar, x) {
      # Loglik upar: unbounded parameter; x: data
      tl <- x[, 1]
      tr <- x[, 2]
      par  <- c(upar[1], exp(upar[2]))
      Ftr  <- stats::plnorm(q = tr, meanlog = par[1], sdlog = par[2])
      Ftl  <- stats::plnorm(q = tl, meanlog = par[1], sdlog = par[2])
      llv  <- sum(log(Ftr - Ftl))
      return(llv)
    }
    originsc <- function(upar) c(upar[1], exp(upar[2]))
    features <- function(par) {
      meanln <- exp(par[1] + 0.5 * par[2]^2)
      varln <- exp(2 * par[1] + par[2]^2) * (exp(par[2]^2) - 1)
      sdln <- sqrt(varln)
      qln <- stats::qlnorm(p = qfeat, meanlog = par[1], sdlog = par[2])
      sstats[1, ] <- c(meanln, varln, sdln, qln)
      return(sstats)
    }
    parspace <- c("par[1]: location (real)", "par[2]: scale > 0")
  } else if (family == "weibull") {
    ll <- function(upar, x) {
      # Loglik upar: unbounded parameter; x: data
      tl <- x[, 1]
      tr <- x[, 2]
      par  <- exp(upar)
      Ftr  <- stats::pweibull(q = tr, shape = par[1], scale = par[2])
      Ftl  <- stats::pweibull(q = tl, shape = par[1], scale = par[2])
      llv  <- sum(log(Ftr - Ftl))
      return(llv)
    }
    originsc <- function(upar) exp(upar)
    features <- function(par) {
      meanw <- par[2] * gamma(1 + 1 / par[1])
      varw <- par[2]^2 * (gamma(1 + 2 / par[1]) - (gamma(1 + 1 / par[1])^2))
      sdw <- sqrt(varw)
      qw <- stats::qweibull(p = qfeat, shape = par[1], scale = par[2])
      sstats[1, ] <- c(meanw, varw, sdw, qw)
      return(sstats)
    }
    parspace <- c("par[1]: shape > 0", "par[2]: scale > 0")
  } else if (family == "gaussian") {
    ll <- function(upar, x) {
      # Loglik upar: unbounded parameter; x: data
      tl <- x[, 1]
      tr <- x[, 2]
      par  <- c(upar[1], exp(upar[2]))
      Ftr  <- stats::pnorm(q = tr, mean = par[1], sd = par[2])
      Ftl  <- stats::pnorm(q = tl, mean = par[1], sd = par[2])
      llv  <- sum(log(Ftr - Ftl))
      return(llv)
    }
    originsc <- function(upar) c(upar[1], exp(upar[2]))
    features <- function(par) {
      meangauss <- par[1]
      vargauss <- par[2]^2
      sdgauss <- par[2]
      qgauss <- stats::qnorm(p = qfeat, mean = par[1], sd = par[2])
      sstats[1, ] <- c(meangauss, vargauss, sdgauss, qgauss)
      return(sstats)
    }
    parspace <- c("par[1]: location (real)", "par[2]: scale > 0")
  }

  outlist <- list(ll = ll, originsc = originsc, features = features,
                  parspace = parspace, sstats = sstats)

  return(outlist)
}









