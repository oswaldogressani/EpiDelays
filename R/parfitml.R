#' Parametric estimation with maximum likelihood for interval-censored data
#'
#' @description
#' This routine fits commonly used parametric models to interval-censored data.
#'
#' @usage parfitml(x, family = "gamma", boot = FALSE, Bboot = 200, pgbar = FALSE)
#'
#' @param x A data frame with n rows and two columns indicating the left and
#'  right bound, respectively, of the interval-censored observations. \code{NA}
#'  values are not allowed.
#'
#' @param family The name of the parametric family used in modeling the data.
#'  Must be one of: "gamma", "lognormal", "weibull" or "gaussian".
#'
#' @param boot Should the bootstrap be implemented to obtain standard errors
#'  and 90\% and 95\% confidence intervals for different features of the
#'  distribution? Default is \code{FALSE}.
#'
#' @param Bboot The number of bootstrap replications. Default is 200.
#'
#' @param pgbar Should a progress bar be displayed while the bootstrap is
#'  ongoing? Default is \code{FALSE}.
#'
#' @details The \code{parfitml} routine provides point estimates and confidence
#' intervals (if bootstrap is called) for different features of the chosen
#' parametric distribution. The following features are considered: mean,
#' variance, standard deviation, and 1\%, 5\%, 25\%, 50\%, 75\%, 95\% and 99\%
#' percentiles. Maximum likelihood estimates of model parameters are computed
#' with the \code{optim} function with the Nelder and Mead (1965) method.
#' Initial parameter values are obtained via a moment matching approach, where
#' the theoretical mean and variance of the chosen parametric model are matched
#' with the empirical mean and variance of the midpoint imputed values computed
#' from the interval-censored observations.
#'
#' @return A list containing parameter estimates, feature estimates, AIC, BIC,
#' and convergence checks.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} and
#'  Dongxuan Chen.
#'
#' @references Nelder, J. A. and Mead, R. (1965). A simplex algorithm for
#' function minimization. \emph{Computer Journal}, \strong{7}, 308–313.
#'
#' @examples
#' # Add examples here
#'
#' @export

parfitml <- function(x, family = "gamma", boot = FALSE, Bboot = 200, pgbar = FALSE) {

  #--- Extract info from data
  n <- nrow(x)                                  # sample size
  model <- logliksingle(x = x, family = family) # parametric model features

  #--- Maximum likelihood estimation
  maxlik <- function(x){
    uparinit <- midmom(x = x, family = family)  # Initial values via moments
    ml <- stats::optim(par = uparinit, fn = model$ll, x = x,
                       control = list(fnscale = -1))
    mlconv <- ml$convergence
    mlpar <- model$originsc(ml$par)
    epifeat <- model$features(mlpar)
    npar <- length(mlpar)
    llmax <- model$ll(ml$par, x)
    aic <- 2 * (npar - llmax)
    bic <- npar * log(n) - 2 * llmax
    outlist <- list(family = family, parspace = model$parspace, mlpar = mlpar,
                    mlconv = mlconv, aic = aic, bic = bic, epifeat = epifeat)
    return(outlist)
  }
  mle <- maxlik(x)

  #--- Nonparametric bootstrap
  if(isTRUE(boot)) {
    sstatsboot <- model$sstats
    convboot <- 0
    rowid <- seq_len(n)
    if(isTRUE(pgbar)) {
    cat(paste0("Nonparametric bootstrap (MLE-", family, " family) \n",
               "Bootstrap sample: B=", Bboot, " \n",
               "Sample size: n=", n," \n"))
    progbar <- utils::txtProgressBar(min = 1, max = Bboot, initial = 1,
                                     style = 3, char ="*")
    }
    for(b in 1:Bboot){
      bootrow <- sample(rowid, size = n, replace = TRUE)
      xboot <- x[bootrow, ]
      mleboot <- maxlik(xboot)
      if(mleboot$mlconv == 0) convboot <- convboot + 1
      sstatsboot[b, ] <- mleboot$epifeat
      if(isTRUE(pgbar)) utils::setTxtProgressBar(progbar, b)
    }
    if(isTRUE(pgbar)) close(progbar)
    se <- apply(sstatsboot, 2, stats::sd)
    CI90L <- apply(sstatsboot, 2, stats::quantile, prob = 0.05)
    CI90R <- apply(sstatsboot, 2, stats::quantile, prob = 0.95)
    CI95L <- apply(sstatsboot, 2, stats::quantile, prob = 0.025)
    CI95R <- apply(sstatsboot, 2, stats::quantile, prob = 0.975)
    epifeatboot <- t(rbind(mle$epifeat, se, CI90L, CI90R, CI95L, CI95R))
    colnames(epifeatboot) <- c("estim", "se", "CI90L", "CI90R",
                               "CI95L", "CI95R")
    epifeatboot <- as.data.frame(epifeatboot)
  } else{
    epifeatboot <- NULL
    Bboot <- NULL
    convboot <- NULL
  }

  #--- Routine output
  outlist <- c(mle, list(Bboot = Bboot, convboot = convboot,
                         epifeatboot = epifeatboot))
  return(outlist)
}
