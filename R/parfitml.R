#' Parametric estimation with maximum likelihood for single or doubly
#' interval-censored data
#'
#' @description
#' Fit parametric models to single or doubly
#' interval-censored data with maximum likelihood. Data input must be a
#' data frame \code{x} with either two columns (named \code{xl} and \code{xr})
#' representing the left and right bound, respectively, of the delay variable,
#' or four columns (named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r})
#' representing the left and right bound, respectively, of the primary and
#' secondary events of the delay variable. The naming convention of the data
#' columns is strict and different namings are not allowed. When data frame
#' \code{x} has two columns, a single interval-censored likelihood function is
#' used. Note that the left bound should be strictly smaller than the right
#' bound, i.e. \code{xl < xr} must be satisfied for all rows in \code{x}.
#' When data frame \code{x} has four columns, a doubly interval-censored
#' likelihood function is used. In that case \code{x1l < x1r} and
#' \code{x2l < x2r} must hold for all rows in \code{x}. \code{NA} values are
#' not allowed in data frame \code{x}.
#'
#' @details This routine computes maximum likelihood estimates of model
#' parameters as well as different features of the fitted delay variable
#' distribution (mean, variance, standard deviation and selected quantiles). The
#' nonparametric bootstrap is used to compute standard errors and confidence
#' intervals for the model parameters and different features. Maximum likelihood
#' estimates are computed with the \code{optim} function using the
#' Nelder and Mead (1965) method. Initial parameter values are computed via a
#' moment matching approach.
#'
#' @param x A data frame with either two columns named \code{xl} and \code{xr},
#' or four columns named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}. See
#' description for constraints imposed on the columns.
#' @param family A character string specifying the name of the parametric
#' family. Can be one of the following: \code{"gaussian"}, \code{"gamma"},
#' \code{"lognormal"}, \code{"weibull"}, or \code{"skewnorm"}.
#' @param ci The engine for computing confidence intervals. Default is
#' \code{"npboot"} for nonparametric bootstrap. Other options are \code{"pboot"}
#' for the parametric bootstrap (currently under construction) or \code{"sbnorm"}
#' which relies on a simulation-based approach to compute confidence intervals
#' using asymptotic normality of the MLE estimator following Mandel (2013).
#' @param ... Further arguments. Specifying e.g. \code{Bboot = 1000} permits
#' to fix the bootstrap sample size to 1000 (the default is 100 here). Likewise
#' specifying \code{ns = 1000} permits to fix the number of simulated samples
#' from the Gaussian distribution of the MLE in the Mandel (2013) method (the
#' default is 100 here). Also, specifying \code{pgbar = TRUE} shows a progress
#' bar for the bootstrap.
#'
#' @return A list containing detailed information on the fitted parametric
#' model. The \code{summary()} function can be used to see further details.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} (original
#' code writing and implementation) and
#' Dongxuan Chen (code editing and testing).
#'
#' @references Nelder, J. A. and Mead, R. (1965). A simplex algorithm for
#' function minimization. \emph{Computer Journal}, \strong{7}, 308–313.
#'
#' @references Mandel, M. (2013). Simulation-based confidence intervals
#' for functions with complicated derivatives. \emph{The American Statistician},
#' \strong{67}(2), 76-81.
#'
#' @export

parfitml <- function(x, family, ci = c("npboot", "pboot", "sbnorm"),...){
  tic <- proc.time()
  m <- kerlikelihood(x = x, family = family)
  n <- nrow(x)
  np <- m$npars
  v0 <- parfitmom(x = x, family = family, incheck = FALSE)$mompoint_ub
  maxs <- list(fnscale = -1)
  mle <- stats::optim(par = v0, fn = m$loglik, x = x, control = maxs,
                      hessian = TRUE)
  mlepar <- as.numeric(m$originscale(mle$par))
  mlefeat <- as.numeric(kerfeats(family = family, par = mlepar))
  mleconv <- (mle$convergence == 0)
  cimethod <- match.arg(ci)
  parl <- stats::setNames(vector("list", np), paste0("par", 1:np))
  feats <- c("mean", "var", "sd", paste0("q", c(1, 5, 25, 50, 75, 95, 99)))
  featsl <- stats::setNames(vector("list", length(feats)), feats)
  if(cimethod == "npboot"){
    if ("Bboot" %in% ...names()) {
      Bboot <- list(...)$Bboot
      if (!(is.numeric(Bboot) && Bboot > 0))
        stop("Bboot must be a positive integer.")
      Bboot <- round(Bboot)
    } else{
      Bboot <- 100
    }
    if ("pgbar" %in% ...names()) {
      if (isTRUE(list(...)$pgbar)) {
        pgbar <- TRUE
        cat(paste0("Fitting parametric model (", family, ") \n",
                   "Bootstrap progress (Bboot=", Bboot, "): \n"))
        progbar <- utils::txtProgressBar(min = 1, max = Bboot, initial = 1,
                                         style = 3, char ="*")
      }
    } else{
      pgbar <- FALSE
    }
    pboot <- matrix(0, nrow = Bboot, ncol = np)
    fboot <- matrix(0, nrow = Bboot, ncol = length(mlefeat))
    bootdiscard <- 0
    for(b in 1:Bboot) {
      bootbconv <- 1
      while (bootbconv != 0) {
        xb <- kerboot(x)
        mleboot <- stats::optim(par = mle$par, fn = m$loglik, x = xb,
                                control = maxs)
        if (mleboot$convergence == 0) {
          bootbconv <- 0
        } else {
          bootdiscard <- bootdiscard + 1
        }
      }
      mleparboot <- as.numeric(m$originscale(mleboot$par))
      pboot[b,] <- mleparboot
      fboot[b, ] <- as.numeric(kerfeats(family = family, par = mleparboot))
      if(isTRUE(pgbar)){
        utils::setTxtProgressBar(progbar, b)
      }
    }
    if(isTRUE(pgbar)){
      close(progbar)
    }
    parfit <- kerstats(slist = parl, pestim = mlepar, method = "boot",
                       boot = pboot)
    delayfit <- kerstats(slist = featsl, pestim = mlefeat, method = "boot",
                         boot = fboot)
    ns <- NULL
  } else if(cimethod == "sbnorm"){
    if ("ns" %in% ...names()) {
      ns <- list(...)$ns
      if (!(is.numeric(ns) && ns > 0))
        stop("Number of samples for ci must be a positive integer.")
      ns <- round(ns)
    } else{
      ns <- 100
    }
    obsFishermle <- (-mle$hessian)
    invobsFishermle <- chol2inv(chol(obsFishermle))
    Jmle <- m$J(mle$par)
    sigmamle <- crossprod(Jmle, invobsFishermle %*% Jmle)
    semle <- sqrt(diag(sigmamle))
    psim <- rmvnorm(n = ns, mean = mlepar, sigma = sigmamle)$sim
    err2feats <- list()
    for(j in 1:ns) {
      err2feats[[j]] <- (kerfeats(family = family, par = psim[j, ]) - mlefeat)^2
      cat(sprintf("\r Simulation-based ci: %d/%d.", j, ns))
    }
    sefeats <- sqrt(colMeans(do.call(rbind, err2feats)))
    parfit <- kerstats(parl, mlepar, method = "norm", se = semle)
    delayfit <- kerstats(featsl, mlefeat, method = "norm", se = sefeats)
    Bboot <- NULL
    bootdiscard <- NULL
  }
  aic <- 2 * (np - mle$value)
  bic <- np * log(n) - 2 * mle$value
  toc <- proc.time() - tic

  o <- c(m[!names(m) %in% c("loglik", "originscale")],
         list(n = n, Bboot = Bboot, parfit = parfit, delayfit = delayfit,
              aic = aic, bic = bic, mleconv = mleconv,
              bootdiscard = bootdiscard, cimethod = cimethod, ns = ns,
              xmin = m$xmin, xmax = m$xmax,
              elapsed = toc[3]))
  attr(o, "class") <- "parfitml"
  return(o)
}



