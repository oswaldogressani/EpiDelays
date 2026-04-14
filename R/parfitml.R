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
#' intervals for the model parameters and different features. By default, the
#' number of bootstrap samples is fixed at 1000. Maximum likelihood estimates
#' are computed with the \code{optim} function using the Nelder and Mead (1965)
#' method. Initial parameter values are computed via a moment matching approach.
#'
#' @param x A data frame with either two columns named \code{xl} and \code{xr},
#' or four columns named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}. See
#' description for constraints imposed on the columns.
#' @param family A character string specifying the name of the parametric
#' family. Can be one of the following: \code{"gaussian"}, \code{"gamma"},
#' \code{"lognormal"}, \code{"weibull"}, or \code{"skewnorm"}.
#' @param Bboot Number of bootstrap samples. Default is 1000.
#' @param pgbar Should a progress bar be displayed in console? Default is TRUE.
#' @param L Lower truncation point of the underlying delay distribution.
#' Default is \code{0}. Must satisfy \code{L >= 0} and \code{L < D}. Data are
#' assumed to be drawn conditional on the underlying delay lying inside
#' \code{[L, D]}; the fitted parameters describe the unconditional delay
#' distribution.
#' @param D Upper truncation point of the underlying delay distribution.
#' Default is \code{Inf}.
#' @param dprimary Primary event density function. See
#' \code{\link{kerlikelihood}} for details. Defaults to \code{stats::dunif}
#' (uniform primary onset), which reproduces the behaviour of earlier
#' EpiDelays versions. Non-uniform choices include
#' \code{primarycensored::dexpgrowth}.
#' @param dprimary_args A named list of additional arguments passed to
#' \code{dprimary}. Defaults to an empty list.
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
#' @export

parfitml <- function(x, family, Bboot = 1000, pgbar = TRUE,
                     L = 0, D = Inf,
                     dprimary = stats::dunif, dprimary_args = list()) {
  # Validate truncation bounds up front, mirroring
  # primarycensored::.check_truncation_bounds.
  if (!is.numeric(L) || length(L) != 1L || is.na(L) || L < 0) {
    stop("L must be a non-negative scalar.")
  }
  if (!is.numeric(D) || length(D) != 1L || is.na(D) || L >= D) {
    stop("L must be less than D.")
  }
  tic <- proc.time()
  m <- kerlikelihood(x = x, family = family, L = L, D = D,
                     dprimary = dprimary,
                     dprimary_args = dprimary_args) # Model specification
  # kerlikelihood drops rows incompatible with [L, D] and returns the
  # surviving frame. Reuse it so the bootstrap, MoM seed, and reported
  # sample size all agree on which rows the fit was conditioned on.
  x <- m$x
  n <- nrow(x)
  np <- m$npars
  # parfitmom ignores dprimary and L/D; the moment-matching seed is known to
  # be biased under non-uniform primary or non-trivial truncation, but
  # Nelder-Mead recovers from moderate bias and keeping MoM primary- and
  # truncation-unaware avoids family/primary-specific corrections that would
  # add complexity disproportionate to the benefit.
  v0 <- parfitmom(x = x, family = family, incheck = FALSE,
                  L = L, D = D)$mompoint_ub
  maxs <- list(fnscale = -1)
  mle <- stats::optim(par = v0, fn = m$loglik, x = x, control = maxs)
  mlepar <- as.numeric(m$originscale(mle$par))
  mlefeat <- as.numeric(kerfeats(family = family, par = mlepar))
  mleconv <- (mle$convergence == 0)
  if(isTRUE(pgbar)) {
    cat(paste0("Fitting parametric model (", family, ") \n",
               "Bootstrap progress (Bboot=", Bboot, "): \n"))
    progbar <- utils::txtProgressBar(min = 1, max = Bboot, initial = 1,
                                     style = 3, char ="*")
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
    if(isTRUE(pgbar)) utils::setTxtProgressBar(progbar, b)
  }
  if(isTRUE(pgbar)) close(progbar)
  parfit <- stats::setNames(vector("list", np), paste0("par", 1:np))
  feats <- c("mean", "var", "sd", paste0("q", c(1, 5, 25, 50, 75, 95, 99)))
  delayfit <- stats::setNames(vector("list", length(feats)), feats)
  parfit <- kerstats(slist = parfit, pestim = mlepar, boot = pboot)
  delayfit <- kerstats(slist = delayfit, pestim = mlefeat, boot = fboot)
  aic <- 2 * (np - mle$value)
  bic <- np * log(n) - 2 * mle$value
  toc <- proc.time() - tic

  o <- c(m[!names(m) %in% c("loglik", "originscale")],
         list(n = n, Bboot = Bboot, parfit = parfit, delayfit = delayfit,
              aic = aic, bic = bic, mleconv = mleconv, bootdiscard = bootdiscard,
              elapsed = toc[3], L = L, D = D))
  attr(o, "class") <- "parfitml"
  return(o)
}
