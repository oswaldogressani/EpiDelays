#' Nonparametric estimation for single interval-censored data based on uniform mixtures
#'
#' @description
#' This routine uses a nonparametric methodology based on uniform mixtures
#' (Gressani and Hens, 2025) to estimate different distributional features from
#' interval-censored data. It is particularly useful in an epidemiological delay
#' modeling context to extract information from coarse delay data without
#' imposing any distributional assumption on the underlying data generating
#' mechanism. The routine requires as main input a data frame \code{x} with two
#' columns (named \code{xl} and \code{xr}) indicating the left and right bound,
#' respectively, of the interval-censored observations. The naming convention of
#' data frame columns is strict and different namings are not allowed.
#' Note that the left bound should be strictly smaller than the right
#' bound, i.e. \code{xl < xr} must be satisfied for all rows in \code{x}.
#' The routine computes point estimates and confidence intervals
#' (via the nonparametric bootstrap) for different features of the underlying
#' random variable being modeled. The following features are considered: mean,
#' variance, standard deviation, and 1\%, 5\%, 25\%, 50\%, 75\%, 95\% and 99\%
#' percentiles.
#'
#' @param x A data frame with n rows and two columns indicating the left and
#'  right bound, respectively, of the interval-censored observations. \code{NA}
#'  values are not allowed.
#' @param Bboot Number of bootstrap samples. Default is 1000.
#' @param pgbar Should a progress bar be displayed in console? Default is TRUE.
#'
#' @return A list containing detailed information on the fitted nonparametric
#' model. The \code{summary()} function can be used to see further details.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references Gressani, O. and Hens, N. (2025). Nonparametric serial interval
#' estimation with uniform mixtures. \emph{Plos Computational Biology},
#'  \strong{21}(8): e1013338.
#'
#' @examples
#' # Add examples here
#'
#' @export

nonparfit <- function(x, Bboot = 1000, pgbar = TRUE){
  tic <- proc.time()
  n <- nrow(x)
  nc <- ncol(x)
  if (nc != 2) {
    stop("Data frame must have 2 columns")
  }
  dfck <- kerdata_check(x = x) # data frame check
  if (dfck$result == "fail") {
    stop(dfck$message)
  }
  xl <- x[, 1]
  xr <- x[, 2]
  ninv <- 1 / n
  pfeats <- c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)
  pointestim <- function(tl, tr) {   # Nonparametric point estimation
    tmid <- 0.5 * (tl + tr)
    tw <- tr - tl
    Fhat <- function(t) ninv * sum((t - tl) / tw * (t >= tl & t <= tr) + (t > tr))
    tord <- sort(c(tl, tr))
    Fhattord <- sapply(tord, Fhat)
    qfun <- function(p) {# Estimation of p-quantiles
      pcub <- which(p <= Fhattord)[1]
      t1 <- tord[pcub - 1]
      t2 <- tord[pcub]
      Fhatt1 <- Fhat(t1)
      Fhatt2 <- Fhat(t2)
      dFinv <- 1 / (Fhatt2 - Fhatt1)
      val <- (t1 * (Fhatt2 - p) + t2 * (p - Fhatt1)) * dFinv
      return(val)
    }
    mu <- mean(tmid)
    sd <- sqrt(mean((tl^2 + tl * tr + tr^2) / 3) - mu^2)
    qp <- sapply(pfeats, qfun)
    o <- list(mu = mu, sd = sd, qp = qp)
    return(o)
  }
  npp <- pointestim(xl, xr)
  npfeat <- c(npp$mu, npp$sd^2, npp$sd, npp$qp)
  if(isTRUE(pgbar)) {
    cat(paste0("Nonparametric fit \n",
               "Bootstrap progress (Bboot=", Bboot, "): \n"))
    progbar <- utils::txtProgressBar(min = 1, max = Bboot, initial = 1,
                                     style = 3, char ="*")
  }
  fboot <- matrix(0, nrow = Bboot, ncol = length(npfeat))
  for (b in 1:Bboot) {
    xb <- kerboot(x)
    bootest <- pointestim(tl = xb[, 1], tr = xb[, 2])
    fboot[b, ] <- c(bootest$mu, bootest$sd^2, bootest$sd, bootest$qp)
    if(isTRUE(pgbar)) utils::setTxtProgressBar(progbar, b)
  }
  if(isTRUE(pgbar)) close(progbar)
  feats <- c("mean", "var", "sd", paste0("q", c(1, 5, 25, 50, 75, 95, 99)))
  delayfit <- stats::setNames(vector("list", length(feats)), feats)
  delayfit <- kerstats(slist = delayfit, pestim = npfeat, method = "boot",
                       boot = fboot)
  toc <- proc.time() - tic
  o <- list(x = x, n = n, Bboot = Bboot, delayfit = delayfit,
            censtype = "single", elapsed = toc[3])
  attr(o, "class") <- "nonparfit"
  return(o)
}
