#' Nonparametric estimation for interval-censored data based on uniform mixtures
#'
#' @description
#' This routine uses a nonparametric methodology based on uniform mixtures
#' (Gressani and Hens, 2025) to estimate different distributional features from
#' interval-censored data. It is particularly useful in an epidemiological delay
#' modeling context to extract information from coarse delay data without
#' imposing any distributional assumption on the underlying data generating
#' mechanism. The routine requires as main input a data frame with two columns
#' indicating the left and right bound, respectively, of the interval-censored
#' observations. It provides point estimates and confidence intervals
#' (via the nonparametric bootstrap) for different features of the underlying
#' random variable being modeled. The following features are considered: mean,
#' variance, standard deviation, and 1\%, 5\%, 25\%, 50\%, 75\%, 95\% and 99\%
#' percentiles.
#'
#' @usage nonparfit(x, boot = FALSE, Bboot = 200, pgbar = FALSE)
#'
#' @param x A data frame with n rows and two columns indicating the left and
#'  right bound, respectively, of the interval-censored observations. \code{NA}
#'  values are not allowed.
#'
#' @param boot Should the nonparametric bootstrap be implemented to obtain
#'  standard errors and 90\% and 95\% confidence intervals for different
#'  features of the distribution? Default is \code{FALSE}.
#'
#' @param Bboot The number of bootstrap replications. Default is 200.
#'
#' @param pgbar Should a progress bar be displayed while the bootstrap is
#'  ongoing? Default is \code{FALSE}.
#'
#' @return A list containing different feature estimates of the underlying
#'  distribution.
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


nonparfit <- function(x, boot = FALSE, Bboot = 200, pgbar = FALSE){

n <- nrow(x)
tl <- x[, 1]
tr <- x[, 2]
ninv <- 1 / n
qfeat <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)
sstats <- data.frame(mean = numeric(0), var = numeric(0), sd = numeric(0),
                     q1 = numeric(0), q5 = numeric(0), q25 = numeric(0),
                     median = numeric(0), q75 = numeric(0), q95 = numeric(0),
                     q99 = numeric(0))

#--- Point estimation
pointestim <- function(tl, tr) {

  tmid <- 0.5 * (tl + tr)
  tw <- tr - tl

  # Nonparametric point estimates
  Fhat <- function(t) ninv * sum((t - tl) / tw * (t >= tl & t <= tr) + (t > tr))
  tord <- sort(c(tl, tr))
  Fhattord <- sapply(tord, Fhat)

  qfun <- function(p) {# Function to estimate p-quantiles
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
  qp <- sapply(qfeat, qfun)

  # Output
  outlist <- list(mu = mu, sd = sd, qp = qp)
  return(outlist)
}

thetahat <- pointestim(tl, tr)
epifeat <- sstats
epifeat[1, ] <- c(thetahat$mu, thetahat$sd^2, thetahat$sd, thetahat$qp)

#--- Nonparametric bootstrap
if(isTRUE(boot)) {
  sstatsboot <- sstats
  rowid <- seq_len(n)
  if(isTRUE(pgbar)) {
  cat(paste0("Nonparametric bootstrap \n",
             "Bootstrap sample: B=", Bboot, " \n",
             "Sample size: n=", n," \n"))
  progbar <- utils::txtProgressBar(min = 1, max = Bboot, initial = 1,
                                   style = 3, char ="*")
  }

  for (b in 1:Bboot) {
    bootrow <- sample(rowid, size = n, replace = TRUE)
    xboot <- x[bootrow, ]
    bootest <- pointestim(tl = xboot[, 1], tr = xboot[, 2])
    sstatsboot[b, ] <- c(bootest$mu, bootest$sd^2, bootest$sd, bootest$qp)
    if(isTRUE(pgbar)) utils::setTxtProgressBar(progbar, b)
  }
  if(isTRUE(pgbar)) close(progbar)
  se <- apply(sstatsboot, 2, stats::sd)
  CI90L <- apply(sstatsboot, 2, stats::quantile, prob = 0.05)
  CI90R <- apply(sstatsboot, 2, stats::quantile, prob = 0.95)
  CI95L <- apply(sstatsboot, 2, stats::quantile, prob = 0.025)
  CI95R <- apply(sstatsboot, 2, stats::quantile, prob = 0.975)
  epifeatboot <- t(rbind(epifeat, se, CI90L, CI90R, CI95L, CI95R))
  colnames(epifeatboot) <- c("estim", "se",
                             "CI90L", "CI90R", "CI95L", "CI95R")
  epifeatboot <- as.data.frame(epifeatboot)
} else{
  epifeatboot <- NULL
  Bboot <- NULL
}

#--- Routine output
outlist <- list(epifeat = epifeat, Bboot = Bboot, epifeatboot = epifeatboot)

return(outlist)
}
