#' Nonparametric serial interval estimation
#'
#' @description
#' This routine uses a nonparametric methodology to estimate the serial interval
#' distribution based on interval censored serial interval data. It requires as
#' input a data frame with two columns (the left and right bound of the obseved
#' serial interval window). The routine outputs point estimates and confidence
#' intervals for selected features of the serial interval that are often reported.
#'
#'
#' @usage estimSI(x, nboot = 2000)
#'
#' @param x A data frame with n rows (number of transmission pairs) and two
#'  columns: the left and right bound, respectively, of the observed serial
#'  interval window.
#' @param nboot The bootstrap sample size to construct confidence intervals.
#'
#' @details The \code{estimSI} routine provides point estimates and confidence
#' intervals for the mean, standard deviation, and 5th, 25th,
#' 50th (median), 75th and 95th quantiles of the serial interval distribution.
#'
#' @return A list containing the parameter estimates.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @examples
#' # Add examples here
#' # estimSI...
#'
#' @export


estimSI <- function(x, nboot = 2000){

n <- nrow(x)
sl <- x[, 1]
sr <- x[, 2]
ninv <- 1 / n
p <- c(0.05,0.25,0.50,0.75,0.95) # quantile evaluation
pairidx <- seq_len(n)

#---------- Point estimation

pointestim <- function(sl, sr) {

  smid <- 0.5 * (sl + sr)
  sw <- sr - sl

  # Nonparametric point estimates
  Fhat <- function(s) ninv * sum((s - sl) / sw * (s >= sl & s <= sr) + (s > sr))
  sord <- sort(c(sl, sr))
  Fhatsord <- sapply(sord, Fhat)

  qfun <- function(p) {# Function to estimate p-quantiles
    pcub <- which(p <= Fhatsord)[1]
    s1 <- sord[pcub - 1]
    s2 <- sord[pcub]
    Fhats1 <- Fhat(s1)
    Fhats2 <- Fhat(s2)
    dFinv <- 1 / (Fhats2 - Fhats1)
    val <- (s1 * (Fhats2 - p) + s2 * (p - Fhats1)) * dFinv
    return(val)
  }

  mu <- mean(smid)
  sd <- sqrt(mean((sl^2 + sl * sr + sr^2) / 3) - mu^2)
  qp <- sapply(p, qfun)

  # Output
  outlist <- list(mu = mu, sd = sd, qp = qp)
  return(outlist)
}

thetahat <- pointestim(sl, sr)

#---------- Nonparametric bootstrap
mu_boot <- c()
sd_boot <- c()
qp_boot <- matrix(0, nrow = nboot, ncol = length(p))

for(b in 1:nboot) {
  resample <- sample(pairidx, size = n, replace = TRUE)
  xboot <- data.frame(as.matrix(x)[resample, ])
  slboot <- xboot[, 1]
  srboot <- xboot[, 2]
  boothat <- pointestim(sl = slboot, sr = srboot)
  mu_boot[b] <- boothat$mu
  sd_boot[b] <- boothat$sd
  qp_boot[b, ] <- boothat$qp
}

#---------- Nonparametric results
npestim <- matrix(0, nrow = 6, ncol = length(p) + 2)
rownames(npestim) <- c("point", "se", "ci90l", "ci90r", "ci95l", "ci95r")
colnames(npestim) <- c("mean", "sd", paste0("q",p))
npestim[,1] <- c(thetahat$mu, stats::sd(mu_boot),
                 stats::quantile(mu_boot,
                                 probs = c(0.05, 0.95, 0.025, 0.975)))
npestim[,2] <- c(thetahat$sd, stats::sd(sd_boot),
                 stats::quantile(sd_boot,
                                 probs = c(0.05, 0.95, 0.025, 0.975)))
for(j in 1:length(p)){
npestim[,j + 2] <- c(thetahat$qp[j], stats::sd(qp_boot[,j]),
                 stats::quantile(qp_boot[,j],
                                 probs = c(0.05, 0.95, 0.025, 0.975)))
}
npestim <- data.frame(npestim)

# Output estimates in list
outlist <- list(x = x, n = n, nboot = nboot, npestim = npestim)

return(outlist)
}
