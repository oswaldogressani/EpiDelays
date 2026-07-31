#' Likelihood function for single and doubly interval-censored data
#'
#' @description
#' The likelihood function plays a key role in maximum likelihood estimation and
#' Bayesian inference. This routine admits as data input a data frame \code{x}
#' with either two columns (named \code{xl} and \code{xr}) representing the left
#' and right bound, respectively, of the delay variable, or four columns
#' (named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}) representing the left
#' and right bound, respectively, of the primary and secondary events of the
#' delay variable. The naming convention of the columns is strict and different
#' namings are not allowed. When data frame \code{x} has two columns, a single
#' interval-censored likelihood function is used. Note that the left bound
#' should be strictly smaller than the right bound, i.e. \code{xl < xr} must
#' be satisfied for all rows in \code{x}. When data frame \code{x} has four
#' columns, a doubly interval-censored likelihood function is used. In that case
#' \code{x1l < x1r} and \code{x2l < x2r} must hold for all rows in \code{x}.
#' Moreover, \code{NA} values are not allowed in data frame \code{x}.
#'
#' @details
#' Returns the log-likelihood function of parameters that have been transformed
#' to live in an unbounded parameter space. This makes the
#' maximization/evaluation of the likelihood function numerically more stable.
#' The doubly interval-censored likelihood function requires to evaluate an
#' inner integral. By default, this is implemented via numerical integration
#' (ni) with \code{likapprox = "ni"}, where the integral is evaluated with the
#' \code{integrate} routine. Another option is \code{likapprox = "mc"}, which
#' uses crude Monte Carlo sampling to evaluate the inner integral in the
#' doubly interval-censored likelihood function. Note that Monte
#' Carlo sampling results in much longer computation times to evaluate the
#' likelihood function. The number of Monte Carlo samples is fixed at 1000.
#'
#' @param x A data frame with either two columns named \code{xl} and \code{xr},
#' or four columns named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}. See
#' description for constraints imposed on the columns.
#' @param family A character string specifying the name of the parametric
#' family.
#' @param likapprox Approximation of the likelihood function in case of doubly
#' interval-censored data. Default is \code{"ni"} for numerical integration.
#' Another possibility is \code{"mc"} for crude Monte Carlo sampling with 1000
#' samples.
#'
#' @return A list containing information on the chosen parametric family,
#' the log-likelihood function, a function that transforms back the parameters
#' in their original scale, the censoring type and the approximation of the
#' likelihood function used in case of doubly interval-censored data. It also
#' returns the Jacobian of the inverse of the function that transforms a
#' bounded domain to an unbounded domain.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerlikelihood <- function(x, family, likapprox = "ni") {
  if(!(likapprox %in% c("ni", "mc"))) {
    stop("likapprox should either be ni or mc")
  }
  # Input checks
  dfck <- kerdata_check(x = x) # data frame check
    if (dfck$result == "fail") {
      stop(dfck$message)
    }
  famck <- kerfamily_check(x = family) # family check
    if (famck$result == "fail") {
      stop(famck$message)
    }
  domck <- kerdomain_check(x = x, family = family) # Domain check
  if (domck$result == "fail") {
    stop(domck$message)
  }
  fset <- kerfamilies()
  fnames <- sapply(fset, "[[", "fname")
  famdesc <- fset[[match(family, fnames)]]
  n <- nrow(x)
  nc <- ncol(x)
  if(nc == 2) {
    censtype <- "single"
    likapprox <- "None"
    xmin <- min(x$xl)
    xmax <- max(x$xr)
  } else if (nc == 4) {
    censtype <- "double"
    xmin <- min(x$x2l - x$x1r)
    xmax <- max(x$x2r - x$x1r)
    logdx1 <- log(x$x1r - x$x1l)
    M <- 1000 #  Monte carlo samples for doubly interval-censored likelihood
  }
  if (family == "gaussian") {
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        Fl  <- stats::pnorm(q = x$xl, mean = par1, sd = par2)
        Fr  <- stats::pnorm(q = x$xr, mean = par1, sd = par2)
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::pnorm(q = x$x2l[i] - x1s , mean = par1, sd = par2)
            Fr  <- stats::pnorm(q = x$x2r[i] - x1s , mean = par1, sd = par2)
            z <- z + log(mean(Fr - Fl))
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            h <- function(x1) {
              Fl  <- stats::pnorm(q = x$x2l[i] - x1 , mean = par1, sd = par2)
              Fr  <- stats::pnorm(q = x$x2r[i] - x1 , mean = par1, sd = par2)
              hval <- Fr - Fl
              hval
            }
            logint <- log(stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value)
            z <- z + logint - logdx1[i]
          }
          z
        }
      }
    }
    originscale <- function(v) {
      z <- data.frame(v[1], exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      z
    }
    J <- function(v) {
      o <- diag(c(1, exp(v[2])))
      o
    }
  }else if (family == "skewnorm") {
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        par3 <- v[3]
        Fl <- pskewnorm(q = x$xl, par1 = par1, par2 = par2, par3 = par3)
        Fr <- pskewnorm(q = x$xr, par1 = par1, par2 = par2, par3 = par3)
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          par3 <- v[3]
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl <- pskewnorm(q = x$x2l[i] - x1s, par1 = par1, par2 = par2, par3 = par3)
            Fr <- pskewnorm(q = x$x2r[i] - x1s, par1 = par1, par2 = par2, par3 = par3)
            z <- z + log(mean(Fr - Fl))
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          par3 <- v[3]
          z <- 0
          for(i in 1:n) {
            h <- function(x1) {
              Fl <- pskewnorm(q = x$x2l[i] - x1, par1 = par1, par2 = par2, par3 = par3)
              Fr <- pskewnorm(q = x$x2r[i] - x1, par1 = par1, par2 = par2, par3 = par3)
              hval <- Fr - Fl
              hval
            }
            logint <- log(stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value)
            z <- z + logint - logdx1[i]
          }
          z
        }
      }
    }
    originscale <- function(v) {
      z <- data.frame(v[1], exp(v[2]), v[3])
      colnames(z) <- c(famdesc$par1, famdesc$par2, famdesc$par3)
      z
    }
    J <- function(v) {
      o <- diag(c(1, exp(v[2]), 1))
      o
    }
  } else if (family == "gamma") {
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- exp(v[1])
        par2 <- exp(v[2])
        Fl  <- stats::pgamma(q = x$xl, shape = par1, rate = par2)
        Fr  <- stats::pgamma(q = x$xr, shape = par1, rate = par2)
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::pgamma(q = x$x2l[i] - x1s, shape = par1, rate = par2)
            Fr  <- stats::pgamma(q = x$x2r[i] - x1s, shape = par1, rate = par2)
            z <- z + log(mean(Fr - Fl))
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            h <- function(x1) {
              Fl  <- stats::pgamma(q = x$x2l[i] - x1, shape = par1, rate = par2)
              Fr  <- stats::pgamma(q = x$x2r[i] - x1, shape = par1, rate = par2)
              hval <- Fr - Fl
              hval
            }
            logint <- log(stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value)
            z <- z + logint - logdx1[i]
          }
          z
        }
      }
    }
    originscale <- function(v) {
      z <- data.frame(exp(v[1]), exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      z
    }
    J <- function(v) {
      o <- diag(c(exp(v[1]), exp(v[2])))
      o
    }
  } else if (family == "lognormal") {
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        Fl  <- stats::plnorm(q = x$xl, meanlog = par1, sdlog = par2)
        Fr  <- stats::plnorm(q = x$xr, meanlog = par1, sdlog = par2)
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::plnorm(q = x$x2l[i] - x1s, meanlog = par1, sdlog = par2)
            Fr  <- stats::plnorm(q = x$x2r[i] - x1s, meanlog = par1, sdlog = par2)
            z <- z + log(mean(Fr - Fl))
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- v[1]
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            h <- function(x1) {
              Fl  <- stats::plnorm(q = x$x2l[i] - x1, meanlog = par1, sdlog = par2)
              Fr  <- stats::plnorm(q = x$x2r[i] - x1, meanlog = par1, sdlog = par2)
              hval <- Fr - Fl
              hval
            }
            logint <- log(stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value)
            z <- z + logint - logdx1[i]
          }
          z
        }
      }
    }
    originscale <- function(v) {
      z <- data.frame(v[1], exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      z
    }
    J <- function(v) {
      o <- diag(c(1, exp(v[2])))
      o
    }
  } else if (family == "weibull") {
    if(nc == 2) {
      loglik <- function(v, x) {
        # v: unbounded parameter
        par1 <- exp(v[1])
        par2 <- exp(v[2])
        Fl  <- stats::pweibull(q = x$xl, shape = par1, scale = par2)
        Fr  <- stats::pweibull(q = x$xr, shape = par1, scale = par2)
        z   <- sum(log(Fr - Fl))
        z
      }
    } else if(nc == 4) {
      if(likapprox == "mc") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            x1s <- stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])
            Fl  <- stats::pweibull(q = x$x2l[i] - x1s, shape = par1, scale = par2)
            Fr  <- stats::pweibull(q = x$x2r[i] - x1s, shape = par1, scale = par2)
            z <- z + log(mean(Fr - Fl))
          }
          z
        }
      } else if(likapprox == "ni") {
        loglik <- function(v, x) {
          # v: unbounded parameter
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          z <- 0
          for(i in 1:n) {
            h <- function(x1) {
              Fl  <- stats::pweibull(q = x$x2l[i] - x1, shape = par1, scale = par2)
              Fr  <- stats::pweibull(q = x$x2r[i] - x1, shape = par1, scale = par2)
              hval <- Fr - Fl
              hval
            }
            logint <- log(stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value)
            z <- z + logint - logdx1[i]
          }
          z
        }
      }
    }
    originscale <- function(v) {
      z <- data.frame(exp(v[1]), exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      z
    }
    J <- function(v) {
      o <- diag(c(exp(v[1]), exp(v[2])))
      o
    }
  }
  o <- c(famdesc, list(loglik = loglik, originscale = originscale,
           censtype = censtype, J = J, likapprox = likapprox, xmin = xmin,
           xmax = xmax))
  return(o)
}
