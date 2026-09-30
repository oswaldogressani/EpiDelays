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
#' inner integral. By default, this is implemented via numerical integration,
#' where the integral is evaluated with the \code{integrate} routine. This
#' routine also returns a function (logliki) that evaluates the log-likelihood
#' pointwise at each datum.
#'
#' @param x A data frame with either two columns named \code{xl} and \code{xr},
#' or four columns named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}. See
#' description for constraints imposed on the columns.
#' @param family A character string specifying the name of the parametric
#' family.
#'
#' @return A list containing information on the chosen parametric family,
#' the log-likelihood function, a function that transforms back the parameters
#' in their original scale, and the censoring type. It also
#' returns the Jacobian of the inverse of the function that transforms a
#' bounded domain to an unbounded domain.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerlikelihood <- function(x, family) {
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
    xmin <- min(x$xl)
    xmax <- max(x$xr)
  } else if (nc == 4) {
    censtype <- "double"
    d1 <- x$x2r - x$x1l
    d2 <- x$x2r - x$x1r
    d3 <- x$x2l - x$x1l
    d4 <- x$x2l - x$x1r
    xmin <- min(d4)
    xmax <- max(d1)
    logdx1 <- log(x$x1r - x$x1l)
  }
  if (family == "gaussian") {
    if(nc == 2) {
      logliki <- function(v, x) {
        par1 <- v[1]
        par2 <- exp(v[2])
        Fl  <- stats::pnorm(q = x$xl, mean = par1, sd = par2)
        Fr  <- stats::pnorm(q = x$xr, mean = par1, sd = par2)
        z <- log(Fr - Fl)
        return(z)
      }
    } else if(nc == 4) {
      logliki <- function(v, x){
        par1 <- v[1]
        par2 <- exp(v[2])
        G <- function(u){
          z <- (u - par1) / par2
          o <- (u - par1) * stats::pnorm(z) + par2 * stats::dnorm(z)
          return(o)
        }
        I <- G(d1) - G(d2) - G(d3) + G(d4)
        intfback <- which(I <= 0)
        if(length(intfback) > 0L){
          for(i in intfback) {
            h <- function(x1) {
              Fl  <- stats::pnorm(q = x$x2l[i] - x1 , mean = par1, sd = par2)
              Fr  <- stats::pnorm(q = x$x2r[i] - x1 , mean = par1, sd = par2)
              hval <- Fr - Fl
              return(hval)
            }
            I[i]<- stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
          }
        }
        z <- log(I) - logdx1
        return(z)
      }
    }
    originscale <- function(v) {
      z <- data.frame(v[1], exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      return(z)
    }
    J <- function(v) {
      o <- diag(c(1, exp(v[2])))
      return(o)
    }
  } else if (family == "skewnorm") {
    if(nc == 2) {
      logliki <- function(v, x) { # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        par3 <- v[3]
        Fl <- pskewnorm(q = x$xl, par1 = par1, par2 = par2, par3 = par3)
        Fr <- pskewnorm(q = x$xr, par1 = par1, par2 = par2, par3 = par3)
        z   <- log(Fr - Fl)
        return(z)
      }
    } else if(nc == 4) {
      rpi <- sqrt(2 / pi)
      logliki <- function(v, x) {
          par1 <- v[1]
          par2 <- exp(v[2])
          par3 <- v[3]
          G <- function(u){
            z1 <- (u - par1) / par2
            z3 <- sqrt(1 + par3^2)
            z2 <- par3 / z3
            o <- (u - par1) * pskewnorm(u, par1 = par1, par2 = par2, par3 = par3) +
              2 * par2 * stats::dnorm(z1) * stats::pnorm(par3 * z1) -
              par2 * rpi * z2 * stats::pnorm(z3 * z1)
            return(o)
          }
          I <- G(d1) - G(d2) - G(d3) + G(d4)
          intfback <- which(I <= 0)
          if(length(intfback) > 0L){
            for(i in intfback) {
              h <- function(x1) {
              Fl <- pskewnorm(q = x$x2l[i] - x1, par1 = par1, par2 = par2, par3 = par3)
              Fr <- pskewnorm(q = x$x2r[i] - x1, par1 = par1, par2 = par2, par3 = par3)
              hval <- Fr - Fl
              return(hval)
              }
              I[i] <- stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
            }
          }
          z <- log(I) - logdx1
          return(z)
      }
    }
    originscale <- function(v) {
      z <- data.frame(v[1], exp(v[2]), v[3])
      colnames(z) <- c(famdesc$par1, famdesc$par2, famdesc$par3)
      return(z)
    }
    J <- function(v) {
      o <- diag(c(1, exp(v[2]), 1))
      return(o)
    }
  } else if (family == "gamma") {
    if(nc == 2) {
      logliki <- function(v, x) {
        par1 <- exp(v[1])
        par2 <- exp(v[2])
        Fl  <- stats::pgamma(q = x$xl, shape = par1, rate = par2)
        Fr  <- stats::pgamma(q = x$xr, shape = par1, rate = par2)
        z   <- log(Fr - Fl)
        return(z)
      }
    } else if(nc == 4) {
      logliki <- function(v, x) { # v: unbounded parameter
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          G <- function(u){
            o <- numeric(length(u))
            sp <- (u > 0)
            if(any(sp)){
              usp <- u[sp]
              o[sp] <- usp * stats::pgamma(usp, shape = par1, rate = par2) -
                (par1/par2) * stats::pgamma(usp, shape = par1 + 1, rate = par2)
            }
            return(o)
          }
          I <- G(d1) - G(d2) - G(d3) + G(d4)
          intfback <- which(I <= 0)
          if(length(intfback) > 0L){
            for(i in intfback) {
              h <- function(x1) {
                Fl  <- stats::pgamma(q = x$x2l[i] - x1, shape = par1, rate = par2)
                Fr  <- stats::pgamma(q = x$x2r[i] - x1, shape = par1, rate = par2)
                hval <- Fr - Fl
                return(hval)
              }
              I[i] <- stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
            }
          }
          z <- log(I) - logdx1
          return(z)
      }
    }
    originscale <- function(v) {
      z <- data.frame(exp(v[1]), exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      return(z)
    }
    J <- function(v) {
      o <- diag(c(exp(v[1]), exp(v[2])))
      return(o)
    }
  } else if (family == "lognormal") {
    if(nc == 2) {
      logliki <- function(v, x) { # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        Fl  <- stats::plnorm(q = x$xl, meanlog = par1, sdlog = par2)
        Fr  <- stats::plnorm(q = x$xr, meanlog = par1, sdlog = par2)
        z   <- log(Fr - Fl)
        return(z)
      }
    } else if(nc == 4) {
      logliki <- function(v, x) {
          par1 <- v[1]
          par2 <- exp(v[2])
          G <- function(u){
            o <- numeric(length(u))
            sp <- (u > 0)
            if(any(sp)){
              usp <- u[sp]
              z1 <- (log(usp) - par1) / par2
              z2 <- (log(usp) - par1 - par2^2) / par2
              o[sp] <- usp * stats::pnorm(z1) - exp(par1 + 0.5 * par2^2) *
                stats::pnorm(z2)
            }
            return(o)
          }
          I <- G(d1) - G(d2) - G(d3) + G(d4)
          intfback <- which(I <= 0)
          if(length(intfback) > 0L){
            for(i in intfback) {
              h <- function(x1) {
                Fl  <- stats::plnorm(q = x$x2l[i] - x1, meanlog = par1, sdlog = par2)
                Fr  <- stats::plnorm(q = x$x2r[i] - x1, meanlog = par1, sdlog = par2)
                hval <- Fr - Fl
                return(hval)
              }
              I[i] <- stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
            }
          }
          z <- log(I) - logdx1
          return(z)
      }
    }
    originscale <- function(v) {
      z <- data.frame(v[1], exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      return(z)
    }
    J <- function(v) {
      o <- diag(c(1, exp(v[2])))
      return(o)
    }
  } else if (family == "weibull") {
    if(nc == 2) {
      logliki <- function(v, x) { # v: unbounded parameter
        par1 <- exp(v[1])
        par2 <- exp(v[2])
        Fl  <- stats::pweibull(q = x$xl, shape = par1, scale = par2)
        Fr  <- stats::pweibull(q = x$xr, shape = par1, scale = par2)
        z   <- log(Fr - Fl)
        return(z)
      }
    } else if(nc == 4) {
      logliki <- function(v, x) {
          par1 <- exp(v[1])
          par2 <- exp(v[2])
          G <- function(u){
            o <- numeric(length(u))
            sp <- (u > 0)
            if(any(sp)){
              usp <- u[sp]
              z1 <- 1 + 1 / par1
              z2 <- (usp / par2)^par1
              o[sp] <- usp * stats::pweibull(usp, shape = par1, scale = par2) -
                par2 * gamma(z1) * stats::pgamma(z2, shape = z1)
            }
            return(o)
          }
          I <- G(d1) - G(d2) - G(d3) + G(d4)
          intfback <- which(I <= 0)
          if(length(intfback) > 0L){
            for(i in intfback) {
               h <- function(x1) {
                Fl  <- stats::pweibull(q = x$x2l[i] - x1, shape = par1, scale = par2)
                Fr  <- stats::pweibull(q = x$x2r[i] - x1, shape = par1, scale = par2)
                hval <- Fr - Fl
                return(hval)
              }
              I[i] <- stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
            }
          }
          z <- log(I) - logdx1
          return(z)
      }
    }
    originscale <- function(v) {
      z <- data.frame(exp(v[1]), exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      return(z)
    }
    J <- function(v) {
      o <- diag(c(exp(v[1]), exp(v[2])))
      return(o)
    }
  }
  loglik <- function(v, x) {
    o <- sum(logliki(v, x))
    return(o)
  }
  o <- c(famdesc, list(logliki = logliki, loglik = loglik,
                       originscale = originscale, censtype = censtype, J = J,
                       xmin = xmin, xmax = xmax))
  return(o)
}


