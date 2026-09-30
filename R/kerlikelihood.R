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
#' The doubly interval-censored likelihood function marginalises out the
#' primary event time. This marginalisation is delegated to
#' \code{primarycensored::dprimarycensored()}, which uses an analytical
#' solution where one is available for the given CDF and primary event
#' distribution and numerical integration otherwise. This routine also
#' returns a function (logliki) that evaluates the log-likelihood pointwise
#' at each datum.
#'
#' @param x A data frame with either two columns named \code{xl} and \code{xr},
#' or four columns named \code{x1l}, \code{x1r}, \code{x2l}, \code{x2r}. See
#' description for constraints imposed on the columns.
#' @param family A character string specifying the name of the parametric
#' family.
#' @param L Lower truncation point applied to the underlying delay
#' distribution. Defaults to \code{-Inf} (no left truncation), matching
#' \code{primarycensored}. When \code{L} is finite the per-row
#' contribution is rescaled by the truncated CDF.
#' @param D Upper truncation point applied to the underlying delay
#' distribution. Defaults to \code{Inf}. Together with \code{L}, rows whose
#' observed interval is incompatible with \code{[L, D]} are rejected at
#' input time.
#' @param dprimary Primary event density function. Used only when \code{x} has
#' four columns (doubly interval-censored data).
#' Must take a vector \code{x} and the arguments \code{min} and \code{max},
#' return a density normalised to integrate to 1 on \code{[min, max]}, and be
#' compatible with \code{primarycensored::dprimarycensored()}. Defaults to
#' \code{stats::dunif} (uniform primary onset within the observation window),
#' which reproduces the behaviour of earlier EpiDelays versions. Non-uniform
#' choices include \code{primarycensored::dexpgrowth} for exponential growth
#' during an outbreak. Ignored when \code{x} has two columns (single
#' interval-censored data) but must equal the default in that case.
#' @param dprimary_args A named list of additional arguments passed to
#' \code{dprimary}, mirroring the \code{dprimary_args} argument of
#' \code{primarycensored::dprimarycensored()}. Example: \code{list(r = 0.1)}
#' for \code{primarycensored::dexpgrowth}. Defaults to an empty list. Must be
#' empty unless \code{x} has four columns.
#'
#' @return A list containing information on the chosen parametric family,
#' the log-likelihood function, its pointwise version, a function that
#' transforms back the parameters in their original scale, the censoring
#' type, and the primary event density settings. It also returns the
#' Jacobian of the inverse of the function that transforms a bounded domain
#' to an unbounded domain.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerlikelihood <- function(x, family, # nolint: cyclocomp_linter.
                          L = -Inf, D = Inf,
                          dprimary = stats::dunif,
                          dprimary_args = list()) {
  # Truncation bound validation mirrors
  # primarycensored::.check_truncation_bounds.
  if (!is.numeric(L) || length(L) != 1L || is.na(L)) {
    stop("L must be a numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(D) || length(D) != 1L || is.na(D) || L >= D) {
    stop("L must be less than D.", call. = FALSE)
  }
  if (!is.function(dprimary)) {
    stop("dprimary must be a function", call. = FALSE)
  }
  if (!is.list(dprimary_args) ||
     (length(dprimary_args) > 0 && is.null(names(dprimary_args)))) {
    stop("dprimary_args must be a named list", call. = FALSE)
  }
  dprimary_default <- identical(dprimary, stats::dunif) &&
    length(dprimary_args) == 0
  # Input checks
  dfck <- kerdata_check(x = x) # data frame check
    if (dfck$result == "fail") {
      stop(dfck$message, call. = FALSE)
    }
  famck <- kerfamily_check(x = family) # family check
    if (famck$result == "fail") {
      stop(famck$message, call. = FALSE)
    }
  domck <- kerdomain_check(x = x, family = family) # Domain check
  if (domck$result == "fail") {
    stop(domck$message, call. = FALSE)
  }
  fset <- kerfamilies()
  fnames <- sapply(fset, "[[", "fname")
  famdesc <- fset[[match(family, fnames)]]
  nc <- ncol(x)
  if (nc == 2) {
    censtype <- "single"
    # Non-uniform primary has no meaning without a primary event window.
    if (!dprimary_default) {
      stop(
        "dprimary only applies to doubly interval-censored data (four ",
        "columns); x has two columns so no primary event window is ",
        "modelled",
        call. = FALSE
      )
    }
    # Drop any row whose observed interval is not fully inside [L, D]. The
    # interval-censored likelihood treats each row's window as an indivisible
    # quantum, so a row that straddles a truncation boundary cannot be
    # interpreted under the truncated model. Warning the caller and removing
    # the row keeps the fit going while making the loss explicit; the caller
    # can re-run after narrowing the offending intervals.
    if (is.finite(L) || is.finite(D)) {
      keep <- x$xl >= L & x$xr <= D
      if (!all(keep)) {
        warning(
          sum(!keep), " row(s) of x straddle or fall outside the truncation ",
          "bounds [L, D] and have been dropped. Each row's interval ",
          "[xl, xr] must satisfy xl >= L and xr <= D. To retain these ",
          "observations, narrow their intervals before calling the fit.",
          call. = FALSE
        )
        x <- x[keep, , drop = FALSE]
        if (nrow(x) == 0L) {
          stop(
            "No rows of x remain after dropping observations incompatible ",
            "with the truncation bounds [L, D].",
            call. = FALSE
          )
        }
      }
    }
    xmin <- min(x$xl)
    xmax <- max(x$xr)
  } else if (nc == 4) {
    censtype <- "double"
    # Drop rows whose secondary observation window is not fully inside
    # [L, D]. The doubly-interval-censored likelihood passes the row's
    # (lower, swindow) pair into primarycensored::dprimarycensored() with
    # truncation bounds L and D; primarycensored requires the entire window
    # to sit inside [L, D] and aborts otherwise. Rather than splitting a
    # straddling window into a visible subset (which silently changes the
    # observation), warn the caller and drop the row so the modelled and
    # observed windows match.
    if (is.finite(L) || is.finite(D)) {
      lowers <- x$x2l - x$x1l
      uppers <- x$x2r - x$x1l
      keep <- lowers >= L & uppers <= D
      if (!all(keep)) {
        warning(
          sum(!keep), " row(s) of x straddle or fall outside the truncation ",
          "bounds [L, D] and have been dropped. Each row's secondary window ",
          "must satisfy x2l - x1l >= L and x2r - x1l <= D. To retain these ",
          "observations, narrow their secondary windows before calling the ",
          "fit.",
          call. = FALSE
        )
        x <- x[keep, , drop = FALSE]
        if (nrow(x) == 0L) {
          stop(
            "No rows of x remain after dropping observations incompatible ",
            "with the truncation bounds [L, D].",
            call. = FALSE
          )
        }
      }
    }
    # Shortest and longest delays compatible with the observed windows.
    xmin <- min(x$x2l - x$x1r)
    xmax <- max(x$x2r - x$x1l)
  }
  # Helper: per-row log-contributions for the single-interval (nc == 2)
  # branch, applying the truncation correction when L or D is finite.
  # Reduces to log(Fr - Fl) in the default case.
  single_interval_i <- function(Fl, Fr, FL, FD) {
    if (is.infinite(L) && is.infinite(D)) {
      return(log(Fr - Fl))
    }
    num <- pmin(Fr, FD) - pmax(Fl, FL)
    log(num) - log(FD - FL)
  }
  # Pointwise loglik for doubly interval-censored data via primarycensored.
  # Rows are grouped by (pwindow, swindow) and evaluated at each group's
  # unique lower bounds. The grouping is cached per x as optim reuses x.
  build_pc_logliki <- function(pdist, pars_fn) {
    force(pdist)
    force(pars_fn)
    state <- new.env(parent = emptyenv())
    prepare <- function(x) {
      pwindows <- x$x1r - x$x1l
      swindows <- x$x2r - x$x2l
      lowers <- x$x2l - x$x1l
      upw <- unique(pwindows)
      usw <- unique(swindows)
      key <- match(pwindows, upw) +
        length(upw) * (match(swindows, usw) - 1L)
      lapply(split(seq_along(key), key), function(r) {
        lower <- lowers[r]
        ulower <- unique(lower)
        list(
          rows = r, pwindow = pwindows[r[1]], swindow = swindows[r[1]],
          lower = ulower, map = match(lower, ulower)
        )
      })
    }
    function(v, x) {
      pars <- pars_fn(v)
      if (!identical(x, state$x)) {
        assign("prep", prepare(x), envir = state)
        assign("x", x, envir = state)
      }
      if (is.null(state$pcens)) {
        # Check pdist and dprimary once, then reuse the object.
        do.call(primarycensored::check_pdist, c(list(pdist, D = D), pars))
        for (g in state$prep) {
          primarycensored::check_dprimary(dprimary, g$pwindow, dprimary_args)
        }
        obj <- do.call(
          primarycensored::new_pcens,
          c(
            list(
              pdist = pdist, dprimary = dprimary,
              primary_args = dprimary_args
            ),
            pars
          )
        )
      } else {
        obj <- do.call(stats::update, c(list(state$pcens), pars))
      }
      assign("pcens", obj, envir = state)
      z <- numeric(nrow(x))
      for (g in state$prep) {
        z[g$rows] <- primarycensored::pcens_pmf(
          obj, g$lower, g$pwindow,
          swindow = g$swindow, L = L, D = D, log = TRUE
        )[g$map]
      }
      z
    }
  }
  if (family == "gaussian") { # nolint: if_switch_linter.
    # Gaussian uses the raw stats::pnorm CDF on the full real line. With the
    # straddle-row drop above, primarycensored::dprimarycensored() sees only
    # rows whose secondary window sits inside [L, D] and applies its own
    # truncation correction via L and D.
    if (nc == 2) {
      logliki <- function(v, x) {
        par1 <- v[1]
        par2 <- exp(v[2])
        Fl  <- stats::pnorm(q = x$xl, mean = par1, sd = par2)
        Fr  <- stats::pnorm(q = x$xr, mean = par1, sd = par2)
        FL  <- stats::pnorm(q = L, mean = par1, sd = par2)
        FD  <- stats::pnorm(q = D, mean = par1, sd = par2)
        single_interval_i(Fl, Fr, FL, FD)
      }
    } else if (nc == 4) {
      logliki <- build_pc_logliki(
        pdist = stats::pnorm,
        pars_fn = function(v) list(mean = v[1], sd = exp(v[2]))
      )
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
    # Skewnorm uses the package-level pskewnorm() CDF, which is internally
    # clamped to [0, 1] and monotonised over the order of q. Owen's T can
    # otherwise return values a hair outside [0, 1] in saturating regions
    # that optim sometimes visits, and check_pdist() inside primarycensored
    # would then abort the fit.
    if (nc == 2) {
      logliki <- function(v, x) { # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        par3 <- v[3]
        Fl <- pskewnorm(q = x$xl, par1 = par1, par2 = par2, par3 = par3)
        Fr <- pskewnorm(q = x$xr, par1 = par1, par2 = par2, par3 = par3)
        FL <- pskewnorm(q = L, par1 = par1, par2 = par2, par3 = par3)
        FD <- pskewnorm(q = D, par1 = par1, par2 = par2, par3 = par3)
        single_interval_i(Fl, Fr, FL, FD)
      }
    } else if (nc == 4) {
      logliki <- build_pc_logliki(
        pdist = pskewnorm,
        pars_fn = function(v) {
          list(par1 = v[1], par2 = exp(v[2]), par3 = v[3])
        }
      )
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
    if (nc == 2) {
      logliki <- function(v, x) {
        par1 <- exp(v[1])
        par2 <- exp(v[2])
        Fl  <- stats::pgamma(q = x$xl, shape = par1, rate = par2)
        Fr  <- stats::pgamma(q = x$xr, shape = par1, rate = par2)
        FL  <- stats::pgamma(q = L, shape = par1, rate = par2)
        FD  <- stats::pgamma(q = D, shape = par1, rate = par2)
        single_interval_i(Fl, Fr, FL, FD)
      }
    } else if (nc == 4) {
      logliki <- build_pc_logliki(
        pdist = stats::pgamma,
        pars_fn = function(v) list(shape = exp(v[1]), rate = exp(v[2]))
      )
    }
    originscale <- function(v) {
      z <- data.frame(exp(v[1]), exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      return(z)
    }
    J <- function(v) {
      o <- diag(exp(v[1:2]))
      return(o)
    }
  } else if (family == "lognormal") {
    if (nc == 2) {
      logliki <- function(v, x) { # v: unbounded parameter
        par1 <- v[1]
        par2 <- exp(v[2])
        Fl  <- stats::plnorm(q = x$xl, meanlog = par1, sdlog = par2)
        Fr  <- stats::plnorm(q = x$xr, meanlog = par1, sdlog = par2)
        FL  <- stats::plnorm(q = L, meanlog = par1, sdlog = par2)
        FD  <- stats::plnorm(q = D, meanlog = par1, sdlog = par2)
        single_interval_i(Fl, Fr, FL, FD)
      }
    } else if (nc == 4) {
      logliki <- build_pc_logliki(
        pdist = stats::plnorm,
        pars_fn = function(v) list(meanlog = v[1], sdlog = exp(v[2]))
      )
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
    if (nc == 2) {
      logliki <- function(v, x) { # v: unbounded parameter
        par1 <- exp(v[1])
        par2 <- exp(v[2])
        Fl  <- stats::pweibull(q = x$xl, shape = par1, scale = par2)
        Fr  <- stats::pweibull(q = x$xr, shape = par1, scale = par2)
        FL  <- stats::pweibull(q = L, shape = par1, scale = par2)
        FD  <- stats::pweibull(q = D, shape = par1, scale = par2)
        single_interval_i(Fl, Fr, FL, FD)
      }
    } else if (nc == 4) {
      logliki <- build_pc_logliki(
        pdist = stats::pweibull,
        pars_fn = function(v) list(shape = exp(v[1]), scale = exp(v[2]))
      )
    }
    originscale <- function(v) {
      z <- data.frame(exp(v[1]), exp(v[2]))
      colnames(z) <- c(famdesc$par1, famdesc$par2)
      return(z)
    }
    J <- function(v) {
      o <- diag(exp(v[1:2]))
      return(o)
    }
  }
  loglik <- function(v, x) {
    o <- sum(logliki(v, x))
    return(o)
  }
  o <- c(famdesc, list(logliki = logliki, loglik = loglik,
                       originscale = originscale, censtype = censtype, J = J,
                       xmin = xmin, xmax = xmax,
                       dprimary = dprimary, dprimary_args = dprimary_args,
                       x = x))
  return(o)
}
