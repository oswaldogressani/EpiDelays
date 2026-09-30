# Reference implementations of the doubly interval-censored likelihood:
#   (1 / (x1r - x1l)) * int_{x1l}^{x1r} [F(x2r - t) - F(x2l - t)] dt

# nolint start: object_length_linter, line_length_linter.
kerlik_integrate_reference <- function(x, pdist, pars) {
  n <- nrow(x)
  z <- 0
  for (i in seq_len(n)) {
    h <- function(t1) {
      do.call(pdist, c(list(x$x2r[i] - t1), pars)) -
        do.call(pdist, c(list(x$x2l[i] - t1), pars))
    }
    logint <- log(
      stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
    )
    z <- z + logint - log(x$x1r[i] - x$x1l[i])
  }
  z
}

kerlik_loglik_via_dprimarycensored <- function(x, pdist, pars) {
  sum(vapply(seq_len(nrow(x)), function(i) {
    do.call(
      primarycensored::dprimarycensored,
      c(
        list(
          x = x$x2l[i] - x$x1l[i],
          pdist = pdist,
          pwindow = x$x1r[i] - x$x1l[i],
          swindow = x$x2r[i] - x$x2l[i],
          log = TRUE
        ),
        pars
      )
    )
  }, numeric(1)))
}
# nolint end

# Closed form from upstream EpiDelays (0f33ae7), with G(u) = int_0^u F(s) ds.
# Falls back to integrate() where cancellation gives a non-positive value.
# nolint start: object_length_linter, line_length_linter.
kerlik_analytical_G <- function(family, pars) {
  switch(family,
    gaussian = function(u) {
      z <- (u - pars$mean) / pars$sd
      (u - pars$mean) * stats::pnorm(z) + pars$sd * stats::dnorm(z)
    },
    skewnorm = function(u) {
      z1 <- (u - pars$par1) / pars$par2
      z3 <- sqrt(1 + pars$par3^2)
      z2 <- pars$par3 / z3
      (u - pars$par1) *
        pskewnorm(u, par1 = pars$par1, par2 = pars$par2, par3 = pars$par3) +
        2 * pars$par2 * stats::dnorm(z1) * stats::pnorm(pars$par3 * z1) -
        pars$par2 * sqrt(2 / pi) * z2 * stats::pnorm(z3 * z1)
    },
    gamma = function(u) {
      u <- pmax(u, 0)
      u * stats::pgamma(u, shape = pars$shape, rate = pars$rate) -
        (pars$shape / pars$rate) *
        stats::pgamma(u, shape = pars$shape + 1, rate = pars$rate)
    },
    lognormal = function(u) {
      o <- numeric(length(u))
      sp <- u > 0
      lu <- log(u[sp])
      o[sp] <- u[sp] * stats::pnorm((lu - pars$meanlog) / pars$sdlog) -
        exp(pars$meanlog + 0.5 * pars$sdlog^2) *
        stats::pnorm((lu - pars$meanlog - pars$sdlog^2) / pars$sdlog)
      o
    },
    weibull = function(u) {
      u <- pmax(u, 0)
      z1 <- 1 + 1 / pars$shape
      u * stats::pweibull(u, shape = pars$shape, scale = pars$scale) -
        pars$scale * gamma(z1) *
        stats::pgamma((u / pars$scale)^pars$shape, shape = z1)
    }
  )
}

kerlik_analytical_reference <- function(x, family, pdist, pars) {
  G <- kerlik_analytical_G(family, pars)
  integral <- G(x$x2r - x$x1l) - G(x$x2r - x$x1r) -
    G(x$x2l - x$x1l) + G(x$x2l - x$x1r)
  for (i in which(integral <= 0)) {
    h <- function(t1) {
      do.call(pdist, c(list(x$x2r[i] - t1), pars)) -
        do.call(pdist, c(list(x$x2l[i] - t1), pars))
    }
    integral[i] <- stats::integrate(h, lower = x$x1l[i], upper = x$x1r[i])$value
  }
  log(integral) - log(x$x1r - x$x1l)
}
# nolint end
