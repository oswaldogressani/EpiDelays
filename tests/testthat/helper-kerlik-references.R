# Test-local reference oracles for kerlikelihood() equivalence tests.
#
# kerlikelihood()'s doubly-interval-censored ni branch is implemented on top
# of primarycensored::dprimarycensored(). The two oracles below reconstruct
# the same inner integral in two independent ways so the equivalence tests
# pin both the dprimarycensored idiom and the underlying integrate()-based
# definition of the inner expectation:
#
#   inner_i = (1 / (x1r - x1l)) *
#     int_{x1l}^{x1r} [F(x2r - t1) - F(x2l - t1)] dt1
#
# Both oracles are CDF-agnostic: the same family_cases entry drives them,
# keyed on whatever pdist is passed in.

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

# Idiomatic primarycensored oracle. CDF-agnostic: whatever pdist and pars
# are supplied feed straight into dprimarycensored().
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
