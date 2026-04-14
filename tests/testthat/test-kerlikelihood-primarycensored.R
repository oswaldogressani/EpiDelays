# Equivalence tests: kerlikelihood() vs primarycensored.
#
# kerlikelihood()'s doubly-interval-censored ni branch is implemented on top
# of primarycensored::dprimarycensored(). These tests pin that equivalence
# down numerically and lock the current implementation against a direct
# integrate()-based reconstruction of the same inner integral, so that
# future refactors cannot drift silently.
#
# Derivation:
#
#   kerlikelihood inner integral (per row i):
#     (1 / (x1r - x1l)) * int_{x1l}^{x1r} [F(x2r - t1) - F(x2l - t1)] dt1
#
# Substitute s = t1 - x1l, so s ~ Uniform(0, pwindow) with pwindow = x1r - x1l:
#     = E_s[ F((x2r - x1l) - s) - F((x2l - x1l) - s) ]
#
# By definition dprimarycensored(x2l - x1l, pwindow, swindow = x2r - x2l)
# returns the log of that difference, so the total log-likelihood is a single
# vectorised call to dprimarycensored with log = TRUE.

make_double_data <- function(seed = 1L, n = 5L) {
  set.seed(seed)
  # Continuous primary onset calendar baseline plus row-varying primary and
  # secondary observation windows. Delays are drawn via
  # primarycensored::rprimarycensored() so the simulator demonstrates the
  # primarycensored idiom on the sampling side. rprimarycensored quantises
  # the total delay to a multiple of swindow, so the secondary observation
  # interval is [x1l + delay, x1l + delay + swidth]. The evaluation family in
  # the tests below is not required to match the sampling distribution — the
  # simulator just needs to produce plausible doubly-interval-censored rows.
  x1l <- sort(stats::runif(n, 0, 3))
  pwidth <- stats::runif(n, 0.5, 1.5)
  swidth <- stats::runif(n, 0.3, 1.5)
  delays <- vapply(seq_len(n), function(i) {
    primarycensored::rprimarycensored(
      n = 1L,
      rdist = stats::rnorm,
      pwindow = pwidth[i],
      swindow = swidth[i],
      mean = 2.5,
      sd = 1
    )
  }, numeric(1))
  data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )
}

# Reference oracles (kerlik_loglik_via_dprimarycensored,
# kerlik_integrate_reference) are defined in
# tests/testthat/helper-kerlik-references.R. They are CDF-agnostic, so each
# family_cases entry supplies its own `pdist`. Each entry uses the same raw
# CDF that kerlikelihood() itself passes into dprimarycensored, so the
# production path and the oracles share the same underlying delay CDF.
family_cases <- list(
  gaussian = list(
    v = c(1.5, log(0.8)),
    pdist = stats::pnorm,
    pars = list(mean = 1.5, sd = 0.8)
  ),
  gamma = list(
    v = log(c(2, 0.5)),
    pdist = stats::pgamma,
    pars = list(shape = 2, rate = 0.5)
  ),
  lognormal = list(
    v = c(0.5, log(0.4)),
    pdist = stats::plnorm,
    pars = list(meanlog = 0.5, sdlog = 0.4)
  ),
  weibull = list(
    v = log(c(2, 2)),
    pdist = stats::pweibull,
    pars = list(shape = 2, scale = 2)
  ),
  skewnorm = list(
    v = c(1, log(1), 2),
    # pskewnorm uses `x` as its first argument; primarycensored calls pdist
    # with the standard `q = ...` keyword, so the oracle wraps it the same
    # way kerlikelihood() does internally.
    pdist = function(q, par1, par2, par3) {
      pskewnorm(x = q, par1 = par1, par2 = par2, par3 = par3)
    },
    pars = list(par1 = 1, par2 = 1, par3 = 2)
  )
)

for (fam in names(family_cases)) {
  local({
    family <- fam
    case <- family_cases[[family]]

    test_that(sprintf(
      "%s double-interval ni matches dprimarycensored oracle", family
    ), {
      skip_if_no_primarycensored()
      x <- make_double_data()

      m <- kerlikelihood(x = x, family = family, likapprox = "ni")
      ker_value <- m$loglik(case$v, x)

      expected <- kerlik_loglik_via_dprimarycensored(
        x = x, pdist = case$pdist, pars = case$pars
      )

      expect_equal(ker_value, expected, tolerance = 1e-8)
    })

    test_that(sprintf(
      "%s double-interval mc matches dprimarycensored oracle", family
    ), {
      skip_if_no_primarycensored()
      x <- make_double_data()

      expected <- kerlik_loglik_via_dprimarycensored(
        x = x, pdist = case$pdist, pars = case$pars
      )

      # kerlikelihood's mc closure uses stats::runif internally; seeding the
      # global RNG immediately before the call pins the Monte Carlo draws and
      # keeps the test reproducible. With M = 1000 samples per row and 5 rows
      # the summed log-likelihood sits well within 5e-2 of the ni target for
      # all families tested here.
      m <- kerlikelihood(x = x, family = family, likapprox = "mc")
      set.seed(20260413L)
      ker_value <- m$loglik(case$v, x)

      expect_equal(ker_value, expected, tolerance = 5e-2)
    })

    test_that(sprintf(
      "%s ni matches the integrate-based reference", family
    ), {
      # Regression lock: the production dprimarycensored path must agree
      # with the direct integrate()-based reconstruction in
      # kerlik_integrate_reference(), flagging any drift in a future
      # refactor. The oracle calls the same raw CDF that kerlikelihood
      # forwards to dprimarycensored, so both sides see identical underlying
      # delay semantics.
      skip_if_no_primarycensored()
      x <- make_double_data()

      m <- kerlikelihood(x = x, family = family, likapprox = "ni")
      expect_equal(
        m$loglik(case$v, x),
        kerlik_integrate_reference(
          x = x, pdist = case$pdist, pars = case$pars
        ),
        tolerance = 1e-8
      )
    })

    test_that(sprintf(
      "%s single-interval matches F(xr) - F(xl)", family
    ), {
      skip_if_no_primarycensored()
      x_double <- make_double_data()
      # Reuse the existing interval for the secondary event as a single
      # interval-censored observation (xl, xr) — positive throughout so the
      # gamma/lognormal/weibull supports are respected.
      x <- data.frame(xl = x_double$x2l, xr = x_double$x2r)

      m <- kerlikelihood(x = x, family = family, likapprox = "ni")
      ker_value <- m$loglik(case$v, x)

      Fl <- do.call(case$pdist, c(list(x$xl), case$pars))
      Fr <- do.call(case$pdist, c(list(x$xr), case$pars))
      expected <- sum(log(Fr - Fl))

      expect_equal(ker_value, expected, tolerance = 1e-10)
    })
  })
}

test_that("parfitml fits gamma end-to-end on doubly censored data", {
  skip_if_no_primarycensored()
  set.seed(42L)
  n <- 30L
  # Simulate doubly interval-censored gamma data via
  # primarycensored::rprimarycensored(). Narrow primary and secondary windows
  # keep the fit quick with few bootstrap iterations. rprimarycensored needs
  # scalar pwindow/swindow, so draw row-by-row.
  x1l <- sort(stats::runif(n, 0, 2))
  pwidth <- stats::runif(n, 0.2, 0.8)
  swidth <- stats::runif(n, 0.2, 0.8)
  delays <- vapply(seq_len(n), function(i) {
    primarycensored::rprimarycensored(
      n = 1L,
      rdist = stats::rgamma,
      pwindow = pwidth[i],
      swindow = swidth[i],
      shape = 3,
      rate = 1
    )
  }, numeric(1))
  x <- data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )

  fit <- parfitml(x = x, family = "gamma", Bboot = 10L, pgbar = FALSE)
  expect_true(fit$mleconv)
  expect_identical(fit$censtype, "double")
})

test_that("kerlikelihood handles mixed pwindow rows", {
  skip_if_no_primarycensored()
  # Mix two distinct primary windows (1 and 2) so the group-by logic in
  # build_pc_loglik() has to dispatch more than one dprimarycensored call.
  x <- data.frame(
    x1l = c(0, 1, 2, 3, 4, 5),
    x1r = c(1, 2, 3, 5, 6, 7),
    x2l = c(3, 4, 5, 6, 7, 8),
    x2r = c(4, 5, 6, 7, 8, 9)
  )
  stopifnot(length(unique(x$x1r - x$x1l)) == 2L)

  v <- c(1.5, log(0.8))
  m <- kerlikelihood(x = x, family = "gaussian", likapprox = "ni")

  expected <- kerlik_loglik_via_dprimarycensored(
    x = x, pdist = stats::pnorm,
    pars = list(mean = v[1], sd = exp(v[2]))
  )
  expect_equal(m$loglik(v, x), expected, tolerance = 1e-8)
})
