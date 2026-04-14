# Tests for non-uniform primary event distributions (dprimary / dprimary_args)
#
# EpiDelays threads the upstream primarycensored arguments `dprimary` and
# `dprimary_args` through `kerlikelihood()` and `parfitml()` so users can
# fit doubly interval-censored data with a non-uniform primary event
# density (e.g. exponential growth during an outbreak). These tests cover:
#
#   1. No-op equivalence: explicit `dprimary = stats::dunif` must give the
#      same log-likelihood as the default path across all families.
#   2. End-to-end recovery under `dexpgrowth`: parameters come back close
#      to the truth when the same primary density used to simulate is
#      passed to the fitter.
#   3. Mismatched primary is materially biased: fitting the same
#      non-uniform data with the default uniform primary returns fits
#      that differ from the matched fit by a meaningful margin, so the
#      argument is not secretly a no-op.
#   4. Error paths: single-interval data and the Monte Carlo branch both
#      reject non-default primary densities with an informative message.

make_double_data <- function(seed = 1L, n = 5L) {
  set.seed(seed)
  # EpiDelays' domain check rejects rows where x2l - x1r < 0 on
  # non-negative-support families, so narrow the primary window and draw
  # from a gamma centred well away from zero to keep every row inside
  # the permitted region for gamma / lognormal / weibull as well as the
  # real-line-support families.
  x1l <- sort(stats::runif(n, 0, 3))
  pwidth <- stats::runif(n, 0.3, 0.6)
  swidth <- stats::runif(n, 0.3, 0.6)
  delays <- vapply(seq_len(n), function(i) {
    primarycensored::rprimarycensored(
      n = 1L,
      rdist = stats::rgamma,
      pwindow = pwidth[i],
      swindow = swidth[i],
      shape = 6,
      rate = 1
    )
  }, numeric(1))
  data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )
}

# Parameter points at which to evaluate the loglik for each family.
family_points <- list(
  gaussian = c(1.5, log(0.8)),
  gamma = log(c(2, 0.5)),
  lognormal = c(0.5, log(0.4)),
  weibull = log(c(2, 2)),
  skewnorm = c(1, log(1), 2)
)

for (fam in names(family_points)) {
  local({
    family <- fam
    v <- family_points[[family]]

    test_that(sprintf(
      "%s: explicit dprimary = stats::dunif matches default path", family
    ), {
      skip_if_no_primarycensored()
      x <- make_double_data()

      m_default <- kerlikelihood(x = x, family = family, likapprox = "ni")
      m_explicit <- kerlikelihood(
        x = x, family = family, likapprox = "ni",
        dprimary = stats::dunif, dprimary_args = list()
      )

      # A single parameter evaluation pins the forwarding call; two extra
      # perturbed points guard against accidental reliance on optim's
      # starting point.
      points <- list(
        v,
        v + 0.1,
        v - 0.1
      )
      for (p in points) {
        expect_equal(
          m_default$loglik(p, x), m_explicit$loglik(p, x),
          tolerance = 1e-12
        )
      }
    })
  })
}

test_that("parfitml recovers gamma parameters under dexpgrowth primary", {
  skip_if_no_primarycensored()
  set.seed(20260413L)
  # Simulate doubly-interval-censored gamma data with an exponentially
  # growing primary onset. rprimarycensored takes the sampler counterpart
  # `rprimary` / `rprimary_args`; the matching density for fitting is
  # `dprimary` / `dprimary_args`. primarycensored supplies
  # `min = 0, max = pwindow` to both sides automatically, so
  # rprimary_args / dprimary_args must NOT include those.
  n <- 800L
  pwidth <- 3
  swidth <- 1
  r_true <- 1
  shape_true <- 10
  rate_true <- 1
  x1l <- sort(stats::runif(n, 0, 40))
  delays <- vapply(seq_len(n), function(i) {
    primarycensored::rprimarycensored(
      n = 1L,
      rdist = stats::rgamma,
      pwindow = pwidth,
      swindow = swidth,
      rprimary = primarycensored::rexpgrowth,
      rprimary_args = list(r = r_true),
      shape = shape_true,
      rate = rate_true
    )
  }, numeric(1))
  x <- data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )

  fit <- parfitml(
    x = x, family = "gamma", Bboot = 5L, pgbar = FALSE,
    dprimary = primarycensored::dexpgrowth,
    dprimary_args = list(r = r_true)
  )
  expect_true(fit$mleconv)
  expect_identical(fit$censtype, "double")
  expect_identical(fit$dprimary, primarycensored::dexpgrowth)
  expect_identical(fit$dprimary_args, list(r = r_true))
  # Parameter recovery is approximate; 25% tolerance is loose enough to
  # absorb sampling noise at n = 800 but tight enough to catch a silent
  # regression in the forwarding logic.
  shape_hat <- fit$parfit$par1$point
  rate_hat <- fit$parfit$par2$point
  expect_equal(shape_hat, shape_true, tolerance = 0.25)
  expect_equal(rate_hat, rate_true, tolerance = 0.25)

  # Exercise the summary.parfitml dprimary display path (the fn_name
  # namespace-walker that labels non-default primary event densities).
  out_full <- utils::capture.output(summary(fit))
  expect_true(any(grepl("Primary event dist", out_full, fixed = TRUE)))
  expect_true(any(grepl("dexpgrowth", out_full, fixed = TRUE)))
  expect_true(any(grepl("r = ", out_full, fixed = TRUE)))
  out_compact <- utils::capture.output(summary(fit, type = "compact"))
  expect_true(any(grepl("Primary event dist", out_compact, fixed = TRUE)))
  expect_true(any(grepl("dexpgrowth", out_compact, fixed = TRUE)))
})

test_that("dexpgrowth vs dunif give materially different fits", {
  skip_if_no_primarycensored()
  set.seed(20260413L)
  # Same simulation as above but compared against a uniform-primary fit.
  # A user who forgets `dprimary` on exponentially-skewed data would hit
  # this branch silently; the test confirms the mismatched fit drifts
  # away from the matched fit by a meaningful margin and therefore the
  # dprimary argument is doing real work.
  n <- 800L
  pwidth <- 3
  swidth <- 1
  r_true <- 1
  shape_true <- 10
  rate_true <- 1
  x1l <- sort(stats::runif(n, 0, 40))
  delays <- vapply(seq_len(n), function(i) {
    primarycensored::rprimarycensored(
      n = 1L,
      rdist = stats::rgamma,
      pwindow = pwidth,
      swindow = swidth,
      rprimary = primarycensored::rexpgrowth,
      rprimary_args = list(r = r_true),
      shape = shape_true,
      rate = rate_true
    )
  }, numeric(1))
  x <- data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )

  fit_matched <- parfitml(
    x = x, family = "gamma", Bboot = 3L, pgbar = FALSE,
    dprimary = primarycensored::dexpgrowth,
    dprimary_args = list(r = r_true)
  )
  fit_unif <- parfitml(
    x = x, family = "gamma", Bboot = 3L, pgbar = FALSE
  )
  # At least one of (shape, rate) must drift by more than 5% between
  # the two fits. With r = 1 on pwindow = 3 the primary density is
  # strongly weighted toward the late end of the window and on this
  # seed shape drifts by roughly 17%, rate by roughly 10%.
  drift_shape <- abs(
    fit_matched$parfit$par1$point - fit_unif$parfit$par1$point
  ) / fit_matched$parfit$par1$point
  drift_rate <- abs(
    fit_matched$parfit$par2$point - fit_unif$parfit$par2$point
  ) / fit_matched$parfit$par2$point
  expect_gt(max(drift_shape, drift_rate), 0.05)
})

test_that("mc branch errors on non-default dprimary", {
  skip_if_no_primarycensored()
  x <- make_double_data()
  expect_error(
    kerlikelihood(
      x = x, family = "gamma", likapprox = "mc",
      dprimary = primarycensored::dexpgrowth,
      dprimary_args = list(r = 0.2)
    ),
    "likapprox"
  )
  # Default uniform on mc still works.
  expect_silent(
    kerlikelihood(x = x, family = "gamma", likapprox = "mc")
  )
})

test_that("single-interval data errors on non-default dprimary", {
  skip_if_no_primarycensored()
  x <- data.frame(xl = c(0.5, 1.0, 1.5), xr = c(1.5, 2.0, 3.0))
  expect_error(
    kerlikelihood(
      x = x, family = "gamma",
      dprimary = primarycensored::dexpgrowth,
      dprimary_args = list(r = 0.2)
    ),
    "doubly interval-censored"
  )
  expect_error(
    parfitml(
      x = x, family = "gamma", Bboot = 2L, pgbar = FALSE,
      dprimary = primarycensored::dexpgrowth,
      dprimary_args = list(r = 0.2)
    ),
    "doubly interval-censored"
  )
  # Default uniform on nc == 2 still works.
  expect_silent(
    kerlikelihood(x = x, family = "gamma")
  )
})
