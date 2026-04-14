# Tests for left/right truncation support in kerlikelihood() and parfitml().
#
# Semantics: L and D truncate the underlying delay distribution. The per-row
# log-likelihood picks up a `-log(F_cens(D) - F_cens(L))` correction term.
# When L = 0 and D = Inf the defaults must reproduce the existing code path
# to machine precision. Covers design doc notes/truncation-design.md sections
# 1-5.

make_double_data <- function(seed = 1L, n = 5L) {
  set.seed(seed)
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

trunc_family_cases <- list(
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
  )
)

# Default behaviour: specifying L = 0, D = Inf must be numerically identical
# to the existing (no-arg) call. Guards against accidental drift in the
# default code path when L/D are added.
for (fam in names(trunc_family_cases)) {
  local({
    family <- fam
    case <- trunc_family_cases[[family]]
    test_that(sprintf(
      "%s default L=0 D=Inf matches existing ni loglik exactly", family
    ), {
      skip_if_no_primarycensored()
      x <- make_double_data()
      m_default <- kerlikelihood(x = x, family = family, likapprox = "ni")
      m_trunc <- kerlikelihood(
        x = x, family = family, likapprox = "ni", L = 0, D = Inf
      )
      expect_identical(
        m_default$loglik(case$v, x),
        m_trunc$loglik(case$v, x)
      )
    })
  })
}

test_that("nc==4 truncation subtracts per-row log(F_cens(D) - F_cens(L))", {
  skip_if_no_primarycensored()
  x_all <- make_double_data()
  v <- c(1.5, log(0.8))
  L <- 0.1
  D <- 10
  # Keep rows whose (x2l - x1l, x2r - x1l) fit inside [L, D].
  keep <- (x_all$x2l - x_all$x1l) >= L & (x_all$x2r - x_all$x1l) <= D
  x <- x_all[keep, ]
  expect_gt(nrow(x), 0L)

  # Oracle: call dprimarycensored() directly row-by-row with L, D against
  # the same raw stats::pnorm CDF that kerlikelihood() forwards.
  expected <- sum(vapply(seq_len(nrow(x)), function(i) {
    primarycensored::dprimarycensored(
      x = x$x2l[i] - x$x1l[i],
      pdist = stats::pnorm,
      pwindow = x$x1r[i] - x$x1l[i],
      swindow = x$x2r[i] - x$x2l[i],
      L = L,
      D = D,
      log = TRUE,
      mean = v[1],
      sd = exp(v[2])
    )
  }, numeric(1)))

  m <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "ni", L = L, D = D
  )
  expect_equal(m$loglik(v, x), expected, tolerance = 1e-8)
})

test_that("nc==4 ni drops left-straddle rows with a warning", {
  skip_if_no_primarycensored()
  # Mix one fully-inside row with a left-straddle row whose observed lower
  # x2l - x1l = 0.6 sits below L = 1. dprimarycensored cannot evaluate a
  # window that crosses L, so the straddle row must be dropped (with a
  # warning advising the caller) before the inner loop runs. The fit then
  # proceeds on the surviving row alone.
  L <- 1
  D <- 10
  x <- data.frame(
    x1l = c(0, 0),
    x1r = c(0.5, 0.5),
    x2l = c(2, 0.6),
    x2r = c(3, 1.5)
  )
  v <- log(c(2, 0.5))
  expect_warning(
    kerlikelihood(x = x, family = "gamma", likapprox = "ni", L = L, D = D),
    regexp = "straddle"
  )
  m <- suppressWarnings(
    kerlikelihood(x = x, family = "gamma", likapprox = "ni", L = L, D = D)
  )
  expect_identical(nrow(m$x), 1L)
  ker_value <- m$loglik(v, m$x)
  expect_true(is.finite(ker_value))

  expected <- primarycensored::dprimarycensored(
    x = m$x$x2l - m$x$x1l,
    pdist = stats::pgamma,
    pwindow = m$x$x1r - m$x$x1l,
    swindow = m$x$x2r - m$x$x2l,
    L = L, D = D, log = TRUE,
    shape = 2, rate = 0.5
  )
  expect_equal(ker_value, sum(expected), tolerance = 1e-10)
})

test_that("nc==4 ni drops right-straddle rows with a warning", {
  skip_if_no_primarycensored()
  # Symmetric case: x2r - x1l = 11 sits above D = 10. The straddle row is
  # dropped with a warning and the fit continues on the in-bounds row.
  L <- 0.5
  D <- 10
  x <- data.frame(
    x1l = c(0, 0),
    x1r = c(0.5, 0.5),
    x2l = c(2, 5),
    x2r = c(3, 11)
  )
  v <- log(c(2, 0.5))
  expect_warning(
    kerlikelihood(x = x, family = "gamma", likapprox = "ni", L = L, D = D),
    regexp = "straddle"
  )
  m <- suppressWarnings(
    kerlikelihood(x = x, family = "gamma", likapprox = "ni", L = L, D = D)
  )
  expect_identical(nrow(m$x), 1L)
  ker_value <- m$loglik(v, m$x)
  expect_true(is.finite(ker_value))

  expected <- primarycensored::dprimarycensored(
    x = m$x$x2l - m$x$x1l,
    pdist = stats::pgamma,
    pwindow = m$x$x1r - m$x$x1l,
    swindow = m$x$x2r - m$x$x2l,
    L = L, D = D, log = TRUE,
    shape = 2, rate = 0.5
  )
  expect_equal(ker_value, sum(expected), tolerance = 1e-10)
})

test_that("nc==2 truncation matches closed-form truncated CDF", {
  skip_if_no_primarycensored()
  x <- data.frame(
    xl = c(1.0, 2.0, 3.0, 4.0),
    xr = c(2.0, 3.0, 4.0, 5.0)
  )
  v <- log(c(2, 0.5))
  shape <- 2
  rate <- 0.5
  L <- 0.5
  D <- 6

  Fr <- pmin(stats::pgamma(x$xr, shape, rate), stats::pgamma(D, shape, rate))
  Fl <- pmax(stats::pgamma(x$xl, shape, rate), stats::pgamma(L, shape, rate))
  FD <- stats::pgamma(D, shape, rate)
  FL <- stats::pgamma(L, shape, rate)
  expected <- sum(log(Fr - Fl) - log(FD - FL))

  m <- kerlikelihood(
    x = x, family = "gamma", likapprox = "ni", L = L, D = D
  )
  expect_equal(m$loglik(v, x), expected, tolerance = 1e-12)
})

test_that("nc==2 default L=0 D=Inf matches untruncated single-interval", {
  skip_if_no_primarycensored()
  x <- data.frame(
    xl = c(0.5, 1.5, 2.5),
    xr = c(1.5, 2.5, 3.5)
  )
  v <- log(c(2, 0.5))
  m_default <- kerlikelihood(x = x, family = "gamma", likapprox = "ni")
  m_trunc <- kerlikelihood(
    x = x, family = "gamma", likapprox = "ni", L = 0, D = Inf
  )
  expect_identical(m_default$loglik(v, x), m_trunc$loglik(v, x))
})

test_that("parfitml recovers gamma parameters under truncation", {
  skip_if_no_primarycensored()
  set.seed(20260413L)
  n_raw <- 10000L
  shape_true <- 3
  rate_true <- 1
  # Left-truncation at 4 cuts most of the gamma(3, 1) mass, so the naive
  # fit that ignores truncation sees a sample with mean far above 3 and
  # converges to biased parameters. Keeps the upper bound generous so that
  # enough rows survive for a stable truncated fit.
  L <- 4
  D <- 12
  # Draw integer-day doubly censored gamma delays.
  x1l <- sort(stats::runif(n_raw, 0, 5))
  pwidth <- rep(1, n_raw)
  swidth <- rep(1, n_raw)
  delays <- vapply(seq_len(n_raw), function(i) {
    primarycensored::rprimarycensored(
      n = 1L,
      rdist = stats::rgamma,
      pwindow = pwidth[i],
      swindow = swidth[i],
      shape = shape_true,
      rate = rate_true
    )
  }, numeric(1))
  x_all <- data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )
  # Filter to rows whose observed interval sits inside [L, D].
  keep <- (x_all$x2l - x_all$x1l) >= L & (x_all$x2r - x_all$x1l) <= D
  x <- x_all[keep, ]
  expect_gt(nrow(x), 200L)

  fit <- parfitml(
    x = x, family = "gamma", Bboot = 3L, pgbar = FALSE, L = L, D = D
  )
  expect_true(fit$mleconv)
  shape_hat <- fit$parfit$par1$point
  rate_hat <- fit$parfit$par2$point
  # Recovery under a tight truncation window is harder than an untruncated
  # fit; a 15% tolerance is enough to separate success from the badly
  # biased naive fit below.
  expect_lt(abs(shape_hat - shape_true) / shape_true, 0.15)
  expect_lt(abs(rate_hat - rate_true) / rate_true, 0.15)

  # Unadjusted fit on the same truncated sample should be materially biased:
  # the estimator sees a mean pinned near the midpoint of [L, D].
  fit_naive <- parfitml(
    x = x, family = "gamma", Bboot = 3L, pgbar = FALSE
  )
  mean_true <- shape_true / rate_true
  mean_naive <- fit_naive$parfit$par1$point / fit_naive$parfit$par2$point
  expect_gt(abs(mean_naive - mean_true) / mean_true, 0.10)
})

test_that("parfitml single-interval recovers gamma under truncation", {
  skip_if_no_primarycensored()
  set.seed(20260414L)
  n_raw <- 3000L
  shape_true <- 3
  rate_true <- 1
  L <- 1
  D <- 10
  delays <- stats::rgamma(n_raw, shape = shape_true, rate = rate_true)
  keep <- delays >= L & delays <= D
  delays <- delays[keep]
  expect_gt(length(delays), 200L)
  x <- data.frame(
    xl = pmax(floor(delays), L),
    xr = pmin(floor(delays) + 1, D)
  )
  x <- x[x$xl < x$xr, ]

  fit <- parfitml(
    x = x, family = "gamma", Bboot = 5L, pgbar = FALSE, L = L, D = D
  )
  expect_true(fit$mleconv)
  shape_hat <- fit$parfit$par1$point
  rate_hat <- fit$parfit$par2$point
  expect_lt(abs(shape_hat - shape_true) / shape_true, 0.15)
  expect_lt(abs(rate_hat - rate_true) / rate_true, 0.15)
})

test_that("parfitml warns and drops nc==4 rows incompatible with [L, D]", {
  skip_if_no_primarycensored()
  set.seed(20260414L)
  # A handful of in-bounds rows alongside one row whose secondary window
  # falls entirely past D. The straddle-row drop in kerlikelihood() keeps
  # the fit alive on the survivors and surfaces a warning so the caller
  # knows observations were excluded.
  n_keep <- 30L
  x_keep <- data.frame(
    x1l = sort(stats::runif(n_keep, 0, 5)),
    x1r = NA, x2l = NA, x2r = NA
  )
  x_keep$x1r <- x_keep$x1l + 1
  x_keep$x2l <- x_keep$x1l + stats::rgamma(n_keep, 3, 1)
  x_keep$x2r <- x_keep$x2l + 1
  x_keep <- x_keep[
    (x_keep$x2l - x_keep$x1l) >= 1 & (x_keep$x2r - x_keep$x1l) <= 10,
  ]
  x_bad <- data.frame(x1l = 0, x1r = 1, x2l = 12, x2r = 13)
  x <- rbind(x_keep, x_bad)
  expect_warning(
    parfitml(
      x = x, family = "gamma", Bboot = 2L, pgbar = FALSE, L = 1, D = 10
    ),
    regexp = "straddle"
  )
  fit <- suppressWarnings(
    parfitml(
      x = x, family = "gamma", Bboot = 2L, pgbar = FALSE, L = 1, D = 10
    )
  )
  expect_true(fit$mleconv)
  expect_identical(fit$n, nrow(x_keep))
})

test_that("parfitml warns and drops nc==2 rows incompatible with [L, D]", {
  skip_if_no_primarycensored()
  set.seed(20260415L)
  # In-bounds rows plus one row that falls fully below L. The fit should
  # warn, drop the offending row, and converge on the rest.
  delays <- stats::rgamma(40L, shape = 3, rate = 1)
  delays <- delays[delays >= 1 & delays <= 10]
  x_keep <- data.frame(
    xl = pmax(floor(delays), 1),
    xr = pmin(floor(delays) + 1, 10)
  )
  x_keep <- x_keep[x_keep$xl < x_keep$xr, ]
  x_bad <- data.frame(xl = 0.1, xr = 0.4)  # fully below L
  x <- rbind(x_keep, x_bad)
  expect_warning(
    parfitml(
      x = x, family = "gamma", Bboot = 2L, pgbar = FALSE, L = 1, D = 10
    ),
    regexp = "straddle"
  )
  fit <- suppressWarnings(
    parfitml(
      x = x, family = "gamma", Bboot = 2L, pgbar = FALSE, L = 1, D = 10
    )
  )
  expect_true(fit$mleconv)
  expect_identical(fit$n, nrow(x_keep))
})

test_that("boundary parity: x2r - x1l == D yields finite loglik", {
  skip_if_no_primarycensored()
  x <- data.frame(
    x1l = 0,
    x1r = 1,
    x2l = 4,
    x2r = 5
  )
  D <- 5
  L <- 0
  m <- kerlikelihood(
    x = x, family = "gamma", likapprox = "ni", L = L, D = D
  )
  val <- m$loglik(log(c(2, 0.5)), x)
  expect_true(is.finite(val))
})

test_that("parfitml errors when L is negative", {
  skip_if_no_primarycensored()
  x <- data.frame(
    xl = c(1, 2),
    xr = c(2, 3)
  )
  expect_error(
    parfitml(
      x = x, family = "gaussian", Bboot = 2L, pgbar = FALSE, L = -1
    ),
    regexp = "non-negative"
  )
})

test_that("parfitml errors when L >= D", {
  skip_if_no_primarycensored()
  x <- data.frame(
    xl = c(1, 2),
    xr = c(2, 3)
  )
  expect_error(
    parfitml(
      x = x, family = "gaussian", Bboot = 2L, pgbar = FALSE, L = 5, D = 5
    ),
    regexp = "less than"
  )
})

test_that("parfitml smoke test: all five families with L/D set", {
  skip_if_no_primarycensored()
  set.seed(42L)
  n <- 120L
  x1l <- sort(stats::runif(n, 0, 2))
  pwidth <- rep(0.5, n)
  swidth <- rep(0.5, n)
  # Gaussian-drawn delays with a generous mean keep the skewnorm MoM seed
  # (which relies on the third sample moment) well-defined.
  delays <- vapply(seq_len(n), function(i) {
    primarycensored::rprimarycensored(
      n = 1L,
      rdist = stats::rnorm,
      pwindow = pwidth[i],
      swindow = swidth[i],
      mean = 6,
      sd = 1.5
    )
  }, numeric(1))
  x_all <- data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )
  L <- 0.5
  D <- 20
  keep <- (x_all$x2l - x_all$x1l) >= L & (x_all$x2r - x_all$x1l) <= D
  x <- x_all[keep, ]
  families <- c("gaussian", "gamma", "lognormal", "weibull", "skewnorm")
  for (fam in families) {
    fit <- parfitml(
      x = x, family = fam, Bboot = 3L, pgbar = FALSE, L = L, D = D
    )
    expect_true(fit$mleconv, info = fam)
    expect_true(all(is.finite(fit$parfit$par1$se)), info = fam)
  }
})

test_that("summary.parfitml surfaces L and D when non-default", {
  skip_if_no_primarycensored()
  set.seed(5L)
  n <- 30L
  x1l <- sort(stats::runif(n, 0, 2))
  pwidth <- rep(0.5, n)
  swidth <- rep(0.5, n)
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
  x_all <- data.frame(
    x1l = x1l,
    x1r = x1l + pwidth,
    x2l = x1l + delays,
    x2r = x1l + delays + swidth
  )
  L <- 0.25
  D <- 15
  keep <- (x_all$x2l - x_all$x1l) >= L & (x_all$x2r - x_all$x1l) <= D
  x <- x_all[keep, ]
  fit <- parfitml(
    x = x, family = "gamma", Bboot = 3L, pgbar = FALSE, L = L, D = D
  )
  out <- utils::capture.output(summary(fit))
  expect_true(any(grepl("Left truncation", out, fixed = TRUE)))
  expect_true(any(grepl("Right truncation", out, fixed = TRUE)))
})

test_that("kerlikelihood mc branch applies truncation correction", {
  skip_if_no_primarycensored()
  x_all <- make_double_data()
  v <- c(1.5, log(0.8))
  L <- 0.1
  D <- 10
  keep <- (x_all$x2l - x_all$x1l) >= L & (x_all$x2r - x_all$x1l) <= D
  x <- x_all[keep, ]
  # Seed controls the primary draws used inside the closure. The mc and ni
  # corrections share the same numerator, so checking that the mc result
  # sits close to the ni result is a sanity check on the correction term.
  m_ni <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "ni", L = L, D = D
  )
  m_mc <- kerlikelihood(
    x = x, family = "gaussian", likapprox = "mc", L = L, D = D
  )
  set.seed(20260413L)
  v_mc <- m_mc$loglik(v, x)
  v_ni <- m_ni$loglik(v, x)
  expect_equal(v_mc, v_ni, tolerance = 5e-2)
})
