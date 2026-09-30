# Upstream API (logliki, J, xmin/xmax, sbnorm, plot) with our backend.

api_cases <- list(
  gaussian = list(v = c(1.5, log(0.8))),
  gamma = list(v = log(c(2, 0.5))),
  lognormal = list(v = c(0.5, log(0.4))),
  weibull = list(v = log(c(2, 2))),
  skewnorm = list(v = c(1, log(1), 2))
)

# Interleaved primary windows, to check results keep row order.
mixed_window_data <- function() {
  data.frame(
    x1l = c(0, 1, 2, 3, 4, 5),
    x1r = c(1, 3, 3, 5, 5, 7),
    x2l = c(3, 4.5, 5, 6.2, 7, 8.4),
    x2r = c(4, 5.5, 6, 7.2, 8, 9.4)
  )
}

for (fam in names(api_cases)) {
  local({
    family <- fam
    v <- api_cases[[family]]$v

    test_that(sprintf("%s logliki is pointwise and sums to loglik", family), {
      skip_if_no_primarycensored()
      x <- mixed_window_data()
      m <- kerlikelihood(
        x = x, family = family, L = 0.5, D = 20,
        dprimary = primarycensored::dexpgrowth,
        dprimary_args = list(r = 0.2)
      )
      li <- m$logliki(v, x)
      expect_length(li, nrow(x))
      expect_identical(sum(li), m$loglik(v, x))
      by_row <- vapply(seq_len(nrow(x)), function(i) {
        m$loglik(v, x[i, , drop = FALSE])
      }, numeric(1))
      expect_equal(li, by_row, tolerance = 1e-10)
    })

    test_that(sprintf("%s J matches a numerical Jacobian", family), {
      skip_if_no_primarycensored()
      m <- kerlikelihood(x = mixed_window_data(), family = family)
      h <- 1e-6
      num <- vapply(seq_along(v), function(j) {
        e <- replace(numeric(length(v)), j, h)
        (as.numeric(m$originscale(v + e)) -
           as.numeric(m$originscale(v - e))) / (2 * h)
      }, numeric(length(v)))
      expect_equal(m$J(v), num, tolerance = 1e-6)
    })
  })
}

test_that("kerlikelihood returns delay range after truncation drop", {
  skip_if_no_primarycensored()
  x <- mixed_window_data()
  m <- kerlikelihood(x = x, family = "gamma")
  expect_identical(m$xmin, min(x$x2l - x$x1r))
  expect_identical(m$xmax, max(x$x2r - x$x1l))
  x2 <- data.frame(xl = c(1, 2, 3), xr = c(2, 3, 12))
  m2 <- suppressWarnings(kerlikelihood(x = x2, family = "gamma", D = 10))
  expect_identical(c(m2$xmin, m2$xmax), c(1, 3))
})

test_that("parfitml sbnorm CIs compose with truncation and dprimary", {
  skip_if_no_primarycensored()
  set.seed(11L)
  n <- 200L
  delays <- primarycensored::rprimarycensored(
    n, rdist = stats::rgamma, pwindow = 1, swindow = 1, D = 15,
    rprimary = primarycensored::rexpgrowth,
    rprimary_args = list(r = 0.2), shape = 3, rate = 1
  )
  x1l <- sample(0:10, n, replace = TRUE)
  x <- data.frame(
    x1l = x1l, x1r = x1l + 1, x2l = x1l + delays, x2r = x1l + delays + 1
  )
  x <- x[x$x2l - x$x1r >= 0, ]
  utils::capture.output({
    fit <- parfitml(
      x = x, family = "gamma", ci = "sbnorm", ns = 50, D = 15,
      dprimary = primarycensored::dexpgrowth, dprimary_args = list(r = 0.2)
    )
  })
  expect_true(fit$mleconv)
  expect_identical(fit$cimethod, "sbnorm")
  expect_identical(fit$D, 15)
  expect_true(all(is.finite(c(fit$parfit$par1$se, fit$parfit$par2$se))))
  expect_lt(abs(fit$parfit$par1$point - 3), 1)
  s <- utils::capture.output(summary(fit))
  expect_true(any(grepl("Simulated", s, fixed = TRUE)))
  expect_true(any(grepl("Right truncation", s, fixed = TRUE)))
  expect_true(any(grepl("dexpgrowth", s, fixed = TRUE)))
})

test_that("plot.parfitml works on a truncated doubly censored fit", {
  skip_if_no_primarycensored()
  set.seed(12L)
  n <- 60L
  delays <- primarycensored::rprimarycensored(
    n, rdist = stats::rlnorm, pwindow = 1, swindow = 1, D = 12,
    meanlog = 1, sdlog = 0.5
  )
  x1l <- sample(0:10, n, replace = TRUE)
  x <- data.frame(
    x1l = x1l, x1r = x1l + 1, x2l = x1l + delays, x2r = x1l + delays + 1
  )
  x <- x[x$x2l - x$x1r >= 0, ]
  utils::capture.output({
    fit <- parfitml(
      x = x, family = "lognormal", ci = "sbnorm", ns = 20, D = 12
    )
  })
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, target = "cdf"))
})

test_that("pskewnorm takes q and stays a valid CDF when saturated", {
  q <- c(5, -5, 0, 60, -60)
  p <- pskewnorm(q = q, par1 = 50, par2 = 10, par3 = -5)
  expect_true(all(p >= 0 & p <= 1))
  expect_false(is.unsorted(p[order(q)]))
})

test_that("logliki gives fresh results when x changes between calls", {
  skip_if_no_primarycensored()
  # The bootstrap passes resampled frames to the same closure.
  x <- mixed_window_data()
  m <- kerlikelihood(x = x, family = "gamma")
  v <- log(c(2, 0.5))
  first <- m$logliki(v, x)
  xb <- x[c(6, 1, 1, 3), ]
  expect_equal(
    m$logliki(v, xb),
    kerlikelihood(x = xb, family = "gamma")$logliki(v, xb),
    tolerance = 1e-12
  )
  expect_identical(m$logliki(v, x), first)
  expect_equal(
    first, kerlikelihood(x = x, family = "gamma")$logliki(v, x),
    tolerance = 1e-12
  )
})

test_that("logliki handles duplicated rows and parameter changes", {
  skip_if_no_primarycensored()
  x <- mixed_window_data()[c(1, 2, 1, 4, 2, 2, 6, 5, 3), ]
  m <- kerlikelihood(x = x, family = "weibull", L = 1, D = 25)
  for (v in list(log(c(2, 2)), log(c(1.2, 5)), log(c(3, 1.5)))) {
    expected <- vapply(seq_len(nrow(x)), function(i) {
      primarycensored::dprimarycensored(
        x = x$x2l[i] - x$x1l[i], pdist = stats::pweibull,
        pwindow = x$x1r[i] - x$x1l[i], swindow = x$x2r[i] - x$x2l[i],
        L = 1, D = 25, shape = exp(v[1]), scale = exp(v[2]), log = TRUE
      )
    }, numeric(1))
    expect_equal(m$logliki(v, x), expected, tolerance = 1e-12)
  }
})

oracle_cases <- list(
  gaussian = list(
    pdist = stats::pnorm, pars = function(v) list(mean = v[1], sd = exp(v[2]))
  ),
  gamma = list(
    pdist = stats::pgamma,
    pars = function(v) list(shape = exp(v[1]), rate = exp(v[2]))
  ),
  lognormal = list(
    pdist = stats::plnorm,
    pars = function(v) list(meanlog = v[1], sdlog = exp(v[2]))
  ),
  weibull = list(
    pdist = stats::pweibull,
    pars = function(v) list(shape = exp(v[1]), scale = exp(v[2]))
  ),
  skewnorm = list(
    pdist = pskewnorm,
    pars = function(v) list(par1 = v[1], par2 = exp(v[2]), par3 = v[3])
  )
)
bounds <- list(l_only = c(1, Inf), d_only = c(-Inf, 20), both = c(1, 20))
primaries <- list(
  uniform = list(d = stats::dunif, args = list()),
  expgrowth = list(d = primarycensored::dexpgrowth, args = list(r = 0.3))
)

for (fam in names(oracle_cases)) {
  for (b in names(bounds)) {
    for (pr in names(primaries)) {
      local({
        family <- fam
        oc <- oracle_cases[[family]]
        lims <- bounds[[b]]
        prim <- primaries[[pr]]
        test_that(sprintf("%s logliki matches dprimarycensored (%s, %s)",
                          family, b, pr), {
          skip_if_no_primarycensored()
          x <- mixed_window_data()
          v <- api_cases[[family]]$v
          m <- kerlikelihood(
            x = x, family = family, L = lims[1], D = lims[2],
            dprimary = prim$d, dprimary_args = prim$args
          )
          expected <- vapply(seq_len(nrow(x)), function(i) {
            do.call(primarycensored::dprimarycensored, c(list(
              x = x$x2l[i] - x$x1l[i], pdist = oc$pdist,
              pwindow = x$x1r[i] - x$x1l[i], swindow = x$x2r[i] - x$x2l[i],
              L = lims[1], D = lims[2], dprimary = prim$d,
              primary_args = prim$args, log = TRUE
            ), oc$pars(v)))
          }, numeric(1))
          expect_equal(m$logliki(v, x), expected, tolerance = 1e-10)
        })
      })
    }
  }
}
