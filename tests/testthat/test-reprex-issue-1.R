# Runnable version of the smoke-test reprex posted to
# https://github.com/oswaldogressani/EpiDelays/issues/1
#
# Simulates doubly-interval-censored data via
# primarycensored::rprimarycensored() and confirms that kerlikelihood()'s
# ni branch agrees with a single vectorised call to
# primarycensored::dprimarycensored() on the resulting frame. Pinned so
# future refactors cannot silently break the comment-side reprex.

pc_loglik <- function(v, x, pdist, par_map) {
  sum(do.call(
    primarycensored::dprimarycensored,
    c(
      list(
        x = x$x2l - x$x1l,
        pdist = pdist,
        pwindow = 1,
        swindow = 1,
        log = TRUE
      ),
      par_map(v)
    )
  ))
}

test_that("reprex (issue #1): gamma kerlikelihood == dprimarycensored", {
  skip_if_no_primarycensored()

  set.seed(1L)
  n <- 100L
  delays <- primarycensored::rprimarycensored(
    n, rdist = stats::rgamma, pwindow = 1, swindow = 1,
    shape = 3, rate = 1
  )
  x1l <- sample(0:14, n, replace = TRUE)
  x <- data.frame(
    x1l = x1l,
    x1r = x1l + 1L,
    x2l = x1l + delays,
    x2r = x1l + delays + 1L
  )
  # kerlikelihood's gamma domain check rejects rows where x2l - x1r < 0.
  x <- x[x$x2l - x$x1r >= 0, ]

  v <- log(c(3, 1))  # shape = 3, rate = 1

  m <- kerlikelihood(x, family = "gamma", likapprox = "ni")
  expect_equal(
    pc_loglik(v, x, stats::pgamma,
              function(v) list(shape = exp(v[1]), rate = exp(v[2]))),
    m$loglik(v, x),
    tolerance = 1e-8
  )
})

test_that("reprex (issue #1): gaussian kerlikelihood == dprimarycensored", {
  # The gaussian family forwards stats::pnorm into dprimarycensored, so the
  # production loglik must equal a single vectorised pc_loglik() call against
  # the same CDF. kerdomain_check permits negative-supporting families on
  # any data, so no row filter is required here.
  skip_if_no_primarycensored()

  set.seed(2L)
  n <- 100L
  delays <- primarycensored::rprimarycensored(
    n, rdist = stats::rnorm, pwindow = 1, swindow = 1,
    mean = 2.5, sd = 1
  )
  x1l <- sample(0:14, n, replace = TRUE)
  x <- data.frame(
    x1l = x1l,
    x1r = x1l + 1L,
    x2l = x1l + delays,
    x2r = x1l + delays + 1L
  )

  v <- c(2.5, log(1))  # mean = 2.5, sd = 1

  m <- kerlikelihood(x, family = "gaussian", likapprox = "ni")
  expect_equal(
    pc_loglik(v, x, stats::pnorm,
              function(v) list(mean = v[1], sd = exp(v[2]))),
    m$loglik(v, x),
    tolerance = 1e-8
  )
})
