# Regression tests for pskewnorm's Owen's T saturation.
#
# In extreme parameter regions visited mid-optim, pskewnorm can return values
# a hair outside [0, 1] or slightly non-monotone on sorted inputs due to
# floating-point noise in Owen's T. primarycensored::check_pdist rejects
# such CDFs, which previously aborted parfitml() on the skewnorm family for
# unlucky small-sample benchmarks. pskewnorm() clamps to [0, 1] and lifts
# rounding-noise dips with cummax over the order of x; parfitmom()'s
# skewnorm seed uses a signed cube root so negative sample skew does not
# seed optim with NaN.

test_that("kerlikelihood skewnorm ni branch tolerates saturation (narrow)", {
  skip_if_no_primarycensored()
  set.seed(1L)
  n <- 10L
  x <- data.frame(
    x1l = 0:(n - 1L),
    x1r = 1:n,
    x2l = (1:n) + stats::runif(n, 0, 0.2),
    x2r = (1:n) + 1 + stats::runif(n, 0, 0.2)
  )
  m <- kerlikelihood(x = x, family = "skewnorm", likapprox = "ni")
  # Narrow, heavily-skewed regime: tiny scale with a large positive slant
  # drives pskewnorm into the Owen's T saturation region.
  val <- m$loglik(c(0, log(0.1), 5), x)
  expect_true(is.finite(val))
})

test_that("kerlikelihood skewnorm ni branch tolerates saturation (shifted)", {
  skip_if_no_primarycensored()
  set.seed(2L)
  n <- 10L
  x <- data.frame(
    x1l = 0:(n - 1L),
    x1r = 1:n,
    x2l = (1:n) + stats::runif(n, 0, 0.2),
    x2r = (1:n) + 1 + stats::runif(n, 0, 0.2)
  )
  m <- kerlikelihood(x = x, family = "skewnorm", likapprox = "ni")
  # Far-from-data location with a large negative slant. The resulting CDF
  # is evaluated at quantiles deep in the tail where Owen's T noise is
  # worst; without the clamp-and-cummax wrapper check_pdist aborts.
  val <- m$loglik(c(50, log(10), -5), x)
  expect_true(is.finite(val))
})

test_that("kerlikelihood skewnorm ni branch tolerates saturation under L/D", {
  skip_if_no_primarycensored()
  set.seed(3L)
  n <- 12L
  # Every row's secondary window sits safely inside [L, D] so the
  # straddle-row drop does not mask the saturation path. Far-from-data
  # location combined with a large negative slant places the data deep in
  # the skewnorm's left tail where Owen's T returns rounding-noise CDF
  # values; the non-default truncation bounds force dprimarycensored to
  # additionally evaluate F_cens(L) and F_cens(D) on the same saturating
  # CDF. Without pskewnorm's internal clamp-and-cummax check_pdist would
  # abort the call.
  x <- data.frame(
    x1l = 0:(n - 1L),
    x1r = 1:n,
    x2l = (1:n) + 1 + stats::runif(n, 0, 0.2),
    x2r = (1:n) + 1.5 + stats::runif(n, 0, 0.2)
  )
  L <- 0.5
  D <- 100
  m <- kerlikelihood(
    x = x, family = "skewnorm", likapprox = "ni", L = L, D = D
  )
  val <- m$loglik(c(50, log(10), -5), x)
  expect_true(is.finite(val))
})

test_that("parfitml skewnorm converges on a small negative-skew sample", {
  skip_if_no_primarycensored()
  # The original benchmark n = 50 cell DNF'd because the MoM seed was NaN
  # (negative sample skew -> rsn^(2/3) = NaN in R). parfitmom now uses a
  # signed cube root so the seed is always finite; the optimiser should
  # converge on this tiny sample without hitting check_pdist saturation.
  set.seed(50L)
  n <- 50L
  pwindow <- 1
  swindow <- 1
  delays <- primarycensored::rprimarycensored(
    n = n,
    rdist = function(n) {
      EpiDelays::qskewnorm(
        stats::runif(n), par1 = 1, par2 = 1, par3 = 2
      )
    },
    pwindow = pwindow,
    swindow = swindow
  )
  x1l <- sort(stats::runif(n, 0, 5))
  x <- data.frame(
    x1l = x1l, x1r = x1l + pwindow,
    x2l = x1l + delays, x2r = x1l + delays + swindow
  )
  fit <- parfitml(
    x = x, family = "skewnorm", Bboot = 3L, pgbar = FALSE
  )
  expect_true(fit$mleconv)
  expect_true(all(is.finite(fit$parfit$par1$point)))
})
