# Helper to skip tests cleanly when primarycensored is not installed.
skip_if_no_primarycensored <- function() {
  testthat::skip_if_not_installed("primarycensored")
}
