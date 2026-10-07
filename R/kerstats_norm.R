#' Computes summary statistics
#'
#' @param x  A list of entries for which the statistics are required.
#' @param pe A vector of point estimates.
#' @param se Estimated standard error.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerstats_norm <- function(x, pe, se) {
  z095 <- stats::qnorm(p = 0.95)
  z0975 <- stats::qnorm(p = 0.975)
  o <- mapply(function(l, point, se, ci90l, ci90r, ci95l, ci95r) {
    c(l, list(point = point, se = se, ci90l = ci90l, ci90r = ci90r,
              ci95l = ci95l, ci95r = ci95r))},
    x, pe, se, pe - z095 * se, pe + z095 * se,
    pe - z0975 * se, pe + z0975 * se,
    SIMPLIFY = FALSE)
  return(o)
}


