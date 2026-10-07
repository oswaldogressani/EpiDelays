#' Computes summary statistics
#'
#' @param x  A list of entries for which the statistics are required.
#' @param pe A vector of point estimates.
#' @param boot Bootstraped statistics.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @keywords internal

kerstats_boot <- function(x, pe, boot) {
  o <- mapply(function(l, point, se, ci90l, ci90r, ci95l, ci95r) {
    c(l, list(point = point, se = se, ci90l = ci90l, ci90r = ci90r,
              ci95l = ci95l, ci95r = ci95r))},
    x, pe, apply(boot, 2, "sd"),
    apply(boot, 2, stats::quantile, prob = 0.05),
    apply(boot, 2, stats::quantile, prob = 0.95),
    apply(boot, 2, stats::quantile, prob = 0.025),
    apply(boot, 2, stats::quantile, prob = 0.975),
    SIMPLIFY = FALSE)
  return(o)
}

