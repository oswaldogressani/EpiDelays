#' Computes summary statistics
#'
#' @param slist  A list of entries for which the statistics are required.
#' @param pestim A vector of point estimates.
#' @param method The method with which confidence intervals are computed.
#' @param se     A vector of standard errors.
#' @param boot   Bootstraped statistics.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#'
#' @keywords internal

kerstats <- function(slist, pestim, method = c("norm", "boot"), se = NULL, 
                     boot = NULL) {
  citype <- match.arg(method)
  if(citype == "norm"){
    if(is.null(se))
      stop("Standard errors missing.")
    z095 <- stats::qnorm(p = 0.95)
    z0975 <- stats::qnorm(p = 0.975)
    o <- mapply(function(l, point, se, ci90l, ci90r, ci95l, ci95r) {
      c(l, list(point = point, se = se, ci90l = ci90l, ci90r = ci90r,
                ci95l = ci95l, ci95r = ci95r))},
      slist, pestim, se, pestim - z095 * se, pestim + z095 * se,
      pestim - z0975 * se, pestim + z0975 * se, SIMPLIFY = FALSE)
  }else if(citype == "boot"){
    if(is.null(boot))
      stop("Bootstrap sample is missing.")
    o <- mapply(function(l, point, se, ci90l, ci90r, ci95l, ci95r) {
      c(l, list(point = point, se = se, ci90l = ci90l, ci90r = ci90r,
                ci95l = ci95l, ci95r = ci95r))},
      slist, pestim, apply(boot, 2, "sd"),
      apply(boot, 2, stats::quantile, prob = 0.05),
      apply(boot, 2, stats::quantile, prob = 0.95),
      apply(boot, 2, stats::quantile, prob = 0.025),
      apply(boot, 2, stats::quantile, prob = 0.975),
      SIMPLIFY = FALSE)
  }
  return(o)
}
