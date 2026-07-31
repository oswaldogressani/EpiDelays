#' Summarize nonparametric model fit
#'
#' @description This routine can be used to summarize a nonparametric model fit.
#'
#' @param object An object of class \code{nonparfit}.
#' @param ndigits Number of digits to print in output.
#' @param ... Further arguments to be passed.
#'
#' @return A summary output.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export
summary.nonparfit <- function(object, ndigits = 3, ...) {
  if(!inherits(object, "nonparfit"))
    stop("object must be a nonparfit object")
  dfdel <- round(do.call(rbind, lapply(object$delayfit, as.data.frame)), ndigits)
  cat("---------------------------------------------------- \n")
  cat("Nonparametric fit (method of Gressani and Hens 2025) \n")
  cat("---------------------------------------------------- \n")
  cat("Censoring type         : ", object$censtype, "\n")
  cat("Sample size            : ", object$n, "\n")
  cat("Routine time (seconds) : ", object$elapsed, "\n")
  cat("Bootstrap sample size  : ", object$Bboot, "\n")
  cat("---------------------------------------------------- \n")
  cat("Estimated features:            \n")
  print(dfdel)
  cat("---------------------------------------------------- \n")
  cat("point: point estimate; se: standard error \n")
  cat("ci: confidence interval \n")
  cat("---------------------------------------------------- \n")
}
