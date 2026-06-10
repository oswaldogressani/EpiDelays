#' Summarize parametric model fit with maximum likelihood
#'
#' @description This routine can be used to summarize a parametric model fit
#' obtained with maximum likelihood.
#'
#' @param object An object of class \code{parfitml}.
#' @param ndigits Number of digits to print in output.
#' @param type Type of summary output. Either \code{"full"} or \code{"compact"}.
#' @param ... Further arguments to be passed.
#'
#' @return A summary output.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export

summary.parfitml <- function(object, ndigits = 3, type = "full", ...) {
  if(!inherits(object, "parfitml"))
    stop("object must be a parfitml object")
  if(!(type %in% c("full", "compact"))) {
    stop("type must either be 'full' or 'compact'")
  }
  dfpar <- cbind(name = unlist(object[names(object$parfit)]),
                 round(do.call(rbind, lapply(object$parfit, as.data.frame)), 
                       ndigits))
  dfdel <- round(do.call(rbind, lapply(object$delayfit, as.data.frame)), 
                 ndigits)
  cimethod <- object$cimethod
  if (cimethod == "npboot") {
    cidesc <- "Nonparametric bootstrap"
  } else if (cimethod == "pboot") {
    cidesc <- "Parametric bootstrap"
  } else if (cimethod == "sbnorm") {
    cidesc <- "Simulated (asympt. normality)"
  }
  if(type == "full") { #--- Print output (full)
    cat("------------------------------------------------------------- \n")
    cat("Parametric model fit (maximum likelihood) \n")
    cat("------------------------------------------------------------- \n")
    cat("Parametric family            : ", object$fname,"\n")
    cat("Number of parameters         : ", object$npars, "\n")
    cat("Censoring type               : ", object$censtype, "\n")
    cat("Sample size                  : ", object$n, "\n")
    cat("Routine time (seconds)       : ", object$elapsed, "\n")
    cat("MLE convergence              : ", object$mleconv, "\n")
    if(cimethod == "npboot" || cimethod == "pboot"){
      cat("Bootstrap sample size        : ", object$Bboot, "\n")
      cat("Bootstrap discarded          : ", object$bootdiscard, "\n")
    } else if(cimethod == "sbnorm"){
      cat("No. of draws to compute ci   : ", object$ns, "\n")
    }
    cat("AIC                          : ", object$aic, "\n")
    cat("BIC                          : ", object$bic, "\n")
    cat("Confidence interval (ci)     : ", cidesc, "\n")
    cat("------------------------------------------------------------- \n")
    cat("Parameter estimation:          \n")
    print(dfpar)
    cat("------------------------------------------------------------- \n")
    cat("Estimated features:            \n")
    print(dfdel)
    cat("------------------------------------------------------------- \n")
    cat("point: point estimate; se: standard error \n")
    cat("ci: confidence interval \n")
    cat("------------------------------------------------------------- \n")
  } else { #--- Print output (compact)
    cat("------------------------------------------------------------- \n")
    cat("Parametric model fit (maximum likelihood) \n")
    cat("------------------------------------------------------------- \n")
    cat("Parametric family    : ", object$fname,"\n")
    cat("Number of parameters : ", object$npars, "\n")
    cat("AIC                  : ", object$aic, "\n")
    cat("BIC                  : ", object$bic, "\n")
    cat("------------------------------------------------------------- \n")
    cat("Parameter estimation:          \n")
    print(dfpar)
    cat("------------------------------------------------------------- \n")
    cat("Estimated features:            \n")
    print(dfdel[1:3,])
    cat("------------------------------------------------------------- \n")
    cat("point: point estimate; se: standard error \n")
    cat("ci: confidence interval \n")
    cat("------------------------------------------------------------- \n")
  }
}
