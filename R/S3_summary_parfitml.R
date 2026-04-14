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
                 round(do.call(rbind, lapply(object$parfit, as.data.frame)), ndigits))
  dfdel <- round(do.call(rbind, lapply(object$delayfit, as.data.frame)), ndigits)
  # Describe the primary event density only when it has been overridden; a
  # default uniform primary is the legacy behaviour and need not clutter
  # the summary block.
  dprimary_is_default <- is.null(object$dprimary) ||
    (identical(object$dprimary, stats::dunif) &&
       length(object$dprimary_args) == 0)
  if(!dprimary_is_default) {
    dprimary_label <- deparse(substitute(NULL)) # placeholder
    fn_name <- tryCatch(
      {
        env <- environment(object$dprimary)
        nm <- NULL
        if(!is.null(env)) {
          nsname <- tryCatch(
            getNamespaceName(env), error = function(e) NULL
          )
          if(!is.null(nsname)) {
            candidates <- ls(envir = env, all.names = TRUE)
            for(c in candidates) {
              if(identical(get(c, envir = env), object$dprimary)) {
                nm <- paste0(nsname, "::", c)
                break
              }
            }
          }
        }
        if(is.null(nm)) nm <- "custom"
        nm
      },
      error = function(e) "custom"
    )
    if(length(object$dprimary_args) > 0) {
      arg_bits <- paste0(
        names(object$dprimary_args), " = ",
        vapply(object$dprimary_args, function(a) {
          paste(format(a), collapse = ", ")
        }, character(1))
      )
      dprimary_label <- paste0(
        fn_name, " (", paste(arg_bits, collapse = ", "), ")"
      )
    } else {
      dprimary_label <- fn_name
    }
  }
  # Truncation bounds are only surfaced when they differ from the defaults;
  # the untruncated path continues to print its original summary block.
  has_trunc <- isTRUE(object$L > 0) ||
    (!is.null(object$D) && is.finite(object$D))
  if(type == "full") { #--- Print output (full)
    cat("---------------------------------------------------- \n")
    cat("Parametric model fit (maximum likelihood) \n")
    cat("---------------------------------------------------- \n")
    cat("Parametric family      : ", object$fname,"\n")
    cat("Number of parameters   : ", object$npars, "\n")
    cat("Censoring type         : ", object$censtype, "\n")
    cat("Sample size            : ", object$n, "\n")
    cat("Routine time (seconds) : ", object$elapsed, "\n")
    cat("MLE convergence        : ", object$mleconv, "\n")
    cat("Bootstrap sample size  : ", object$Bboot, "\n")
    cat("Bootstrap discarded    : ", object$bootdiscard, "\n")
    if(!dprimary_is_default) {
      cat("Primary event dist     : ", dprimary_label, "\n")
    }
    if(has_trunc) {
      cat("Left truncation (L)    : ", object$L, "\n")
      cat("Right truncation (D)   : ", object$D, "\n")
    }
    cat("AIC                    : ", object$aic, "\n")
    cat("BIC                    : ", object$bic, "\n")
    cat("---------------------------------------------------- \n")
    cat("Parameter estimation:          \n")
    print(dfpar)
    cat("---------------------------------------------------- \n")
    cat("Estimated features:            \n")
    print(dfdel)
    cat("---------------------------------------------------- \n")
    cat("point: point estimate; se: standard error \n")
    cat("ci: confidence interval \n")
    cat("---------------------------------------------------- \n")
  } else { #--- Print output (compact)
    cat("---------------------------------------------------- \n")
    cat("Parametric model fit (maximum likelihood) \n")
    cat("---------------------------------------------------- \n")
    cat("Parametric family    : ", object$fname,"\n")
    cat("Number of parameters : ", object$npars, "\n")
    if(!dprimary_is_default) {
      cat("Primary event dist   : ", dprimary_label, "\n")
    }
    if(has_trunc) {
      cat("Left truncation (L)  : ", object$L, "\n")
      cat("Right truncation (D) : ", object$D, "\n")
    }
    cat("AIC                  : ", object$aic, "\n")
    cat("BIC                  : ", object$bic, "\n")
    cat("---------------------------------------------------- \n")
    cat("Parameter estimation:          \n")
    print(dfpar)
    cat("---------------------------------------------------- \n")
    cat("Estimated features:            \n")
    print(dfdel[1:3,])
    cat("---------------------------------------------------- \n")
    cat("point: point estimate; se: standard error \n")
    cat("ci: confidence interval \n")
    cat("---------------------------------------------------- \n")
  }
}
