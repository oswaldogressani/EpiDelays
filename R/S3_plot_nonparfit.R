#' Plot the cumulative distribution function of a fitted nonparametric model
#'
#' @description Can be used to plot the cumulative distribution function of a
#' nonparametric model fitted with the \code{nonparfit} routine.
#'
#' @param x An object of class \code{nonparfit}.
#' @param xlim A numeric vector of size 2 giving the x-axis range. Default left
#' bound is given by the minimum of the observed left bound and default right
#' bound is given by the maximum of the observed right bound.
#' @param grid The number of grid points on which to evaluate the fitted
#' distribution.
#' @param legend Add a legend?
#' @param ... Further graphical parameters to be passed to the plot function.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export

plot.nonparfit <- function(x, xlim = NULL, grid = 500L, legend = TRUE,...){
  if(!inherits(x, "nonparfit"))
    stop("x must be a nonparfit object")
  xl <- x$x[, 1]
  xr <- x$x[, 2]
  xw <- xr - xl
  if(is.null(xlim)) { # Default domain
    xmin <- min(xl)
    xmax <- max(xr)
  } else{ # User-defined domain
    xmin <- xlim[1]
    xmax <- xlim[2]
  }
  xx <- seq(xmin, xmax, length = grid)
  ninv <- 1 / x$n
  Fhat <- function(t) ninv * sum((t - xl) / xw * (t >= xl & t <= xr) + (t > xr))
  dots <- list(...)
  if(!is.null(dots$col)){
    legcol <- dots$col
  } else{
    legcol <- "black"
  }
  if(!is.null(dots$lty)){
    leglty <- dots$lty
  } else{
    leglty <- 1
  }
  if(!is.null(dots$lwd)){
    leglwd <- dots$lwd
  } else{
    leglwd <- 1
  }
  graphics::plot(x = xx, y = sapply(xx, Fhat),
       ylab = "Cumulative distribution function",...)
  if(isTRUE(legend)){
    graphics::legend("topleft", col = legcol, lty = leglty, lwd = leglwd,
           "Nonparametric", bty = "n")
  }
}

