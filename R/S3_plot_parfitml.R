#' Plot the density or cumulative distribution function of a fitted parametric model
#'
#' @description Can be used to plot the probability density function
#' or the cumulative distribution function of a parametric model fitted with
#' the \code{parfitml} routine.
#'
#' @param x An object of class \code{parfitml}.
#' @param xlim A numeric vector of size 2 giving the x-axis range. Default left
#' bound is given by the minimum of the observed left bound and default right
#' bound is given by the maximum of the observed right bound.
#' @param grid The number of grid points on which to evaluate the fitted
#' distribution.
#' @param target The fitted target function to be plotted, `pdf` is the default
#' and gives the probability density function. Specifying `cdf` plots the fitted
#' cumulative distribution function.
#' @param legend Add a legend?
#' @param ... Further graphical parameters to be passed to the plot function.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export

plot.parfitml <- function(x, xlim = NULL, grid = 500L,
                          target = c("pdf", "cdf"), legend = TRUE, ...) {
  if(!inherits(x, "parfitml"))
    stop("x must be a parfitml object")
  tartype <- match.arg(target)
  if(is.null(xlim)) { # Default domain
    xmin <- x$xmin
    xmax <- x$xmax
  } else{ # User-defined domain
    xmin <- xlim[1]
    xmax <- xlim[2]
  }
  xx <- seq(xmin, xmax, length = grid)
  parvals <- unname(as.list(sapply(x$parfit, `[[`, 1)))
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
  if(tartype == "pdf"){
    tarname <- x$dname
    eval <- c(list(x = xx), parvals)
    legpos <- "topright"
  } else{
    tarname <- x$pname
    eval <- c(list(q = xx), parvals)
    legpos <- "topleft"
  }
  if(x$fname == "skewnorm"){
    graphics::plot(xx, do.call(get(paste0(tarname)), eval),...)
  } else{
    graphics::plot(xx, do.call(get(paste0(tarname), asNamespace("stats")), eval),...)
  }
  if(isTRUE(legend)){
    graphics::legend(legpos, col = legcol, lty = leglty, lwd = leglwd,
           x$fname, bty = "n")
  }
}
