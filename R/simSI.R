#' Simulation of artificial serial interval data
#'
#' @description
#' This routine simulates artificial serial interval data with serial
#' interval windows having a width of at least two days. It assumes that the
#' target serial interval (SI) distribution is Gaussian. The width of the
#' serial interval windows can be controlled by specifying a discrete censoring
#' distribution.
#'
#' @usage simSI(muS = 3, sdS = 2, maxcoarse = 3, probcoarse = c(0.8,0.15,0.05), n = 10)
#'
#' @param muS The mean of the target serial interval distribution.
#' @param sdS The standard deviation of the serial interval distribution.
#' @param maxcoarse The largest degree of coarseness allowed (in days).
#' @param probcoarse The (discrete) censoring probability distribution (must sum to one).
#' @param n The sample size (number of transmission pairs).
#'
#' @details The \code{simSI} routine outputs a data frame containing the
#'  generated serial interval values (which are in practice unobserved), the left
#'  and right bound of the serial windows (observed) and the width of the
#'  serial interval windows.
#'
#' @return A list with the following components:
#' \itemize{
#'  \item s: The true serial interval value.
#'  \item sl: Left bound of the serial interval window.
#'  \item sr: Right bound of the serial interval window.
#'  \item sw: Width of the serial interval window.
#' }
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @examples
#' # Simulate serial interval windows for n=15 transmission pairs
#' sigen <- simSI(n = 15)
#'
#' @export

simSI <- function(muS = 3, sdS = 2, maxcoarse = 3, probcoarse = c(0.8,0.15,0.05), n = 10){
 s <- stats::rnorm(n = n, mean = muS, sd = sdS)
 coarseness <- seq_len(maxcoarse)
 c <- sample(coarseness, size = n, replace = TRUE, prob = probcoarse)
 u <- stats::runif(n, min = 0, max = 1)
 sl <- floor((s - u * c))
 sr <- ceiling((s + (1 - u) * c))
 sw <- sr - sl
 x <- data.frame(s = s, sl = sl, sr = sr, sw = sw)
 return(x)
}






