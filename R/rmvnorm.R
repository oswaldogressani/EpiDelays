#' Simulation from a multivariate Gaussian distribution
#'
#' @description
#' Simulates random draws from a multivariate normal distribution with mean 
#' vector \code{mean} and covariance matrix \code{sigma}. The number of draws
#' is specfied by \code{n} (a positive integer). The routine first tries to 
#' compute the Cholesky factorization of the covariance matrix by assuming that
#' the latter is a symmetric positive-definite matrix. If this approach fails,
#' the routine relies on an eigendecomposition of the covariance matrix as a 
#' fallback option.
#' 
#' @param n A positive integer indicating the number of draws.
#' @param mean The mean vector of the multivariate Gaussian.
#' @param sigma The covariance matrix of the multivariate Gaussian distribution.
#'
#' @return A list containing the simulated draws and the method used for 
#' covariance matrix decomposition.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}.
#'
#' @references Ripley, B. D. (2009). Stochastic simulation. \emph{John Wiley & Sons}.
#' @references Rubinstein, R. Y. and Kroese, D. P. (2016). Simulation and
#' the Monte Carlo method. \emph{John Wiley & Sons}.
#'
#' @export

rmvnorm <- function(n = 1, mean, sigma){
  if(!is.numeric(n) || length(n)!=1 || n<=0)
    stop("Number of draws n must be a positive integer.")
  if(!isSymmetric(sigma))
    stop("Matrix sigma must be symmetric.")
  d <- length(mean)
  if(d != nrow(sigma))
    stop("The dimension of mean and covar do not match.")
  NAmat <- matrix(NA, nrow = d, ncol = d)
  # Try Cholesky decomposition and compute lower triangular matrix
  method <- "cholesky"
  S <- tryCatch({t(chol(sigma))}, error = function(e){NAmat})
  if(anyNA(S)){ # Fallback to eigendecomposition
    eigdec <- eigen(sigma, symmetric = TRUE)
    eigval <- eigdec$values
    if (any(eigval < 0)) {
      stop("Matrix sigma must be positive semi-definite.")
    }
    S <- eigdec$vectors %*% diag(sqrt(eigval), d)
    method <- "eigendecomposition"
  }
  Z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  o <- list(sim = sweep(Z %*% t(S), 2, mean, "+"), method = method)
  return(o)
}





