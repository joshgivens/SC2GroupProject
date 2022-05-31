#' Normalising Constant Calculator
#' 
#' This function Calculates the normalising constant of a distribution from SMC samples from said distribution.
#'
#' Where \eqn{w_t^{(j)}} represents the un-normalised weight from particle \eqn{j} and time \eqn{t}, 
#' this function calculates the normalising constant of the associated distribution as:
#' \deqn{\prod_{t=1}^T\frac{1}{N}\sum_{j=1}^N w_t^{(j)}.}
#' (where \eqn{T} is our final time-point and we have \eqn{N} total particles.)
#' 
#' @param weights A matrix containing all the weights from the SMC process. 
#' Each row represent a particle and each column represents a time-point.
#'
#' @return A single value which is the estimate of the normalising constant. 
#' @export
#'
basic_norm <- function(weights){
  prod(colMeans(weights))
}

#' Harmonic Mean Normalising Constant
#'
#' This function calculates the normalising constant using the harmonic mean. 
#' This should only be used when the proposal sampler is identical to the density.
#' 
#' If \eqn{x_t^{i},~w_t^{(j)}} represent the sample and unnormalised weight from 
#' particle \eqn{j} at time \eqn{t} the this functio calculates the normalising constant as:
#' \deqn{\frac{\sum_{j=1}^N |g(x_T^{(j)})|^{1-\varphi}\prod_{t=1}^T w_t^{(j)}}
#' {\sum_{j=1}^N |g(x_T^{(j)})|^{-\varphi}\prod_{t=1}^T w_t^{(j)}}.}
#' 
#' This function is specifically for use with samples from \code{diff_SMC} or \code{diff_SMC_sin} 
#' only in the case where we use the same transition kernel as the true transition kernel.
#' 
#' @param g A vector containing the final timepoint for each particle transformed by 
#' the function which you are trying to estimate the expectation of.
#' @param weights A matrix containing all the weights from the SMC process. 
#' Each row represent a particle and each column represents a time-point.
#' 
#' @param varphi The constant between 0 and 1 used in the SMC procedure
#'
#' @return A single value giving the estimated normalising constant.
#' @export
#'
harmonic_norm <- function(g,weights,varphi){
  weight_prod <- apply(weights,1,prod)
  (sum(abs(g)^(1-varphi)*weight_prod))/
    (sum(abs(g)^(-varphi)*weight_prod))
}