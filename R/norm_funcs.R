#' Normalising Constant Calculator
#' 
#' This function Calculates the normalising constant of a distribution from SMC samples from said distribution.
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