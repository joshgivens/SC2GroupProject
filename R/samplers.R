#' Brownian Motion Step Sampler
#'
#' @param n Number of points to sample
#' @param x_0 Starting point
#' @param t Time between starting point and sampled point
#' @param sigma The value of the parameter Sigma in the Brownian Motion (see Details)
#'
#' @return A vector containing samples from the Brownian Motion
#' @export
#'
brown_samp <- function(n,x_0,t,sigma){
  sd=sqrt(t)*sigma
  rnorm(n,mean=x_0,sd=sd)
}

#' Ornstein-Uhlenbeck Step Sampler
#'
#' @param n Number of points to sample
#' @param x_0 Starting point
#' @param t Time between starting point and sampled point
#' @param alpha The value of the parameter alpha in the Ornstein-Uhlenbeck process (see Details) 
#' @param sigma The value of the parameter sigma in the Ornstein-Uhlenbeck process (see Details) 
#'
#' @return A vector containing samples from the Ornstein-Uhlenbeck Process
#' @export
#'
ornuhl_samp <- function(n,x_0,t,sigma,alpha){
  sd=sqrt(sigma*(1-exp(-2*alpha*t))/alpha)
  mean=x_0*exp(-alpha*t)
  rnorm(n,mean=mean,sd=sd)
}

#' Sinusoidal Brownian Motion Step Sampler
#'
#' @param n Number of points to sample
#' @param x_0 Starting point
#' @param t Time between starting point and sampled point
#' @param omega The value of the parameter omega in the Sinusoidal Brownian Motion (see Details) 
#'
#' @return A vector containing samples from the Sinusoidal Brownian Motion
#' @export
#'
sin_samp <- function(n,x_0,t,omega){
  mean <- sin(x_0-omega)+x_0*t
  rnorm(n,mean=mean,sd=t)
}