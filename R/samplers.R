#' Brownian Motion Step Sampler
#'
#' Generates a sample after fixed time from Brownian motion given starting point.
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
#'
#' This function generates samples from a step in an Ornstein-Uhlenbeck Process after time `t` given starting point `x_0. This process is of the form
#' \eqn{\mathrm{d}X_t=\exp{-\alpha X_t}\mathrm{d}t+\sigma\mathrm{d}Wt}. Which has transition kernel after time t given by
#' \deqn{p_t(x,x_0)=\sqrt{\frac{\alpha}{2\pi\sigma(1-\exp(-2\alpha t))}}\exp\{\frac{-\alpha(x-x_0\exp(-\alpha t))^2}{2\sigma(1-\exp(-2\alpha t))}\}}
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
#' This function samples an approximate step from a sinusoidal Brownian Motion.
#' 
#' The true diffusion process we are attempting to sample from is \eqn{\mathrm{d}X_t=\sin(X_t-\omega)\mathrm{d}t+\mathrm{d}W_t}. 
#' We approximate this at time `t` by sampling with transition kernel: 
#' \deqn{p_t(x,x_0)=\sqrt{\frac{1}{2\pi t^2}}\exp\{\frac{-(x-\sin(x_0-\omega)t+x_0)^2}{2t^2}\}}
#'
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
  mean <- sin(x_0-omega)*t+x_0
  rnorm(n,mean=mean,sd=t)
}