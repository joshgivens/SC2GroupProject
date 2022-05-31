#' Brownian Motion Step Likelihood
#' 
#' This gives the exact likelihood of a step in Brownian Motion conditional on the starting point
#'
#' @param x The proposed value of x
#' @param x_0 the current value of x
#' @param t the time between steps
#' @param sigma the variance parameter for the Brownian Motion
#' @param log Logical indicating whether the log-likelihood should be returned
#'
#' @return A single number giving the likelihood of the step
#' @export
#'
brown_lik <- function(x,x_0,t,sigma,log=FALSE){
  sd <- sqrt(t)*sigma
  return(dnorm(x,mean=x_0,sd=sd,log=log))
}

#' Ornstein-Uhlenbeck Step Likelihood
#' 
#' This gives the exact likelihood of a step in Brownian Motion conditional on the starting point.
#' 
#' This function calculates the likelihood of a step in an Ornstein-Uhlenbeck Process. This process is of the form
#' \eqn{\mathrm{d}X_t=\exp{-\alpha X_t}\mathrm{d}t+\sigma\mathrm{d}Wt}. Which has transition kernel after time t given by
#' \deqn{p_t(x,x_0)=\sqrt{\frac{\alpha}{\sigma(1-\exp(-2\alpha t))}}\exp\{\frac{-\alpha(x-x_0\exp(-\alpha t))^2}{2\sigma(1-\exp(-2\alpha t))}\}}
#'
#' @param x The proposed value of x
#' @param x_0 the current value of x
#' @param t the time between steps
#' @param alpha The value of the parameter alpha in the Ornstein-Uhlenbeck process (see Details) 
#' @param sigma The value of the parameter sigma in the Ornstein-Uhlenbeck process (see Details) 
#' @param log Logical indicating whether the log-likelihood should be returned
#'
#' @return A single number giving the likelihood of the step
#' @export
#'
ornuhl_lik <- function(x,x_0,t,alpha,sigma,log=F){
  sd=sqrt(sigma*(1-exp(-2*alpha*t))/alpha)
  mean=x_0*exp(-alpha*t)
  return(dnorm(x,mean=mean,sd=sd,log = log))
}

#' Sinusoidal Step likelihood
#' 
#' This gives the approximate likelihood of a step in Sinusoidal Brownian Motion conditional on the starting point
#'
#' The true diffusion process we are attempting to obtain likelihoods from is \eqn{\mathrm{d}X_t=\sin(X_t-\omega)\mathrm{d}t+\mathrm{d}W_t}. 
#' We approximate this at time `t` with transition kernel: 
#' \deqn{p_t(x,x_0)=\sqrt{\frac{1}{2\pi t^2}}\exp\{\frac{-(x-\sin(x_0-\omega)t+x_0)^2}{2t^2}\}}
#' 
#' @param x The proposed value of x
#' @param x_0 the current value of x
#' @param t the time between steps
#' @param omega The value of the parameter omega in the Sinusoidal Brownian Motion (see Details) 
#' @param log Logical indicating whether the log-likelihood should be returned
#'
#' @return A single number giving the likelihood of the step
#' @export
#'
sin_lik <- function(x,x_0,t,omega,log=F){
  mean <- sin(x_0-omega)*t+x_0
  return(dnorm(x,mean=mean,sd=t,log=log))
}





