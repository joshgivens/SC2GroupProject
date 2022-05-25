#' Brownian Motion Step Likelihood
#' 
#' This gives the exact likelihood of a step in Brownian Motion conditional on the starting point
#'
#' @param x The proposed value of x
#' @param x_0 the current value of x
#' @param t the time between steps
#' @param sigma the variance parameter
#' @param log Logical indicating whether the log-likelihood should be returned
#'
#' @return A single number giving the likelihood of the step
#' @export
#'
brown_lik <- function(x,x_0,t,sigma,log=F){
  sd <- sqrt(t)*sigma
  return(dnorm(x,mean=x_0,sd=sd,log=log))
}

#' Ornstein-Uhlenbeck Step Likelihood
#' 
#' This gives the exact likelihood of a step in Brownian Motion conditional on the starting point
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
#' This gives the exact likelihood of a step in Sinusoidal Brownian Motion conditional on the starting point
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





