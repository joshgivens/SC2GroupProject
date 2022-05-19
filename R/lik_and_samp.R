ornuhl_lik <- function(x,x_0,t,alpha,sigma,log=F){
  sd=alpha/(sigma*(1-exp(-2*alpha*t)))
  mean=x_0*exp(-alpha*t)
  dnorm(x,mean=mean,sd=sd,log = log)
}

brown_lik <- function(x,x_0,t,sigma,log=F){
  sd <- sqrt(t)*sigma
  dnorm(x,mean=x_0,sd=sd,log=log)
}

brown_samp <- function(n,x_0,t,sigma){
  sd=sqrt(t)*sigma
  rnorm(n,mean=x_0,sd=sd)
}