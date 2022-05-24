brown_lik <- function(x,x_0,t,sigma,log=F){
  sd <- sqrt(t)*sigma
  return(dnorm(x,mean=x_0,sd=sd,log=log))
}

ornuhl_lik <- function(x,x_0,t,alpha,sigma,log=F){
  sd=sqrt(sigma*(1-exp(-2*alpha*t))/alpha)
  mean=x_0*exp(-alpha*t)
  return(dnorm(x,mean=mean,sd=sd,log = log))
}

sin_lik <- function(x,x_0,t,omega,log=F){
  mean <- sin(x_0-omega)*t+x_0
  return(dnorm(x,mean=mean,sd=t,log=log))
}





