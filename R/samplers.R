brown_samp <- function(n,x_0,t,sigma){
  sd=sqrt(t)*sigma
  rnorm(n,mean=x_0,sd=sd)
}

ornuhl_sample <- function(n,x_0,t,sigma,alpha){
  sd=sqrt(sigma*(1-exp(-2*alpha*t))/alpha)
  mean=x_0*exp(-alpha*t)
  rnorm(n,mean=mean,sd=sd)
}