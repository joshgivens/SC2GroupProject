discrete_sampler <- function(sampler,x_0=NULL,timesteps=100,n=1000,
                             init_sampler=NULL){
  Xmat <- matrix(NA,nrow=n,ncol=timesteps+1)
  if(!is.null(x_0)){
    Xmat[,1] <- x_0
  }
  else if (!is.null(init_sampler)){
    Xmat[,1] <- init_sampler(n)
  }
  for (t in 2:timesteps+1){
    for(i in 1:n){
      Xmat[i,t] <- sampler(Xmat[i,t-1])
    }
  }
  return(Xmat)
}