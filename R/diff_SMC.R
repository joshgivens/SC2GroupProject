diff_SMC <- function(g,p_delta,Mn_samplers,Mn_lik,n=1000,timesteps=100,varphi=1,x_0){
  # Set up matrix for samples
  Xmat <- matrix(NA,n,timesteps)
  # Set up matrix for weights
  Wmat <- matrix(NA,n,timesteps)
  gmat <- matrix(NA,n,timesteps)
  w_new <- rep(NA,n)
  # 1st sample points from initial distribution
  
  #Calculate init g-value  
  g_0=abs(g(x_0))^varphi
  for (i in 1:n){
    Xmat[i,1] <- Mn_samplers[[1]](x_0)
    gmat[i,1] <- abs(g(Xmat[i,1]))^varphi
    Wmat[i,1] <- gmat[i,1]*p_delta(Xmat[i,1],x_0)/
      (g_0*Mn_lik[[1]](Xmat[i,1],x_0))
  }
  #renormalise weights
  Wmat[,1] <- Wmat[,1]/sum(Wmat[,1])
  
  for (t in 2:timesteps){
    # resample previous iter to get new points
    x_sampled <- sample(Xmat[,j-1],size=n,replace=TRUE,prob=Wmat[,j-1])
    
    for(i in 1:n){
      #Get new sample
      Xmat[i,j]=Mn_samplers[[t]](x_sampled[i])
      #Get weights associated with new sample
      gmat[i,t] <- abs(g(Xmat[i,t]))^varphi
      Wmat[i,t] <- gmat[i,t]*p_delta(Xmat[i,t],Xmat[i,t-1])/
        (gmat[i,t-1]*Mn_lik[[1]](Xmat[i,t],Xmat[i,t-1]))
    }
    #Normalise weights
    Wmat[,t]=Wmat[,t]/sum(Wmat[,t])
  }
  
  
}