#' Title
#'
#' @param g Function which we are taking expectation over
#' @param p_delta true transition kernel between discrete time steps
#' @param M_T_samplers list of samplers for T approximating transition kernels
#' @param M_T_lik list of likelihood for T approximating kernels
#' @param n Number of particles to simulate
#' @param timesteps Number of time-steps to simulate
#' @param varphi varphi function used in approximation
#' @param x_0 starting value
#'
#' @return A list containing the simulated particles and their associated weights
#' @export
#'
diff_SMC <- function(g,p_delta,M_T_samplers,M_T_lik,n=1000,timesteps=100,varphi=1,x_0){
  # Set up matrix for samples
  Xmat <- matrix(NA,n,timesteps)
  # Set up matrix for weights
  Wmat <- matrix(NA,n,timesteps)
  gmat <- matrix(NA,n,timesteps)
  gmat_transf <- matrix(NA,n,timesteps)
  w_new <- rep(NA,n)
  # 1st sample points from initial distribution
  
  #Calculate init g-value  
  for (i in 1:n){
    Xmat[i,1] <- M_T_samplers[[1]](x_0)
    gmat[i,1] <- g(Xmat[i,1])
    gmat_transf[i,1] <- abs(gmat[i,1])^varphi 
    Wmat[i,1] <- gmat_transf[i,1]*p_delta(Xmat[i,1],x_0)/
      (M_T_lik[[1]](Xmat[i,1],x_0))
  }
  
  for (t in 2:timesteps){
    # resample previous iter to get new points
    x_sampled <- sample(Xmat[,t-1],size=n,replace=TRUE,
                        prob=Wmat[,t-1]/sum(Wmat[,t-1]))
    
    for(i in 1:n){
      #Get new sample
      Xmat[i,t]=M_T_samplers[[t]](x_sampled[i])
      #Get weights associated with new sample
      gmat[i,t] <- g(Xmat[i,t])
      gmat_transf[i,t] <- abs(gmat[i,t])^varphi
      Wmat[i,t] <- gmat_transf[i,t]*p_delta(Xmat[i,t],Xmat[i,t-1])/
        (gmat_transf[i,t-1]*M_T_lik[[t]](Xmat[i,t],Xmat[i,t-1]))
    }
  }
  return(list(Xmat=Xmat, Wmat=Wmat, gmat=gmat))
}
