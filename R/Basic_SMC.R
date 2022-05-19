#' SMC sampler
#'
#' @param init_sample sampler function for the initial distribution. 
#' This function takes 1 argument which is the number of samples to take
#' @param proposal_sample sampler function for conditional proposal
#' @param proposal_lik likelihood for conditional proposal
#' @param true_lik true likelihood (up to a constant)
#' @param g_lik likelihood of observed y given x
#'
#' @return A list containing the simulated particles and their associated weights
#' @export
#'
Basic_SMC <- function(init_sample,proposal_sample, proposal_lik, true_lik, g_lik=NULL, y=NULL, n=1000, timesteps=100){
  # Set up matrix for samples
  Xmat <- matrix(NA,n,timesteps)
  # Set up matrix for weights
  Wmat <- matrix(NA,n,timesteps)
  w_new <- rep(NA,n)
  # 1st sample points from initial distribution
  Xmat[,1]=init_sample(n)
  
  # DO first sample
  for(i in 1:n){
    Xmat[i,2]=proposal_sample(Xmat[i,1])
    Wmat[i,2] <- true_lik(Xmat[i,2],Xmat[i,1])*g_lik(y[t],Xmat[i,2])/
      proposal_lik(Xmat[i,2],Xmat[i,1])
  }
  Wmat[,2]=Wmat[,2]/sum(Wmat[,2])
  
  
  for (t in 3:timesteps){
    # resample previous iter to get new points
    x_sampled <- sample(Xmat[,t-1],size=n,replace=TRUE,prob=Wmat[,t-1])
    
    for(i in 1:n){
      #Get new sample
      Xmat[i,t]=proposal_sample(x_sampled[i])
      #Get weights associated with new sample
      Wmat[i,t] <- true_lik(Xmat[i,t],Xmat[i,t-1])*g_lik(y[t],Xmat[i,t])/
        proposal_lik(Xmat[i,t],Xmat[i,t-1])
    }
    #Normalise weights
    Wmat[,t]=Wmat[,t]/sum(Wmat[,t])
  }
  return(list(Xmat=Xmat, Wmat=Wmat))
}