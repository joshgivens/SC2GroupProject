#' SMC sampler
#'
#' @param init_sample sampler function for the initial distibution. This function takes 1 argument which is the nuber of samples to taje
#' @param proposal_sample 
#' @param proposal_lik 
#' @param true_lik 
#' @param g_lik
#'
#' @return
#' @export
#'
#' @examples
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
    Xmat[i,2]=proposal_sample(Xmat[i,t-1])
    Wmat[i,2] <- true_lik(x_new[i],Xmat[i,t-1])*g(y[t],x_new[i])/
      proposal_lik(x_new[i],Xmat[i,t-1])
  }
  Wmat[,2]=Wmat[,2]/sum(Wmat[,2])
  
  
  for (t in 3:timesteps){
    # resample previous iter to get new points
    x_sampled <- sample(Xmat[,j-1],size=n,replace=TRUE,prob=Wmat[,j-1])
    
    for(i in 1:n){
      #Get new sample
      Xmat[i,t]=proposal_sample(x_sampled[i])
      #Get weights associated with new sample
      Wmat[i,t] <- true_lik(Xmat[i,t],Xmat[i,t-1])*g(y[t],Xmat[i,t])/
        proposal_lik(Xmat[i,t],Xmat[i,t-1])
    }
    #Normalise weights
    Wmat[,t]=Wmat[,t]/sum(Wmat[,t])
  }
  
}