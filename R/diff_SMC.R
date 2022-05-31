#' SMC Sampler 2
#' 
#' Samples using SMC as in the framework laid out in \insertCite{doucet}{SC2GroupProject}.
#'
#' @param g Function which we are taking expectation over
#' @param p_delta true transition kernel between discrete time steps
#' @param M_T_samplers list of samplers for T approximating transition kernels
#' @param M_T_lik list of likelihood for T approximating kernels
#' @param n Number of particles to simulate
#' @param timesteps Number of time-steps to simulate
#' @param varphi varphi function used in approximation
#' @param x_0 starting value
#' @param resample Logical indicating whether to resample or not
#'
#' @return A list with matrices containing various information from each observation, each entry gives an individual observation with 
#' each column being a separate time-point and reach row being a separate  particle
#' \item{Xmat - }{A matrix containing the simulated samples.}
#' \item{Wmat - }{A matrix containing the weights associated with each sample}
#' \item{gmat - }{A matrix containing the value of the function g evaulated at each sample.}
#'
#'
#' @references
#' \insertRef{doucet}{SC2GroupProject}
#'
#' @examples
#' # Set-up to process described in section 5.2 of ...
#' # Set our final timepoints
#' Final_T=1
#' #Set number of timesteps we will use
#' timesteps=20
#' #Get the size of each timestep
#' step=Final_T/timesteps
#' 
#' #Create sampling function
#' sampler <- function(x_0){
#'   sin_samp(1,x_0,step,pi)
#' }
#' 
#' # Create Likelihood Function
#' true_lik <- function(x,x_0){
#'   sin_lik(x,x_0,step,pi)
#' }
#' 
#' # Create g function
#' g <- function(x){
#'   zeta <- 100
#'   omega2 <- 1e-30
#'   diff=sqrt(zeta)-omega2
#'   if(x>-diff & x<diff){
#'     return(
#'       exp(1/(zeta-x^2))
#'     )
#'   }
#'   else{
#'     return(0)
#'   }
#' }
#' 
#' # Now run our SMC sampler
#' out <- diff_SMC(g=g,p_delta=true_lik,M_T_samplers = rep(list(sampler),timesteps),
#'     M_T_lik = rep(list(true_lik),timesteps), n = 1000, timesteps = timesteps,
#'     varphi = 0.9,x_0=5,resample = TRUE)
#' 
#' @export
#'
diff_SMC <- function(g,p_delta,M_T_samplers,M_T_lik,n=1000,timesteps=100,varphi=1,x_0,resample=T){
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
    if (resample){
      #resample previous iter to get new points
      x_sampled <- sample(Xmat[,t-1],size=n,replace=TRUE,
                          prob=Wmat[,t-1]/sum(Wmat[,t-1]))
      Wmat[,t] <- 1
    }
    else{
      x_sampled <- Xmat[,t-1] 
      Wmat[,t] <- Wmat[,t-1]
    }
    
    for(i in 1:n){
      #Get new sample
      Xmat[i,t]=M_T_samplers[[t]](x_sampled[i])
      #Get weights associated with new sample
      gmat[i,t] <- g(Xmat[i,t])
      gmat_transf[i,t] <- abs(gmat[i,t])^varphi
      Wmat[i,t] <- Wmat[i,t]*gmat_transf[i,t]*p_delta(Xmat[i,t],Xmat[i,t-1])/
        (gmat_transf[i,t-1]*M_T_lik[[t]](Xmat[i,t],Xmat[i,t-1]))
    }
  }
  return(list(Xmat=Xmat, Wmat=Wmat, gmat=gmat))
}
