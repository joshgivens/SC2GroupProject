#' A function to do discrete time sampling from your distribution.
#' 
#' Does step wise discrete sampling for a given number of particles and a given discrete time transitiona kernel
#'
#' @param sampler The function that will be used to sample at each timestep. 
#' Takes in one argument which is the current value of x .
#' @param x_0 The starting point. If `NULL` then `init_sampler` must be specified.
#' @param timesteps The number of steps to take
#' @param n The number of particles to sample
#' @param init_sampler A function that will be used to produce the initial sample.
#' Take in one argument which is the number of points to sample. Only to be specified if `x_0=NULL`
#'
#' @examples 
#' #Try taking true samples from Ornstein-Uhlenbeck process at 100 timepoints of step size 1
#' true_sampler <- function(x_0){
#'     ornuhl_samp(1,x_0,t=1,sigma=1,alpha=1)
#' }
#' 
#' out_true <- discrete_sampler(true_sampler,x_0=0.01,timesteps = 100,n=100)
#' # Now plot the first 10 particles
#' plot(out_true[1,],type="l",ylim=c(-0.4,0.4))
#' for (i in 1:10){
#'   lines(out_true[i,],col=i)
#' }
#' @return A matrix containing the sampled particles. Each row is a particle and each column is a time-step.
#' @export
#'
discrete_sampler <- function(sampler,x_0=NULL,timesteps=100,n=1000,
                             init_sampler=NULL){
  Xmat <- matrix(NA,nrow=n,ncol=timesteps+1)
  if(!is.null(x_0)){
    Xmat[,1] <- x_0
  }
  else if (!is.null(init_sampler)){
    Xmat[,1] <- init_sampler(n)
  }
  for (t in 1:timesteps+1){
    for(i in 1:n){
      Xmat[i,t] <- sampler(Xmat[i,t-1])
    }
  }
  return(Xmat)
}



