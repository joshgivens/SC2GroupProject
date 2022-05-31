#' SMC sampler
#' 
#' For a given process, number of time-steps, and number of particles, 
#' performs Sequential Monte-Carlo (SMC) sampling. Gives back the sampled points as well as their weights.
#' 
#' 
#' For \eqn{q_n} the proposal transition (evaluated with \code{proposal_lik[[n]]} and sampled from by \code{proposal_samples[[n]]})
#' ,\eqn{p_n} the true transition (evaluated with \code{true_lik}) and \eqn{\gamma_n} the marginals (evaluated with \code{g_lik}), performs the following algorithm
#' \itemize{
#'    \item \bold{Sample} \eqn{x_0^{(i)}\sim q_0} for \eqn{i = 1,\ldots,N} and let \eqn{T} be a threshold
#'    \item \bold{Evaluate} \eqn{w_0(x_0^{(i)})=\gamma_0(x_0^{(i)})/q_0(x_0^{(i)})}
#'    \item \bold{Normalise} \eqn{W_0^{(i)}=w_0^{(i)}/\sum_{j=1}^N w_0^{(j)}}
#'    \item  For \eqn{n\geq 1}, perform the following steps for all \eqn{i\in\{1,\dotsc,N\}}
#'    \itemize{
#'       \item \bold{Sample} \eqn{a_{n-1}^{(i)}} w.p. \eqn{\mathbb{P}(a_{n-1}^{(i)}=k)=W_{n-1}^{(k)}} independently 
#'       \item \bold{Draw} \eqn{x_n^{(i)}\sim q_n(x_{n-1}^{(a_{n-1}^{(i)})}, \cdot)}
#'       \item \bold{Evaluate} \eqn{w^{(i)}_n(x^{(i)}_{n-1},x^{(i)}_n)=\frac{\gamma_n(x^{(i)}_n)p_{n}(x^{(i)}_{n-1},x^{(i)}_{n})}{\gamma_{n-1}(x^{(i)}_{n-1})q_n(x^{(i)}_{n-1},x^{(i)}_n)}}
#'       \item \bold{Normalise} \eqn{W_n^{(i)}=W_{n-1}^{(i)}w^{(i)}_n(x^{(i)}_{n-1},x^{(i)}_n)/\sum_{j=1}^N W_{n-1}^{(j)}w^{(j)}_n(x^{(j)}_{n-1},x^{(j)}_n)}
#'    }
#'}
#'   
#' @param init_sample sampler function for the initial distribution. 
#' This function takes 1 argument which is the number of samples to take
#' @param proposal_sample sampler function for conditional proposal
#' @param proposal_lik likelihood for conditional proposal
#' @param true_lik true likelihood (up to a constant)
#' @param g_lik likelihood of observed y given x
#' @param y A vector containing our observed samples
#' @param n The number of particles to simulate
#' @param timesteps The number of time-steps to model for each particle
#'
#' @return A list with matrices containing various information from each observation, each entry gives an individual observation with 
#' each column being a separate time-point and reach row being a separate  particle
#' \item{Xmat - }{A matrix containing the simulated samples.}
#' \item{Wmat - }{A matrix containing the weights associated with each sample}
#' 
#' @examples 
#' Final_T=1
#' timesteps=100
#' step=Final_T/timesteps
#' 
#' init_sample <- function(n){
#'   brown_samp(n,1,step,sigma=1)
#' }
#' 
#' sampler <- function(x_0){
#'   brown_samp(1,x_0,step,sigma=1)
#' }
#' 
#' true_sampler <- function(x_0){
#'   ornuhl_samp(1,x_0,step,alpha=1,sigma=1)
#' }
#' 
#' prop_lik <- function(x,x_0){
#'   brown_lik(x,x_0,step,sigma=1)
#' }
#' 
#' true_lik <- function(x,x_0){
#'   ornuhl_lik(x,x_0,step,alpha = 1,sigma=1)
#' }
#' 
#' #### Generate y using exponential with rate parameter x^2 ###
#' x_0 <- 1
#' hidden_x <- rep(NA,timesteps)
#' hidden_x[1] <- true_sampler(x_0)
#' # Generate hidden x sample
#' for(i in 2:100){
#'   hidden_x[i] <- true_sampler(hidden_x[i-1])
#' }
#' 
#' y=rexp(timesteps,rate=hidden_x^2)
#' 
#' g_lik <- function(y,x){
#'   dexp(y,rate=x^2)
#' }
#' 
#' out <- Basic_SMC(init_sample=init_sample, 
#'                  proposal_sample=sampler, 
#'                  proposal_lik = prop_lik,
#'                  true_lik=true_lik,
#'                  g_lik=g_lik,
#'                  y=y,
#'                  n=1000,
#'                  timesteps=timesteps)
#' @export
#'
Basic_SMC <- function(init_sample,proposal_sample, proposal_lik, true_lik, g_lik=NULL, y=NULL, n=1000, timesteps=100){
  # Set up matrix for samples
  Xmat <- matrix(NA,n,timesteps)
  # Set up matrix for weights
  Wmat <- matrix(NA,n,timesteps)
  w_new <- rep(NA,n)
  # 1st sample points from initial distribution
  x_0=init_sample(n)
  
  # DO first sample
  for(i in 1:n){
    Xmat[i,1] <- proposal_sample(x_0[i])
    Wmat[i,1] <- true_lik(Xmat[i,1],x_0[i])*g_lik(y[1],Xmat[i,1])/
      proposal_lik(Xmat[i,1],x_0[i])
  }
  Wmat[,1]=Wmat[,1]/sum(Wmat[,1])
  
  
  for (t in 2:timesteps){
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