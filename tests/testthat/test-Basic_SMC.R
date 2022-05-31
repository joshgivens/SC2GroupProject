test_that("Basic_SMC produces correct dimensions", {
  #Set up situation
  Final_T=1
  timesteps=100
  step=Final_T/timesteps
  
  init_sample <- function(n){
    brown_samp(n,1,step,sigma=1)
  }
  
  sampler <- function(x_0){
    brown_samp(1,x_0,step,sigma=1)
  }
  
  true_sampler <- function(x_0){
    ornuhl_samp(1,x_0,step,alpha=1,sigma=1)
  }
  
  prop_lik <- function(x,x_0){
    brown_lik(x,x_0,step,sigma=1)
  }
  
  true_lik <- function(x,x_0){
    ornuhl_lik(x,x_0,step,alpha = 1,sigma=1)
  }
  ## Generate observed sample
  # Generate y using exponential with shape parameter x^2
  set.seed(1)
  x_0 <- 1
  hidden_x <- rep(NA,timesteps)
  hidden_x[1] <- true_sampler(x_0)
  # Generate hidden x sample
  for(i in 2:100){
    hidden_x[i] <- true_sampler(hidden_x[i-1])
  }
  
  y=rexp(timesteps,rate=hidden_x^2)
  
  g_lik <- function(y,x){
    dexp(y,rate=x^2)
  }
  #Run SMC
  out <- Basic_SMC(init_sample=init_sample, 
                   proposal_sample=sampler, 
                   proposal_lik = prop_lik,
                   true_lik=true_lik,
                   g_lik=g_lik,
                   y=y,
                   n=1000,
                   timesteps=timesteps)
  
  expect_equal(dim(out$Wmat),dim(out$Xmat))
  expect_equal(dim(out$Wmat),c(1000,timesteps))
})
