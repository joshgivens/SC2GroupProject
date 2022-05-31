test_that("diff_SMC produces correct dimensions", {
  #Set our final timepoints
  Final_T=1
  #Set number of timesteps we will use
  timesteps=20
  #Get the size of each timestep
  step=Final_T/timesteps
  
  #Create sampling function
  sampler <- function(x_0){
    sin_samp(1,x_0,step,pi)
  }
  
  # Create Likelihood Function
  true_lik <- function(x,x_0){
    sin_lik(x,x_0,step,pi)
  }
  
  # Create g function
  g <- function(x){
    zeta <- 100
    omega2 <- 1e-30
    diff=sqrt(zeta)-omega2
    if(x>-diff & x<diff){
      return(
        exp(1/(zeta-x^2))
      )
    }
    else{
      return(0)
    }
  }
  #Simulate 1,000 particles
  out <- diff_SMC(g=g,p_delta=true_lik,M_T_samplers = rep(list(sampler),timesteps),
                  M_T_lik = rep(list(true_lik),timesteps), n = 1000, timesteps = timesteps,
                  varphi = 0.9,x_0=1,resample = T)
  
  expect_equal(dim(out$Wmat),dim(out$Xmat))
  expect_equal(dim(out$Wmat),dim(out$gmat))
  expect_equal(dim(out$Wmat),c(1000,timesteps))
})
