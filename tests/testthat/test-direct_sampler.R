test_that("Test discrete_sampler gives correct dimensions", {
  true_sampler <- function(x_0){
    brown_samp(1,x_0,t=1,sigma=2)
  }
  out_true <- discrete_sampler(true_sampler,x_0=0.01,timesteps = 100,n=10)
  expect_equal(dim(out_true),c(10,101))
})
