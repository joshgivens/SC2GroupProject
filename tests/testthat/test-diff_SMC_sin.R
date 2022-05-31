test_that("diff_SMC_sin produces correct dimensions", {
  timesteps=100
  step=1/timesteps
  out <- diff_SMC_sin(x_0=1,t=step,n=10000,timesteps=timesteps,varphi=0.9,TRUE)
  #Simulate 10,000 particles
  
  expect_equal(dim(out$Wmat),dim(out$Xmat))
  expect_equal(dim(out$Wmat),dim(out$gmat))
  expect_equal(dim(out$Wmat),c(10000,timesteps))
})
