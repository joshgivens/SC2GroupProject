% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/direct_sampler.R
\name{discrete_sampler}
\alias{discrete_sampler}
\title{A function to do discrete time sampling from your distribution.}
\usage{
discrete_sampler(
  sampler,
  x_0 = NULL,
  timesteps = 100,
  n = 1000,
  init_sampler = NULL
)
}
\arguments{
\item{sampler}{The function that will be used to sample at each timestep. 
Takes in one argument which is the current value of x .}

\item{x_0}{The starting point. If `NULL` then `init_sampler` must be specified.}

\item{timesteps}{The number of steps to take}

\item{n}{The number of particles to sample}

\item{init_sampler}{A function that will be used to produce the initial sample.
Take in one argument which is the number of points to sample. Only to be specified if `x_0=NULL`}
}
\value{
A matrix containing the sampled particles. Each row is a particle and each column is a time-step.
}
\description{
Does step wise discrete sampling for a given number of particles and a given discrete time transitiona kernel
}
\examples{
#Try taking true samples from Ornstein-Uhlenbeck process at 100 timepoints of step size 1
true_sampler <- function(x_0){
    ornuhl_samp(1,x_0,t=1,sigma=1,alpha=1)
}

out_true <- discrete_sampler(true_sampler,x_0=0.01,timesteps = 100,n=100)
# Now plot the first 10 particles
plot(out_true[1,],type="l",ylim=c(-0.4,0.4))
for (i in 1:10){
  lines(out_true[i,],col=i)
}
}
