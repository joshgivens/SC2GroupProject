\name{SC2GroupProject-package}
\alias{SC2GroupProject-package}
\alias{SC2GroupProject}
\docType{package}
\title{
  Performs SMC algorithms.
}
\description{
This package implements various SMC algorithms specifically for use in diffusion processes. Alongside this there are samplers and transition probability function for various common Feynman-Kac processes. 
}
\details{
This package should main be used to implement standard SMC approaches for diffusion processes or exact SMC approaches. The basic SMC algorithm can be used in the function \code{Basic_SMC}.

To use the SMC process given in \insertCite{doucet}{SC2GroupProject} one can use the function \code{diff_SMC}.

Finally a specific C++ implementation of \code{diff_SMC} for the case given in Example 1 of \insertCite{doucet}{SC2GroupProject} can be found in the function \code{diff_SMS_sin}. This version provides 30 times speed up over the R version.
}
\author{
Josh Givens, Ettore Fincato.

Maintainer: Josh Givens <josh.givens@bristol.ac.uk>
}
\references{
\insertRef{doucet}{SC2GroupProject}
}
\keyword{ package }

\examples{
  \dontrun{
    # Generate 10,000 samples and at 20 time-steps up to time 1
    out <- diff_SMC_sin(x_0=5, t=1/20, n=10000, timesteps=20, varphi=0.9,TRUE)  
  }
}
