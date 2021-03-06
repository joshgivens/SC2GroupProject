% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplers.R
\name{ornuhl_samp}
\alias{ornuhl_samp}
\title{Ornstein-Uhlenbeck Step Sampler}
\usage{
ornuhl_samp(n, x_0, t, sigma, alpha)
}
\arguments{
\item{n}{Number of points to sample}

\item{x_0}{Starting point}

\item{t}{Time between starting point and sampled point}

\item{sigma}{The value of the parameter sigma in the Ornstein-Uhlenbeck process (see Details)}

\item{alpha}{The value of the parameter alpha in the Ornstein-Uhlenbeck process (see Details)}
}
\value{
A vector containing samples from the Ornstein-Uhlenbeck Process
}
\description{
This function generates samples from a step in an Ornstein-Uhlenbeck Process after time `t` given starting point `x_0. This process is of the form
\eqn{\mathrm{d}X_t=\exp{-\alpha X_t}\mathrm{d}t+\sigma\mathrm{d}Wt}. Which has transition kernel after time t given by
\deqn{p_t(x,x_0)=\sqrt{\frac{\alpha}{2\pi\sigma(1-\exp(-2\alpha t))}}\exp\{\frac{-\alpha(x-x_0\exp(-\alpha t))^2}{2\sigma(1-\exp(-2\alpha t))}\}}
}
