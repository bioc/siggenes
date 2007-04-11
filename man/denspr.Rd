\name{denspr}
\alias{denspr}
\title{Density Estimation}
\description{
  Estimates the density of a vector of observations by a Poisson regression
  fit to histogram counts.
}
\usage{
  denspr(x, n.interval = NULL, df = 7)
}
\arguments{
  \item{x}{a numeric vector containing the observations for which the density
     should be estimated}
  \item{n.interval}{an integer specifying the number of cells for the histogram.
     If \code{NULL}, \code{n.interval} is set to the squared root of the number
     of observations}
  \item{df}{integer specifying the degrees of freedom of the natural cubic spline
     used in the Poisson regression fit}
}

\value{
  A numeric vector of the same length as \code{x} containing the estimated density
  for each of the observations
}

\references{
  Efron, B., and Tibshirani, R. (1996). Using specially designed exponential 
  families for density estimation. \emph{Annals of Statistics}, 24, 2431--2461. 
}

\keyword{utilities}
