\name{quantiles}
\alias{quantiles}
\title{Empirical Quantiles}
\description{Computes the desired quantiles of a sample and returns them.}
\usage{quantiles(x,prob)}

\arguments{
    \item{x}{vector of data}
    \item{prob}{vector of desired probability levels. Values must be
        between 0 and 1 inclusive.}
}

\value{
    a vector containing the desired quantiles of the sample
}

\details{
Consider a sample containing \eqn{n} observations. If \eqn{\alpha n} is a non-integer, then
the \eqn{\alpha} quantile will be given by the \eqn{k}th smallest observation \eqn{x(k)}, where
\eqn{k} is the smallest integer larger than \eqn{\alpha n}. If \eqn{\alpha n} is an integer, then
the \eqn{\alpha} quantile will be computed by \eqn{0.5(x(k)+x(k+1))}.}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\keyword{univar}
