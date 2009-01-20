\name{rowWilcoxon}
\alias{rowWilcoxon}
\alias{rowRanksWilc}

\title{Rowwise Wilcoxon Rank Sum Statistics}
\description{
  Computes either the Wilcoxon Rank Sum or Signed Rank Statistics for all
  rows of a matrix simultaneously.
}
\usage{
  rowWilcoxon(X, cl, rand = NA)
}

\arguments{
  \item{X}{a matrix in which each row corresponds to a variable, and each
     column to an observation/sample.}
  \item{cl}{a numeric vector consisting of ones and zeros. The length of
     \code{cl} must be equal to the number of observations. If \code{cl}
     consists of zeros and ones, Wilcoxon Rank Sums are computed. If
     \code{cl} contains only ones, Wilcoxon Signed Rank Statistics are
     calculated.}
  \item{rand}{Sets the random number generator into a reproducible state.
     Ignored if Wilcoxon rank sums are computed, or \code{X} contains
     no zeros.}
}
\details{
  If there are ties, then the ranks of the observations belonging to the same
  group of tied observations will be set to the maximum rank available for the
  corresponding group.
}
\value{
  A numeric vector containing Wilcoxon rank statistics for each row of \code{X}.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{wilc.stat}},\code{\link{wilc.ebam}}
}  

\keyword{htest}
