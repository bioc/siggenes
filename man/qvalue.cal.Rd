\name{qvalue.cal}
\alias{qvalue.cal}
\title{Computation of the q-value}
\description{
  Computes the q-values of a given set of p-values.
}
\usage{
  qvalue.cal(p, p0, version = 1)
}
\arguments{
  \item{p}{a numeric vector containing the p-values.}
  \item{p0}{a numeric value specifying the prior probability
     that a gene is not differentially expressed.}
  \item{version}{If \code{version=2}, the original version of
     the q-value, i.e. min\{pFDR\}, will be computed. if
     \code{version=1}, min\{FDR\} will be used in the computation
     of the q-value.}
}
\details{
  Using \code{version = 1} in \code{qvalue.cal} corresponds to setting
  \code{robust = FALSE} in the function \code{qvalue} of John Storey's
  \R package \pkg{qvalue}, while \code{version = 2} corresponds to
  \code{robust = TRUE}.
}
\value{
  A vector of the same length as \code{p} containing the q-values
  corresponding to the p-values in \code{p}.
}
\references{
  Storey, J.D. (2003). The positive False Discovery Rate: A Bayesian
  Interpretation and the q-value. \emph{Annals of Statistics}, 31,
  2013-2035.

  Storey, J.D., and Tibshirani, R. (2003). Statistical Significance for 
  Genome-wide Studies. \emph{PNAS}, 100, 9440-9445.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
  \code{\link{pi0.est}},\code{\link{SAM-class}},\code{\link{sam}}
}
\examples{\dontrun{
  # Load the package multtest and the data of Golub et al. (1999)
  # contained in multtest.
  library(multtest)
  data(golub)

  # Perform a SAM analysis.
  sam.out<-sam(golub, golub.cl, B=100, rand=123)

  # Estimate the prior probability that a gene is not significant.
  pi0 <- pi0.est(sam.out@p.value)$p0
  
  # Compute the q-values of the genes.
  q.value <- qvalue.cal(sam.out@p.value, pi0)
}}
\keyword{htest}
