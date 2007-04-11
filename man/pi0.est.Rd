\name{pi0.est}
\alias{pi0.est}
\title{Estimation of the prior probability}
\description{
  Estimates the prior probability that a gene is not differentially
  expressed by the natural cubic splines based method of Storey and
  Tibshirani (2003).
}
\usage{
  pi0.est(p, lambda = seq(0, 0.95, 0.05), ncs.value = "max", 
      ncs.weights = NULL)
}
\arguments{
  \item{p}{a numeric vector containing the p-values of the genes}
  \item{lambda}{a numeric vector or value specifying the \eqn{\lambda}{lambda}
     values used in the estimation of the prior probability}
  \item{ncs.value}{a character string. Only used if \code{lambda} is a
     vector. Either \code{"max"} or \code{"paper"}. For details, see 
     \code{Details}}
  \item{ncs.weights}{a numerical vector of the same length as \code{lambda}
     containing the weights used in the natural cubic spline fit. By default
     no weights are used}
}
\details{
  For each value of \code{lambda}, \eqn{\pi_0(\lambda)}{pi0(lambda)} is
  computed by the number of p-values \code{p} larger than
  \eqn{\lambda}{lambda} divided by \eqn{(1-\lambda)/m}{(1-lambda)\m},
  where \eqn{m} is the length of \code{p}.
  
  If \code{lambda} is a value, \eqn{\pi_0(\lambda)}{pi0(lambda)} is the
  estimate for the prior probabiltity \eqn{\pi_0}{pi0} that a gene is
  not differentially expressed.
  
  If \code{lambda} is a vector, a natural cubic spline \eqn{h} with 3 degrees of
  freedom is fitted through the data points 
  \eqn{(\lambda,\pi_0(\lambda))}{(lambda,pi0(lambda))},
  where each point is weighed by \code{ncs.weights}. \eqn{\pi_0}{pi0} is estimated
  by \eqn{h(v)}{h(v)}, where \eqn{v=\max\{\lambda\}}{v=max\{lambda\}} if 
  \code{ncs.value="max"}, and \eqn{v=1}{v=1} if \code{ncs.value="paper"}.
}
\value{
  \item{p0}{the estimate of the prior probability that a gene is not
     differentially expressed}
  \item{spline.out}{the output of \code{smooth.spline} used in this function}
}
\references{
   Storey, J.D., and Tibshirani, R. (2003). Statistical Significance for 
   Genome-wide Studies. \emph{PNAS}, 100, 9440-9445.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
  \code{\link{SAM-class}},\code{\link{sam}},\code{\link{qvalue.cal}}
}
\examples{\dontrun{
  # Load the package multtest and the data of Golub et al. (1999)
  # contained in multtest.
  library(multtest)
  data(golub)

  # Perform a SAM analysis.
  sam.out<-sam(golub,golub.cl,B=100,rand=123)

  # Estimate the prior probability that a gene is not significant
  pi0.est(sam.out@p.value)
}}
\keyword{htest}
\keyword{smooth}
