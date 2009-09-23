\name{samControl}
\alias{samControl}

\title{Further SAM Arguments}
\description{
  Specifies most of the optional arguments of \code{sam}.
}
\usage{
samControl(delta = NULL, n.delta = 10, p0 = NA, lambda = seq(0, 0.95, 0.05), 
   ncs.value = "max", ncs.weights = NULL, q.version = 1)
}

\arguments{

  \item{delta}{a numeric vector specifying a set of values for the threshold 
     \eqn{\Delta}{Delta} that should be used. If \code{NULL}, \code{n.delta}
     \eqn{\Delta}{Delta} values will be computed automatically.}
  \item{n.delta}{a numeric value specifying the number of \eqn{\Delta}{Delta} values
     that will be computed over the range of all possible values for \eqn{\Delta}{Delta}
     if \code{delta} is not specified.}
  \item{p0}{a numeric value specifying the prior probability \eqn{\pi_0}{pi0} 
     that a gene is not differentially expressed. If \code{NA}, \code{p0} will
     be computed by the function \code{\link{pi0.est}}.}
  \item{lambda}{a numeric vector or value specifying the \eqn{\lambda}{lambda}
     values used in the estimation of the prior probability. For details, see
     \code{\link{pi0.est}}.}
  \item{ncs.value}{a character string. Only used if \code{lambda} is a
     vector. Either \code{"max"} or \code{"paper"}. For details, see \code{\link{pi0.est}}.}
  \item{ncs.weights}{a numerical vector of the same length as \code{lambda}
     containing the weights used in the estimation of \eqn{\pi_0}{pi0}. By default
     no weights are used. For details, see \code{?pi0.est}.}
  \item{q.version}{a numeric value indicating which version of the q-value should
     be computed. If \code{q.version = 2}, the original version of the q-value, i.e.
     min\{pFDR\}, will be computed. If \code{q.version = 1}, min\{FDR\} will be used
     in the calculation of the q-value. Otherwise, the q-value is not computed.
     For details, see \code{\link{qvalue.cal}}.}
}
\details{
  These parameters should only be changed if they are fully understood.
}
\value{
  A list containing the values of the parameters that are used in \code{sam}.
}
\references{   
   Schwender, H., Krause, A., and Ickstadt, K. (2006). Identifying Interesting Genes with siggenes.
   \emph{RNews}, 6(5), 45-50.
   
   Storey, J.D. and Tibshirani, R. (2003). Statistical Significance for Genome-Wide
   Studies. \emph{Proceedings of the National Academy of Sciences}, 100, 9440-9445.

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.   
   
}
\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{limma2sam}}, \code{\link{sam}}}

\keyword{utilities}

