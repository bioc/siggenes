\name{samControl}
\alias{samControl}
\alias{ebamControl}

\title{SAM and EBAM Arguments}
\description{
  Specifies most of the optional arguments of \code{sam} and \code{ebam}.
  
  Currently, only used in \code{limma2sam} and \code{limma2ebam}.
}
\usage{
samControl(delta = NULL, n.delta = 10, p0 = NA, lambda = seq(0, 0.95, 0.05), 
   ncs.value = "max", ncs.weights = NULL, q.version = 1)
   
ebamControl(p0 = NA, p0.estimation = c("splines", "interval", "adhoc"), 
   lambda = NULL, ncs.value = "max", use.weights = FALSE)

}

\arguments{
  \item{delta}{a numeric vector specifying a set of values for the threshold 
     \eqn{\Delta}{Delta} that should be used in \code{\link{sam}}. If \code{NULL}, \code{n.delta}
     \eqn{\Delta}{Delta} values will be computed automatically. Please note that the
     meaning of the argument \code{delta} differs between \code{sam} and \code{ebam}.}
  \item{n.delta}{a numeric value specifying the number of \eqn{\Delta}{Delta} values
     that will be computed over the range of all possible values for \eqn{\Delta}{Delta}
     in \code{sam}, if \code{delta} is not specified.}
  \item{p0}{a numeric value specifying the prior probability \eqn{\pi_0}{pi0} 
     that a gene is not differentially expressed. If \code{NA}, \code{p0} will
     be computed by the function \code{\link{pi0.est}}.}
  \item{p0.estimation}{either \code{"splines"} (default), \code{"interval"}, or \code{"adhoc"}. 
    If \code{"splines"}, the spline based method of Storey and Tibshirani (2003) is used to estimate
    \eqn{p_0}{p0}. If \code{"adhoc"} (\code{"interval"}), the adhoc (interval based) method 
    proposed by Efron et al.\ (2001) is used to estimate \eqn{p_0}{p0}. In the case of \code{samControl},
    \code{p0.estimation = "splines"} is always used. For details on this estimation, see \code{\link{pi0.est}}.}
  \item{lambda}{a numeric vector or value specifying the \eqn{\lambda}{lambda}
     values used in the estimation of the prior probability. For details, see
     \code{\link{pi0.est}}.}
  \item{ncs.value}{a character string. Only used if \code{lambda} is a
     vector. Either \code{"max"} or \code{"paper"}. For details, see \code{\link{pi0.est}}.}
  \item{ncs.weights}{a numeric vector of the same length as \code{lambda}
     containing the weights used in the estimation of \eqn{\pi_0}{pi0}. By default
     no weights are used. For details, see \code{\link{pi0.est}}. Only considered in \code{samControl}.}
  \item{use.weights}{should weights be used in the spline based estimation of \eqn{p_0}{p0}? If
    \code{TRUE}, 1 - \code{lambda} is used as weights. For details, see \code{\link{pi0.est}}.
    Only considered in \code{ebamControl}.}
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
  A list containing the values of the parameters that are used in \code{sam} or \code{ebam},
  respectively.
}
\references{   
   Efron, B., Tibshirani, R., Storey, J.D. and Tusher, V. (2001). Empirical Bayes Analysis
   of a Microarray Experiment. \emph{JASA}, 96, 1151-1160.
   
   Schwender, H., Krause, A., and Ickstadt, K. (2006). Identifying Interesting Genes with siggenes.
   \emph{RNews}, 6(5), 45-50.
   
   Storey, J.D. and Tibshirani, R. (2003). Statistical Significance for Genome-Wide
   Studies. \emph{Proceedings of the National Academy of Sciences}, 100, 9440-9445.

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.   
   
}
\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{limma2sam}}, \code{\link{limma2ebam}}, \code{\link{sam}}, \code{\link{ebam}}}

\keyword{utilities}

