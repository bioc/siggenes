\name{fudge2}
\alias{fudge2}
\title{Fudge Factor}
\description{
  Computes the fudge factor as described by Tusher et al. (2001).
}
\usage{
  fudge2(r, s, alpha = seq(0, 1, 0.05), include.zero = TRUE)
}
\arguments{
  \item{r}{a numeric vector. The numerator of the test statistic computed
     for each gene is represented by one component of this vector}
  \item{s}{a numeric vector. Each component of this vector corresponds
     to the denominator of the test statistic of a gene}
  \item{alpha}{a numeric value or vector specifying quantiles of the
     \code{s} values. If \code{alpha} is numeric, this quantile of \code{s}
     will be used as fudge factor. Otherwise, the \code{alpha} quantile
     of the \code{s} values is computed that is optimal following the criterion
     of Tusher et al. (2001)}
  \item{include.zero}{If \code{TRUE}, \eqn{s_0=0}{s0=0} is also a possible
     choice for the fudge factor}
}


\value{
  \item{s.zero}{the value of the fudge factor \eqn{s_0}{s0}}
  \item{alpha.hat}{the optimal quantile of the \code{s} values. If 
     \eqn{s_0=0}{s0=0}, \code{alpha.hat} will not be returned}
  \item{vec.cv}{the vector of the coefficients of variations. 
     Following Tusher et al. (2001), the optimal \code{alpha} quantile 
     is given by the quantile that leads to the smallest CV of the
     modified test statistics}
  \item{msg}{a character string summarizing the most important 
     information about the fudge factor}
}

\references{
   Tusher, V., Tibshirani, R., and Chu, G. (2001). Significance Analysis of
   Microarrays Applied to the Ionizing Radiation Response. \emph{PNAS}, 98,
   5116-5121.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
   \code{\link{SAM-class}},\code{\link{sam}},\code{\link{sam.dstat}}
}

\keyword{optimize}
