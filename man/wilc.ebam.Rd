\name{wilc.ebam}
\alias{wilc.ebam}
\title{EBAM Analysis Using Wilcoxon Rank Statistics}
\description{
  Generates the required statistics for an Empirical Bayes Analysis of Microarrays
  analysis using standardized Wilcoxon rank statistics. 
  
  Should not be called directly, but via ebam(..., method = wilc.ebam).
}
\usage{
   wilc.ebam(data, cl, approx50 = TRUE, ties.method = c("min", "random", 
       "max"), use.offset = TRUE, df.glm = 5, rand = NA)
}
\arguments{
  \item{data}{a matrix or a data frame. Each row of
    \code{data} must correspond to a variable, and each column to a sample}
  \item{cl}{a numeric vector of length \code{ncol(data)} containing the class
     labels of the samples. In the two class paired case, \code{cl} can also 
     be a matrix with \code{ncol(data)} rows and 2 columns. For details
     on how \code{cl} should be specified, see \code{\link{ebam}}}
  \item{approx50}{if \code{TRUE}, the null distribution will be approximated by
     the standard normal distribution. Otherwise, the exact null distribution is
     computed. This argument will automatically be set to \code{FALSE} if there
     are less than 50 samples in each of the groups}
  \item{ties.method}{either \code{"min"} (default), \code{"random"}, or \code{"max"}. If
     \code{"random"}, the ranks of ties are randomly assigned. If \code{"min"} or \code{"max"},
     the ranks of ties are set to the minimum or maximum rank, respectively. For the handling
     of Zeros, see Details}
  \item{use.offset}{should an offset be used in the Poisson regression employed to estimate
     the density of the observed Wilcoxon rank sums? If \code{TRUE}, the log-transformed
     values of the null density is used as offset}
  \item{df.glm}{integer specifying the degrees of freedom of the natural cubic spline employed
     in the Poisson regression}  
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state}
}
\details{
  Standardized versions of the Wilcoxon rank statistics are computed. This means that
  \eqn{W*=(W-W_{mean})/W_{sd}}{W*=(W-mean(W))/sd(W)} is used as expression 
  score \eqn{z}, where \eqn{W} is the usual Wilcoxon rank sum statistic or Wilcoxon
  signed rank statistic, respectively. 
  
  In the computation of these statistics, the ranks of ties are by default set to the
  minimum rank. In the computation of the Wilcoxon signed rank statistic, zeros are randomly 
  set either to a very small positive or negative value.
  
  If there are less than 50 observations in each of the groups, the exact null distribution
  will be used. If there are more than 50 observations in at least one group, the null
  distribution will by default be approximated by the standard normal distribution. It is,
  however, still possible to compute the exact null distribution by  setting \code{approx50}
  to \code{FALSE}. 
}
\value{
   a list of statistics required by \code{ebam}
}

\references{
   Efron, B., Storey, J.D., Tibshirani, R.\ (2001). Microarrays, empirical Bayes methods, and
   the false discovery rate, \emph{Technical Report}, Department of Statistics, Stanford
   University.

   Schwender, H., Krause, A. and Ickstadt, K. (2003). Comparison of
   the Empirical Bayes and the Significance Analysis of Microarrays.
   \emph{Technical Report}, SFB 475, University of Dortmund, Germany.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{ebam}}
}  

\keyword{htest}
