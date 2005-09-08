\name{sam.snp}
\alias{sam.snp}
\title{SAM Analysis for Categorical Data}
\description{
  Performs a SAM (Significance Analysis of Microarrays) analysis for categorical
  data such a SNP data
}
\usage{
  sam.snp(data, cl, B = 1000, med = FALSE, delta = NULL, n.delta = 10, 
     p0 = NA, lambda = seq(0, 0.95, 0.05), ncs.value = "max", ncs.weights = NULL,
     gene.names = dimnames(data)[[1]], q.version = 1, na.replace = TRUE, 
     check.levels = TRUE, rand = NA)
}
\arguments{
  \item{data}{a matrix or data frame. Each row must correspond to a SNP, and
     each column to a sample}
  \item{cl}{a numeric vector of length \code{ncol(data)} indicating to which class
     a sample belongs. Recommended way of specifying \code{cl} is the use of the
     integers between 1 and \eqn{g}, where \eqn{g} is the number of different groups,
     or in the two-class case the use of 0's and 1's}
  \item{B}{the number of permutations used in the estimation of the null distribution}
  \item{med}{if \code{FALSE} (default), the mean number of falsely called SNPs 
     will be computed. Otherwise, the median number is calculated}
  \item{delta}{a numeric vector specifying a set of values for the threshold 
     \eqn{\Delta}{Delta} that should be used. If \code{NULL}, \code{n.delta}
     \eqn{\Delta}{Delta} values will be computed automatically}
  \item{n.delta}{a numeric value specifying the number of \eqn{\Delta}{Delta} values
     that will be computed over the range of possible values of \eqn{\Delta}{Delta}
     if \code{delta} is not specified}
  \item{p0}{a numeric value specifying the prior probability \eqn{\pi_0}{pi0} 
     that a SNP is not differentially expressed. If \code{NA}, \code{p0} will
     be computed by the function \code{pi0.est}}
  \item{lambda}{a numeric vector or value specifying the \eqn{\lambda}{lambda}
     values used in the estimation of the prior probability. For details, see
     \code{?pi0.est}}
  \item{ncs.value}{a character string. Only used if \code{lambda} is a
     vector. Either \code{"max"} or \code{"paper"}. For details, see \code{?pi0.est}}
  \item{ncs.weights}{a numerical vector of the same length as \code{lambda}
     containing the weights used in the estimation of \eqn{\pi_0}{pi0}. By default
     no weights are used. For details, see \code{?pi0.est}}
  \item{gene.names}{a character vector of length \code{nrow(data)} containing the
     names of the SNPs. By default the row names of \code{data} are used}
  \item{q.version}{a numeric value indicating which version of the q-value should
     be computed. If \code{q.version=2}, the original version of the q-value, i.e.
     min\{pFDR\}, will be computed. If \code{q.version=1}, min\{FDR\} will be used
     in the calculation of the q-value. Otherwise, the q-value is not computed.
     For details, see \code{?qvalue.cal}}
  \item{na.replace}{if \code{TRUE}, the missing values of a SNP will be replaced
     by random draws from the empirical distribution of that SNP}
  \item{check.levels}{if \code{TRUE}, it will be checked if all variables/SNPs have
     the same number of levels/categories}
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state}
}
\details{
  For each SNP, Pearson's Chi-Square statistic is computed to test if the distribution
  of the SNP differs between several groups. Since it is very likely that the assumptions
  for the Chi-square-approximation are not fulfilled a permutation based
  method is used to estimate the null distribution. Since only one null distribution is estimated
  for all SNPs as proposed in the original SAM procedure of Tusher et al. (2001) all SNPs must
  have the same number of levels/categories. 
}
\section{Warning}{This procedure will only work correctly if all SNPs/variables have the same
  number of levels/categories.}


\value{
  an object of class SAM
}
\note{ 
   SAM was deveoped by Tusher et al. (2001).
    
   !!! There is a patent pending for the SAM technology at Stanford University. !!!
}

\references{
   Schwender, H. (2004). Modifying Microarray Analysis Methods for 
   Categorical Data -- SAM and PAM for SNPs. To appear in: \emph{Proceedings
   of the the 28th Annual Conference of the GfKl}.

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{SAM-class}},\code{\link{sam}},\code{\link{sam.dstat}},
  \code{\link{sam.wilc}}
}
\keyword{htest}
