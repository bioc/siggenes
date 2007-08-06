\name{cat.stat}
\alias{sam.snp}
\alias{cat.stat}

\title{SAM Analysis for Categorical Data}
\description{
  Generates the required statistics for a Significance Analysis of Microarrays of categorical
  data such as SNP data.
  
  Should not be called directly, but via sam(..., method = cat.stat).
}

\usage{
  cat.stat(data, cl, B = 100, approx = FALSE, n.split = 1,
    check.levels = TRUE, check.for.NN = FALSE, lev = NULL, 
    B.more = 0.1, B.max = 50000, n.subset = 10, rand = NA)
}


\arguments{
  \item{data}{a matrix or data frame. Each row must correspond to a variable/SNP, and
     each column to a sample}
  \item{cl}{a numeric vector of length \code{ncol(data)} indicating to which class
     a sample belongs. Must consist of the
     integers between 1 and \eqn{c}, where \eqn{c} is the number of different groups}
  \item{B}{the number of permutations used in the estimation of the null distribution,
     and hence, in the computation of the expected \eqn{d}-values}
  \item{approx}{should the null distribution be approximated by the \eqn{\chi^2}{Chisquare}-distribution?}
  \item{n.split}{number of chunks in which the variables are splitted in the computation
     of the values of the test statistic. Currently, only available if \code{approx = TRUE}.
     By default, the test scores of all variables are calculated simultaneously.
     If the number of variables or observations is large, setting \code{n.split} to a
     larger value than 1 can help to avoid memory problems}
  \item{check.levels}{if \code{TRUE}, it will be checked if all variables/SNPs have
     the same number of levels/categories}
  \item{check.for.NN}{if \code{TRUE}, it will be checked if any of the genotypes
     is equal to "NN". Might be very time-consuming when the data set is high-dimensional}
  \item{lev}{numeric or character vector specifying the codings of the levels of the
     variables/SNPs. Must only be specified if the variables are not coded by the
     integers between 1 and the number of levels. Can also be a list. In this case,
     each element of this list must be a numeric or character vector specifying the codings,
     where all elements must have the same length}
  \item{B.more}{a numeric value. If the number of all possible permutations is smaller
     than or equal to (1+\code{B.more})*\code{B}, full permutation will be done. 
     Otherwise, \code{B} permutations are used}
  \item{B.max}{a numeric value. If the number of all possible permutations is smaller
     than or equal to \code{B.max}, \code{B} randomly selected permutations will be used
     in the computation of the null distribution. Otherwise, \code{B} random draws
     of the group labels are used}  
  \item{n.subset}{a numeric value indicating how many permutations are considered
     simultaneously when computing the expected \eqn{d}-values}
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state}
}
\details{
  For each SNP, Pearson's Chi-Square statistic is computed to test if the distribution
  of the SNP differs between several groups.  Since only one null distribution is estimated
  for all SNPs as proposed in the original SAM procedure of Tusher et al. (2001) all SNPs must
  have the same number of levels/categories. 
}
\section{Warning}{This procedure will only work correctly if all SNPs/variables have the same
  number of levels/categories.}


\value{
  a list containing statistics required by \code{sam}
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
  \code{\link{SAM-class}},\code{\link{sam}}
}
\keyword{htest}
