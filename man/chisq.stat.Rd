\name{chisq.stat}
\alias{chisq.stat}
\alias{cat.stat}

\title{SAM Analysis for Categorical Data}
\description{
  Generates the required statistics for a Significance Analysis of Microarrays of categorical
  data such as SNP data.
  
  Should not be called directly, but via sam(..., method = chisq.stat).
  
  Replaces \code{cat.stat}
}
\usage{
chisq.stat(data, cl, approx = NULL, B = 100, n.split = 1, 
   check.for.NN = FALSE, lev = NULL, B.more = 0.1, 
   B.max = 50000, n.subset = 10, rand = NA)
}

\arguments{
  \item{data}{a matrix, data frame, or list. If a matrix or data frame, then each row 
     must correspond to a variable (e.g., a SNP), and each column to a sample (i.e.\ an observation).
     If the number of observations is huge it is better to specify \code{data} as a list consisting
     of matrices, where each matrix represents one group and summarizes
     how many observations in this group show which level at which variable. These matrices can
     be generated using the function \code{rowTables} from the package \pkg{scrime}. For details on
     how to specify this list, see the examples section on this man page, and the help for 
     \code{rowChisqMultiClass} in the package \pkg{scrime}.}
  \item{cl}{a numeric vector of length \code{ncol(data)} indicating to which class
     a sample belongs. Must consist of the integers between 1 and \eqn{c}, where 
     \eqn{c} is the number of different groups. Needs only to be specified if \code{data}
     is a matrix or a data frame.}
  \item{approx}{should the null distribution be approximated by a \eqn{\chi^2}{ChiSquare}-distribution?
     Currently only available if \code{data} is a matrix or data frame. If not specified, 
     \code{approx = FALSE} is used, and the null distribution is estimated by employing a 
     permutation method.}
  \item{B}{the number of permutations used in the estimation of the null distribution,
     and hence, in the computation of the expected \eqn{d}-values.}
  \item{n.split}{number of chunks in which the variables are splitted in the computation
     of the values of the test statistic. Currently, only available if \code{approx = TRUE}
     and \code{data} is a matrix or data frame.
     By default, the test scores of all variables are calculated simultaneously.
     If the number of variables or observations is large, setting \code{n.split} to a
     larger value than 1 can help to avoid memory problems.}
  \item{check.for.NN}{if \code{TRUE}, it will be checked if any of the genotypes
     is equal to "NN". Can be very time-consuming when the data set is high-dimensional.}
  \item{lev}{numeric or character vector specifying the codings of the levels of the
     variables/SNPs. Can only be specified if \code{data} is a matrix or a data frame.
     Must only be specified if the variables are not coded by the
     integers between 1 and the number of levels. Can also be a list. In this case,
     each element of this list must be a numeric or character vector specifying the codings,
     where all elements must have the same length.}
  \item{B.more}{a numeric value. If the number of all possible permutations is smaller
     than or equal to (1+\code{B.more})*\code{B}, full permutation will be done. 
     Otherwise, \code{B} permutations are used.}
  \item{B.max}{a numeric value. If the number of all possible permutations is smaller
     than or equal to \code{B.max}, \code{B} randomly selected permutations will be used
     in the computation of the null distribution. Otherwise, \code{B} random draws
     of the group labels are used.}  
  \item{n.subset}{a numeric value indicating how many permutations are considered
     simultaneously when computing the expected \eqn{d}-values.}
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state.}
}
\details{
  For each SNP (or more general, categorical variable), 
  Pearson's Chi-Square statistic is computed to test if the distribution
  of the SNP differs between several groups.  Since only one null distribution is estimated
  for all SNPs as proposed in the original SAM procedure of Tusher et al. (2001) all SNPs must
  have the same number of levels/categories. 
}
\section{Warning}{This procedure will only work correctly if all SNPs/variables have the same
  number of levels/categories. Therefore, it is stopped when the number of levels differ between
  the variables.
}


\value{
  A list containing statistics required by \code{sam}.
}


\references{
   Schwender, H. (2005). Modifying Microarray Analysis Methods for Categorical Data -- SAM 
   and PAM for SNPs. In Weihs, C. and Gaul, W. (eds.), \emph{Classification -- The Ubiquitous Challenge}. 
   Springer, Heidelberg, 370-377.

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{SAM-class}},\code{\link{sam}, \code{\link{chisq.ebam}}, \code{\link{trend.stat}}}
}
\examples{\dontrun{
  # Generate a random 1000 x 40 matrix consisting of the values
  # 1, 2, and 3, and representing 1000 variables and 40 observations.
  
  mat <- matrix(sample(3, 40000, TRUE), 1000)
  
  # Assume that the first 20 observations are cases, and the
  # remaining 20 are controls.
  
  cl <- rep(1:2, e=20)
  
  # Then an SAM analysis for categorical data can be done by
  
  out <- sam(mat, cl, method=chisq.stat, approx=TRUE)
  out
  
  # approx is set to TRUE to approximate the null distribution
  # by the ChiSquare-distribution (usually, for such a small
  # number of observations this might not be a good idea
  # as the assumptions behind this approximation might not
  # be fulfilled).
  
  # The same results can also be obtained by employing
  # contingency tables, i.e. by specifying data as a list.
  # For this, we need to generate the tables summarizing
  # groupwise how many observations show which level at
  # which variable. These tables can be obtained by
  
  library(scrime)
  cases <- rowTables(mat[, cl==1])
  controls <- rowTables(mat[, cl==2])
  ltabs <- list(cases, controls)
  
  # And the same SAM analysis as above can then be 
  # performed by 
  
  out2 <- sam(ltabs, method=chisq.stat, approx=TRUE)
  out2
}}

\keyword{htest}
