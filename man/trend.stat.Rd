\name{trend.stat}
\alias{trend.stat}
\alias{trend.stat.data.frame}
\alias{trend.stat.list}
\alias{trend.stat.default}

\title{SAM Analysis of Linear Trend}
\description{
  Generates the required statistics for a Significance Analysis of Microarrays for a linear
  trend in (ordinal) data.
  
  In the two-class case, the Cochran-Armitage trend statistic is computed. 
  Otherwise, the statistic for the general test of trend described on page 87 of
  Agresti (2002) is determined.
  
  Should not be called directly, but via sam(..., method = trend.stat).
}
\usage{
\method{trend.stat}{default}(data, cl, catt = TRUE, approx = TRUE, B = 100, 
   B.more = 0.1, B.max = 50000, n.subset = 10, rand = NA, ...)
   
\method{trend.stat}{list}(data, cl, catt = TRUE, approx = TRUE, B = 100, 
   B.more = 0.1, B.max = 50000, n.subset = 10, rand = NA, ...)
}

\arguments{
  \item{data}{either a numeric matrix or data frame, or a list. If a matrix or data frame, then each row 
     must correspond to a variable (e.g., a SNP), and each column to a sample (i.e.\ an observation).
     The values in the matrix or data frame are interpreted as the scores for the different levels
     of the variables.
     
     If the number of observations is huge it is better to specify \code{data} as a list consisting
     of matrices, where each matrix represents one group and summarizes
     how many observations in this group show which level at which variable. The row and column names
     of all matrices must be identical and in the same order. The column names must be interpretable
     as numeric scores for the different levels of the variables. These matrices can, e.g.,
     be generated using the function \code{rowTables} from the package \pkg{scrime}. (It is recommended
     to use this function, as \code{trend.stat} has been made for using the output of \code{rowTables}.) 
     For details on how to specify this list, see the examples section on this man page, and the help for 
     \code{rowChisqMultiClass} in the package \pkg{scrime}.}
  \item{cl}{a numeric vector of length \code{ncol(data)} indicating to which classes
     the samples in the matrix or data frame \code{data} belongs. The values in \code{cl} must be interpretable 
     as scores for the different classes. Must be specified if \code{data} is a matrix or a data frame,
     whereas \code{cl} can but must not be specified if \code{data} is a list. If specified in the latter case,
     \code{cl} must have length \code{data}, i.e.\ one score for each of the matrices, and thus for each of
     the groups. If not specified, \code{cl} will be set to the integers between 1 and \eqn{c}, where \eqn{c}
     is the number of classes/matrices.}
  \item{catt}{should the Cochran-Armitage trend statistic be computed in the two-class case? If \code{FALSE},
     the trend statistic described on page 87 of Agresti (2002) is determined which differs by the factor
     \eqn{(n - 1) / n} from the Cochran-Armitage trend statistic.}
  \item{approx}{should the null distribution be approximated by the \eqn{\chi^2}{Chisquare}-distribution
     with one degree of freedom? If \code{FALSE}, a permutation method is used to estimate the null distribution.
     If \code{data} is a list, \code{approx} must currently be \code{TRUE}.}
  \item{B}{the number of permutations used in the estimation of the null distribution,
     and hence, in the computation of the expected \eqn{d}-values.}
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
  \item{...}{ignored.}
}

\value{
    A list containing statistics required by \code{sam}.
}
\references{
   Agresti, A.\ (2002). \emph{Categorical Data Analysis}. Wiley, Hoboken, NJ. 2nd Edition.

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{SAM-class}},\code{\link{sam}}, \code{\link{chisq.stat}}, \code{\link{trend.ebam}}
}
\examples{\dontrun{
  # Generate a random 1000 x 40 matrix consisting of the values
  # 1, 2, and 3, and representing 1000 variables and 40 observations.
  
  mat <- matrix(sample(3, 40000, TRUE), 1000)
  
  # Assume that the first 20 observations are cases, and the
  # remaining 20 are controls, and that the values 1, 2, 3 in mat
  # can be interpreted as scores for the different levels
  # of the variables represented by the rows of mat.
  
  cl <- rep(1:2, e=20)
  
  # Then an SAM analysis of linear trend can be done by
  
  out <- sam(mat, cl, method=trend.stat)
  out
  
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
  
  out2 <- sam(ltabs, method=trend.stat, approx=TRUE)
  out2
}}

\keyword{htest}
