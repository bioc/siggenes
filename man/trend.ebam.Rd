\name{trend.ebam}
\alias{trend.ebam}
\alias{trend.ebam.data.frame}
\alias{trend.ebam}
\alias{trend.ebam.list}
\alias{trend.ebam.default}

\title{EBAM Analysis of Linear Trend}
\description{
  Generates the required statistics for an Empirical Bayes Analysis of Microarrays for a linear
  trend in (ordinal) data.
  
  In the two-class case, the Cochran-Armitage trend statistic is computed. 
  Otherwise, the statistic for the general test of trend described on page 87 of
  Agresti (2002) is determined.
  
  Should not be called directly, but via ebam(..., method = trend.ebam).
}
\usage{
\method{trend.ebam}{default}(data, cl, catt = TRUE, approx = TRUE, n.interval = NULL,
    df.dens = NULL, knots.mode = NULL, type.nclass = "wand",
    B = 100, B.more = 0.1, B.max = 50000, n.subset = 10, 
    fast = FALSE, df.ratio = 3, rand = NA, ...)
    
\method{trend.ebam}{list}(data, cl, catt = TRUE, approx = TRUE, n.interval = NULL, 
    df.dens = NULL, knots.mode = NULL, type.nclass = "wand", ...)
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
  \item{n.interval}{the number of intervals used in the logistic regression with
     repeated observations for estimating the ratio \eqn{f_0/f}{f0/f} 
     (if \code{approx = FALSE}), or in the Poisson regression used to estimate
     the density of the observed \eqn{z}-values (if \code{approx = TRUE}).
     If \code{NULL}, \code{n.interval} is set to 139 if \code{approx = FALSE},
     and estimated by the method specified by \code{type.nclass} if \code{approx = TRUE}.}
  \item{df.dens}{integer specifying the degrees of freedom of the natural cubic
     spline used in the Poisson regression to estimate the density of the observed
     \eqn{z}-values. Ignored if \code{approx = FALSE}. 
     If \code{NULL}, \code{df.dens} is set to 3 if the degrees of freedom
     of the appromimated null distribution, i.e.\ the \eqn{\chi^2}{ChiSquare}-distribution,
     are less than or equal to 2, and otherwise \code{df.dens} is set to 5.}
  \item{knots.mode}{if \code{TRUE} the \code{df.dens} - 1 knots are centered around the
     mode and not the median of the density when fitting the Poisson regression model.
     Ignored if \code{approx = FALSE}. 
     If not specified, \code{knots.mode} is set to
     \code{TRUE} if the degrees of freedom of the approximated null distribution, i.e.\
     tht \eqn{\chi^2}{ChiSquare}-distribution, are larger than or equal to 3, and otherwise
     \code{knots.mode} is set to \code{FALSE}. For details on this density estimation, 
     see \code{\link{denspr}}.}
  \item{type.nclass}{character string specifying the procedure used to compute the
     number of cells of the histogram. Ignored if \code{approx = FALSE} or 
     \code{n.interval} is specified. Can be either
     \code{"wand"} (default), \code{"scott"}, or \code{"FD"}. For details, see
     \code{\link{denspr}}.}
  \item{B}{the number of permutations used in the estimation of the null distribution,
     and hence, in the computation of the expected \eqn{z}-values.}
  \item{B.more}{a numeric value. If the number of all possible permutations is smaller
     than or equal to (1+\code{B.more})*\code{B}, full permutation will be done. 
     Otherwise, \code{B} permutations are used.}
  \item{B.max}{a numeric value. If the number of all possible permutations is smaller
     than or equal to \code{B.max}, \code{B} randomly selected permutations will be used
     in the computation of the null distribution. Otherwise, \code{B} random draws
     of the group labels are used.} 
  \item{n.subset}{a numeric value indicating in how many subsets the \code{B} 
     permutations are divided when computing the permuted \eqn{z}-values. Please note
     that the meaning of \code{n.subset} differs between the SAM and the EBAM functions.}
  \item{fast}{if \code{FALSE} the exact number of permuted test scores that are
     more extreme than a particular observed test score is computed for each of
     the variables/SNPs. If \code{TRUE}, a crude estimate of this number is used.} 
  \item{df.ratio}{integer specifying the degrees of freedom of the natural cubic
     spline used in the logistic regression with repeated observations. Ignored
     if \code{approx = TRUE}.} 
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state.}
  \item{\dots}{ignored.}
}
\value{
  A list containing statistics required by \code{ebam}.
}


\references{
   Agresti, A.\ (2002). \emph{Categorical Data Analysis}. Wiley, Hoboken, 
   NJ. 2nd Edition.
   
   Efron, B., Tibshirani, R., Storey, J.D., and Tusher, V. (2001). 
   Empirical Bayes Analysis of a Microarray Experiment, \emph{JASA}, 
   96, 1151-1160.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{EBAM-class}},\code{\link{ebam}, \code{\link{trend.stat}}, \code{\link{chisq.ebam}}}
}

\examples{\dontrun{
  # Generate a random 1000 x 40 matrix consisting of the values
  # 1, 2, and 3, and representing 1000 variables and 40 observations.
  
  mat <- matrix(sample(3, 40000, TRUE), 1000)
  
  # Assume that the first 20 observations are cases, and the
  # remaining 20 are controls, and that the values 1, 2, 3 in
  # mat can be interpreted as scores for the different levels
  # of the variables.
  
  cl <- rep(1:2, e=20)
  
  # Then an EBAM analysis of linear trend can be done by
  
  out <- ebam(mat, cl, method=trend.ebam)
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
  
  # And the same EBAM analysis as above can then be 
  # performed by 
  
  out2 <- ebam(ltabs, method=trend.ebam)
  out2
}}

\keyword{htest}

