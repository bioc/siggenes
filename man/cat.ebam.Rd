\name{cat.ebam}
\alias{cat.ebam}

\title{EBAM Analysis for Categorical Data}
\description{
  Generates the required statistics for an Empirical Bayes Analysis of Microarrays (EBAM)
  of categorical data such as SNP data.
  
  Should not be called directly, but via ebam(..., method = cat.ebam).
}

\usage{
  cat.ebam(data, cl, approx = FALSE, B = 100, n.split = 1,
    check.levels = TRUE, check.for.NN = FALSE, lev = NULL, 
    B.more = 0.1, B.max = 50000, n.subset = 10, fast = FALSE, 
    n.interval = NULL, df.ratio = 3, df.dens = NULL, 
    knots.mode = NULL, type.nclass = "wand", rand = NA)
}


\arguments{
  \item{data}{a matrix or data frame. Each row must correspond to a variable/SNP, and
     each column to a sample}
  \item{cl}{a numeric vector of length \code{ncol(data)} indicating to which class
     a sample belongs. Must consist of the
     integers between 1 and \eqn{c}, where \eqn{c} is the number of different groups}
  \item{approx}{should the null distribution be approximated by a \eqn{\chi^2}{Chisquare}-distribution?}
  \item{B}{the number of permutations used in the estimation of the null distribution,
     and hence, in the computation of the expected \eqn{z}-values.}
  \item{n.split}{number of chunks in which the variables are splitted in the computation
     of the values of the test statistic. Currently, only available if \code{approx = TRUE}.
     By default, the test scores of all variables are calculated simultaneously.
     If the number of variables or observations is large, setting \code{n.split} to a
     larger value than 1 can help to avoid memory problems}
  \item{check.levels}{if \code{TRUE}, it will be checked if all variables/SNPs have
     the same number of levels/categories}
  \item{check.for.NN}{if \code{TRUE}, it will be checked if any of the genotypes
     is equal to "NN". Can be very time-consuming when the data set is high-dimensional}
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
  \item{n.subset}{a numeric value indicating in how many subsets the \code{B} 
     permutations are divided when computing the permuted \eqn{z}-values. Please note
     that the meaning of \code{n.subset} differs between the SAM and the EBAM functions}
  \item{fast}{if \code{FALSE} the exact number of permuted test scores that are
     more extreme than a particular observed test score is computed for each of
     the variables/SNPs. If \code{TRUE}, a crude estimate of this number is used}
  \item{n.interval}{the number of intervals used in the logistic regression with
     repeated observations for estimating the ratio \eqn{f_0/f}{f0/f} 
     (if \code{approx = FALSE}), or in the Poisson regression used to estimate
     the density of the observed \eqn{z}-values (if \code{approx = TRUE}).
     If \code{NULL}, \code{n.interval} is set to 139 if \code{approx = FALSE},
     and estimated by the method specified by \code{type.nclass} if \code{approx = TRUE}}
  \item{df.ratio}{integer specifying the degrees of freedom of the natural cubic
     spline used in the logistic regression with repeated observations. Ignored
     if \code{approx = TRUE}} 
  \item{df.dens}{integer specifying the degrees of freedom of the natural cubic
     spline used in the Poisson regression to estimate the density of the observed
     \eqn{z}-values. Ignored if \code{approx = FALSE}. 
     If \code{NULL}, \code{df.dens} is set to 3 if the degrees of freedom
     of the appromimated null distribution, i.e.\ the \eqn{\chi^2}{ChiSquare}-distribution,
     are less than or equal to 2, and otherwise \code{df.dens} is set to 5}
  \item{knots.mode}{if \code{TRUE} the \code{df.dens} - 1 knots are centered around the
     mode and not the median of the density when fitting the Poisson regression model.
     Ignored if \code{approx = FALSE}. 
     If not specified, \code{knots.mode} is set to
     \code{TRUE} if the degrees of freedom of the approximated null distribution, i.e.\
     tht \eqn{\chi^2}{ChiSquare}-distribution, are larger than or equal to 3, and otherwise
     \code{knots.mode} is set to \code{FALSE}. For details on this density estimation, 
     see \code{\link{denspr}}}
  \item{type.nclass}{character string specifying the procedure used to compute the
     number of cells of the histogram. Ignored if \code{approx = FALSE} or 
     \code{n.interval} is specified. Can be either
     \code{"wand"} (default), \code{"scott"}, or \code{"FD"}. For details, see
     \code{\link{denspr}}}
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state}
}
\details{
  For each variable, Pearson's Chi-Square statistic is computed to test if the distribution
  of the variable differs between several groups.  Since only one null distribution is estimated
  for all variables as proposed in the original EBAM application of Efron et al. (2001),
  all variables must have the same number of levels/categories. 
}
\section{Warning}{This procedure will only work correctly if all SNPs/variables have the same
  number of levels/categories.}


\value{
  a list containing statistics required by \code{ebam}
}


\references{
   Efron, B., Tibshirani, R., Storey, J.D. and Tusher, V. (2001). 
   Empirical Bayes Analysis of a Microarray Experiment, \emph{JASA}, 
   96, 1151-1160.
   
   Schwender, H. (2007). Empirical Bayes Analysis of Single Nucleotide
   Polymorphisms. \emph{Technical Report}, Department of Statistics,
   University of Dortmund. To appear soon.
   
   Schwender, H., Krause, A. and Ickstadt, K. (2003). Comparison of
   the Empirical Bayes and the Significance Analysis of Microarrays.
   \emph{Technical Report}, SFB 475, University of Dortmund, Germany.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{SAM-class}},\code{\link{sam}}
}
\keyword{htest}
