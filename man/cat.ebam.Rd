\name{cat.ebam}
\alias{cat.ebam}

\title{EBAM Analysis for Categorical Data}
\description{
  Generates the required statistics for an Empirical Bayes Analysis of Microarrays (EBAM)
  of categorical data such as SNP data.
  
  Should not be called directly, but via ebam(..., method = cat.ebam).
}

\usage{
  cat.ebam(data, cl, approx = FALSE, B = 100, check.levels = TRUE, 
    check.for.NN = FALSE, lev = NULL, B.more = 0.1, B.max = 50000, 
    n.subset = 10, fast = FALSE, n.interval = 139, df.ratio = 3,
    df.glm = NULL, rand = NA)
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
  \item{check.levels}{if \code{TRUE}, it will be checked if all variables/SNPs have
     the same number of levels/categories}
  \item{check.for.NN}{if \code{TRUE}, it will be checked if any of the genotypes
     is equal to "NN". Can be very time-consuming when the data set is high-dimensional}
  \item{lev}{numeric or character vector specifying the codings of the levels of the
     variables/SNPs. Must only be specified if the variables are not coded by the
     integers between 1 and the number of levels}
  \item{B.more}{a numeric value. If the number of all possible permutations is smaller
     than or equal to (1+\code{B.more})*\code{B}, full permutation will be done. 
     Otherwise, \code{B} permutations are used}
  \item{B.max}{a numeric value. If the number of all possible permutations is smaller
     than or equal to \code{B.max}, \code{B} randomly selected permutations will be used
     in the computation of the null distribution. Otherwise, \code{B} random draws
     of the group labels are used}  
  \item{n.subset}{a numeric value indicating how many permutations are considered
     simultaneously when computing the expected \eqn{z}-values}
  \item{fast}{if \code{FALSE} the exact number of permuted test scores that are
     more extreme than a particular observed test score is computed for each of
     the variables/SNPs. If \code{TRUE}, a crude estimate of this number is used}
  \item{n.interval}{the number of intervals used in the logistic regression with
     repeated observations for estimating the ratio \eqn{f_0/f}{f0/f} 
     (if \code{approx = FALSE}), or in the Poisson regression used to estimate
     the density of the observed \eqn{z}-values (if \code{approx = TRUE})}
  \item{df.ratio}{integer specifying the degrees of freedom of the natural cubic
     spline used in the logistic regression with repeated observations. Ignored
     if \code{approx = TRUE}} 
  \item{df.glm}{integer specifying the degrees of freedom of the natural cubic
     spline used in the Poisson regression to estimate the density of the observed
     \eqn{z}-values. If \code{NULL}, \code{df.glm} is set to \code{min\{df.chisq, 5\}},
     where \code{df.chisq} are the degrees of freedom of the \eqn{chi^2}{Chisquare}-distribution.
     Ignored if \code{approx = FALSE}}
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
