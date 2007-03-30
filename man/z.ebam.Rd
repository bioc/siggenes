\name{z.ebam}
\alias{z.ebam}
\alias{z.find}

\title{EBAM analysis Using t- or F-test}
\description{
  Computes the required statistics for an Empirical Bayes Analysis with a modified
  t- or F-test. Should not be called directly, but via 
  \code{ebam(..., method = z.ebam)} or \code{find.a0(..., method = z.find)},
  respectively.
}
\usage{
  z.ebam(data, cl, a0 = NULL, quan.a0 = NULL, B = 100, delta = 0.9, 
      var.equal = FALSE, B.more = 0.1, B.max = 30000, 
      gene.names = dimnames(data)[[1]], n.subset = 10, fast = FALSE, 
      n.interval = 139, df.ratio = NULL, rand = NA)
  
  z.find(data, cl, B = 100, var.equal = FALSE, B.more = 0.1, 
      B.max = 30000)
}

\arguments{
  \item{data}{a matrix, data frame or ExpressionSet object. Each row of \code{x}
    (or \code{exprs(x)}) must correspond to a gene, and each column to a sample}
  \item{cl}{a numeric vector of length \code{ncol(x)} containing the class
     labels of the samples. For details on how \code{cl} should be specified, 
     see \code{\link{ebam}}}
  \item{a0}{a numeric value specifying the fudge factor}
  \item{quan.a0}{a numeric value between 0 and 1 specifying the quantile of the
     standard deviations of the genes that is used as fudge factor}
  \item{B}{an integer indicating how many permutations should be used in the
     estimation of the null distribution}
  \item{delta}{a probability. All genes showing a posterior probability that is
     larger than or equal to \code{delta} are called differentially expressed}
  \item{var.equal}{should the ordinary t-statistic assuming equal group variances
     be computed? If \code{FALSE} (default), Welch's t-statistic will be computed}
  \item{B.more}{a numeric value. If the number of all possible permutations is smaller
     than or equal to (1+\code{B.more})*\code{B}, full permutation will be done. 
     Otherwise, \code{B} permutations are used. This avoids that \code{B} permutations
     will be used -- and not all permutations -- if the number of all possible 
     permutations is just a little larger than \code{B}}
  \item{B.max}{a numeric value. If the number of all possible permutations is smaller
     than or equal to \code{B.max}, \code{B} randomly selected permutations will be used
     in the computation of the null distribution. Otherwise, \code{B} random draws
     of the group labels are used. In the latter way of permuting, it is possible that
     some of the permutations are used more than once}
  \item{gene.names}{a character vector of length \code{nrow(x)} containing the
     names of the genes. By default the row names of \code{x} are used}
  \item{n.subset}{an integer specifying in how many subsets the \code{B} permutations
     should be split when computing the permuted test scores}
  \item{fast}{if \code{FALSE} the exact number of permuted test scores that are
     more extreme than a particular observed test score is computed for each of
     the genes. If \code{TRUE}, a crude estimate of this number is used}
  \item{n.interval}{the number of intervals used in the logistic regression with
     repeated observations for estimating the ratio \eqn{f/f_0}{f/f0}}
  \item{df.ratio}{integer specifying the degrees of freedom of the natural cubic
     spline used in the logistic regression with repeated observations}
  \item{rand}{integer. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state}
}
\value{
  a list of object required by find.a0 or ebam, respectively
}
\references{   
   Efron, B., Tibshirani, R., Storey, J.D. and Tusher, V. (2001). 
   Empirical Bayes Analysis of a Microarray Experiment, \emph{JASA}, 
   96, 1151-1160.
    
   Schwender, H., Krause, A. and Ickstadt, K. (2003). Comparison of
   the Empirical Bayes and the Significance Analysis of Microarrays.
   \emph{Technical Report}, SFB 475, University of Dortmund, Germany.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{\code{\link{ebam}}, \code{\link{find.a0}}}

\keyword{htest}