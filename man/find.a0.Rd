\name{find.a0}
\alias{find.a0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of the Fudge Factor}
\description{
  Suggests an optimal value for the fudge factor in an EBAM analysis as proposed
  by Efron et al. (2001).
}
\usage{
  find.a0(data, cl, fun = z.find, B = 100, delta = 0.9, quan.a0 = (0:5)/5, 
      include.zero = TRUE, gene.names = dimnames(data)[[1]], n.chunk = 5, 
      n.interval = 139, df.ratio = NULL, p0.estimation = c("splines", 
      "adhoc", "interval"), lambda = NULL, ncs.value = "max", 
      use.weights = FALSE, rand = NA, ...)
}

\arguments{
  \item{data}{a matrix, data frame or an ExpressionSet object. 
    Each row of \code{data} (or \code{exprs(data)}, respectively) must 
    correspond to a gene, and each column to a sample}
  \item{cl}{a numeric vector of length \code{ncol(data)} containing the class
     labels of the samples. In the two class paired case, \code{cl} can also 
     be a matrix with \code{ncol(data)} rows and 2 columns. If \code{data} is
     an ExpressionSet object, \code{cl} can also be a character string naming 
     the column of \code{pData(data)} that contains the class labels of the samples.
     
     In the one-class case, \code{cl} should be a vector of 1's. 
     
     In the two class unpaired case, \code{cl} should be a vector containing 0's
     (specifying the samples of, e.g., the control group) and 1's (specifying,
     e.g., the case group). 
     
     In the two class paired case, \code{cl} can be either a numeric vector or 
     a numeric matrix. If it is a vector, then \code{cl} has to consist of the 
     integers between -1 and \eqn{-n/2} (e.g., before treatment group) and between 
     1 and \eqn{n/2} (e.g., after treatment group), where \eqn{n} is the length of
     \code{cl} and \eqn{k} is paired with \eqn{-k}, \eqn{k=1,\dots,n/2}. If \code{cl} 
     is a matrix, one column should contain -1's and 1's specifying, e.g., the before 
     and the after treatment samples, respectively, and the other column should 
     contain integer between 1 and \eqn{n/2} specifying the \eqn{n/2} pairs of 
     observations.
     
     In the multiclass case and if \code{method=cat.stat}, \code{cl} should be a 
     vector containing integers between 1 and \eqn{g}, where \eqn{g} is the number 
     of groups.
        
     For examples of how \code{cl} can be specified, see the manual of \pkg{siggenes}}
  \item{fun}{the name of a function for computing the numerator and the denominator
     of the test statistic of interest, and for specifying other objects required
     for the identification of the fudge factor. The default function \code{z.find}
     provides these objects for t- and F-statistics. It is, however, also possible
     to employ an user-written function. For how to write such a function, see the
     vignette of \pkg{siggenes}}
  \item{B}{the number of permutations used in the estimation of the null distribution}
  \item{delta}{a probability. All genes showing a posterior probability that is
     larger than or equal to \code{delta} are called differentially expressed}
  \item{quan.a0}{a numeric vector indicating over which quantiles of the
     standard deviations of the genes the fudge factor \eqn{a_0}{a0} should be
     optimized}
  \item{include.zero}{should \eqn{a_0=0}{a0 = 0}, i.e. the not-modified test statistic
     also be a possible choice for the fudge factor?}
  \item{gene.names}{a character vector of length \code{nrow(data)} containing the
     names of the genes. By default the row names of \code{data} are used}
  \item{n.chunk}{an integer specifying in how many subsets the \code{B} permutations
     should be split when computing the permuted test scores}
  \item{n.interval}{the number of intervals used in the logistic regression with
     repeated observations for estimating the ratio \eqn{f_0/f}{f0/f}}
  \item{df.ratio}{integer specifying the degrees of freedom of the natural cubic
     spline used in the logistic regression with repeated observations}
  \item{p0.estimation}{either \code{"splines"} (default), \code{"interval"}, or \code{"adhoc"}. 
    If \code{"splines"}, the spline based method of Storey and Tibshirani (2003) is used to estimate
    \eqn{p_0}{p0}. If \code{\"adhoc"} (\code{"interval")}, the adhoc (interval based) 
    method proposed by Efron et al.\ (2001) is used to estimate \eqn{p_0}{p0}}
  \item{lambda}{a numeric vector or value specifying the \eqn{\lambda}{lambda} values used in
    the estimation of \eqn{p_0}{p0}. If \code{NULL}, \code{lambda} is set to \code{seq(0, 0.95, 0.05)}
    if \code{p0.estimation = "splines"}, and to \code{0.5} if \code{p0.estimation = "interval"}.
    Ignored if \code{p0.estimation = "adhoc"}. For details, see \code{\link{pi0.est}}}
  \item{ncs.value}{a character string. Only used if \code{p0.estimation = "splines"} and
    \code{lambda} is a vector. Either \code{"max"} or \code{"paper"}. For details, see
    \code{\link{pi0.est}}}
  \item{use.weights}{should weights be used in the spline based estimation of \eqn{p_0}{p0}? If
    \code{TRUE}, 1 - \code{lambda} is used as weights. For details, see \code{\link{pi0.est}}}
  \item{rand}{integer. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state}
  \item{\dots}{further arguments for the function specified by \code{fun}. For
     further arguments of \code{fun = z.find}, see \code{\link{z.find}}}
}
\details{
  The suggested choice for the fudge factor is the value of \eqn{a_0}{a0} that
  leads to the largest number of genes showing a posterior probability larger
  than \code{delta}. 
  
  Actually, only the genes having a posterior probability larger than \code{delta}
  are called differentially expressed that do not exhibit a test score less extreme
  than the score of a gene whose posterior probability is less than \code{delta}.
  So, let's say, we have done an EBAM analysis with a t-test and we have ordered
  the genes by their t-statistic. Let's further assume that Gene 1 to Gene 5 (i.e.
  the five genes with the lowest t-statistics), Gene 7 and 8, Gene 3012 to 3020, 
  and Gene 3040 to 3051 are the only genes that show a posterior probability larger
  than \code{delta}. Then, Gene 1 to 5, and 3040 to 3051 are called differentially
  expressed, but Gene 7 and 8, and 3012 to 3020 are not called differentially 
  expressed, since Gene 6 and Gene 3021 to 3039 show a posterior probability less
  than \code{delta}.
}
\value{
  an object of class FindA0
}

\note{
   The numbers of differentially expressed genes can differ between \code{find.a0}
   and \code{ebam}, even though the same value of the fudge factor is used, since
   in \code{find.a0} the observed and permuted test scores are monotonically
   transformed such that the observed scores follow a standard normal distribution
   (if the test statistic can take both positive and negative values) and
   an F-distribution (if the test statistic can only take positive values) for each
   possible choice of the fudge factor.
}

\references{  
   Efron, B., Tibshirani, R., Storey, J.D. and Tusher, V. (2001). 
   Empirical Bayes Analysis of a Microarray Experiment, \emph{JASA}, 
   96, 1151-1160.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{\code{\link{ebam}}, \code{\link{FindA0-class}}}
\examples{\dontrun{
  # Load the data of Golub et al. (1999) contained in the package multtest.
  data(golub)
  
  # golub.cl contains the class labels.
  golub.cl
  
  # Obtain the number of differentially expressed genes and the FDR for the
  # default set of values for the fudge factor.
  find.out <- find.a0(golub, golub.cl, rand = 123)
  find.out
  
  # Obtain the number of differentially expressed genes and the FDR when using
  # the t-statistic assuming equal group variances
  find.out2 <- find.a0(golub, golub.cl, var.equal = TRUE, rand = 123)
  
  # Using the Output of the first analysis with find.a0, the number of 
  # differentially expressed genes and the FDR for other values of
  # delta, e.g., 0.95, can be obtained by
  print(find.out, 0.95)
  
  # The logit-transformed posterior probabilities can be plotted by
  plot(find.out)
  
  # To avoid the logit-transformation, set logit = FALSE.
  plot(find.out, logit = FALSE)
}}
\keyword{htest}
