\name{ebam}
\alias{ebam}

\title{Empirical Bayes Analysis of Microarrays}
\description{
  Performs an Empirical Bayes Analysis of Microarrays (EBAM). It is possible to perform
  one and two class analyses using either a modified t-statistic or a (standardized) 
  Wilcoxon rank statistic, and a multiclass analysis using a modified F-statistic. 
  Moreover, this function provides a EBAM procedure for categorical data such as SNP data
  and the possibility to employ an user-written score function.
}
\usage{
  ebam(x, cl, method = z.ebam, delta = 0.9, which.a0 = NULL, 
       control = ebamControl(), gene.names = dimnames(x)[[1]],
       ...)
}

\arguments{
  \item{x}{either a matrix, a data frame or an ExpressionSet object, or the output of \code{\link{find.a0}}, 
     i.e.\ an object of class \code{FindA0}. Can also be a list (if \code{method = chisq.ebam} or 
     \code{method = trend.ebam}). For the latter case, see \code{\link{chisq.ebam}}. 
     If \code{x} is not a FindA0 object, then each row of \code{x} (or
     \code{exprs(x)}, respectively) must correspond to a variable (e.g., a gene or a SNP), and each column to
     a sample.}
  \item{cl}{a specification of the class labels of the samples. Ignored if \code{x} is a \code{FindA0} object.
     Needs not to be specified if \code{x} is a list.
     
     Typically, \code{cl} is specified by a vector of length \code{ncol(x)}.  
     In the two class paired case, \code{cl} can also 
     be a matrix with \code{ncol(x)} rows and 2 columns. If \code{x} is
     an ExpressionSet object, \code{cl} can also be a character string naming the column 
     of \code{pData(x)} that contains the class labels of the samples.
     
     In the one-class case, \code{cl} should be a vector of 1's. 
     
     In the two class unpaired case, \code{cl} should be a vector containing 0's
     (specifying the samples of, e.g., the control group) and 1's (specifying,
     e.g., the case group). 
     
     In the two class paired case, \code{cl} can be either a numeric vector or a numeric matrix. 
     If it is a vector, then \code{cl} has to consist of the integers between -1 and 
     \eqn{-n/2} (e.g., before treatment group) and between 1 and \eqn{n/2} (e.g.,
     after treatment group), where \eqn{n} is the length of \code{cl} and \eqn{k}
     is paired with \eqn{-k}, \eqn{k=1,\dots,n/2}. If \code{cl} is a matrix, one
     column should contain -1's and 1's specifying, e.g., the before and the after
     treatment samples, respectively, and the other column should contain integer
     between 1 and \eqn{n/2} specifying the \eqn{n/2} pairs of observations.
     
     In the multiclass case and if \code{method = chisq.ebam} or \code{method = trend.ebam}, 
     \code{cl} should be a vector containing integers between 1 and \eqn{g}, where \eqn{g} is 
     the number of groups. In the two latter cases, \code{cl} needs not to be specified, if
     \code{x} is a list. For details, see \code{\link{chisq.ebam}}.
      
     For examples of how \code{cl} can be specified, see the manual of \pkg{siggenes}.}
  \item{method}{a character string or name specifying the method or function that should be
     used in the computation of the expression score \eqn{z}. 
     
     If \code{method = z.ebam},
     a modified t- or F-statistic, respectively, will be computed as proposed by Efron et al. (2001).
     
     If \code{method = wilc.ebam}, a (standardized) Wilcoxon sum / signed rank statistic will
     be used as expression score. 
     
     For an analysis of categorical data such as SNP data,
     \code{method} can be set to \code{\link{chisq.ebam}}. In this case, Pearson's Chi-squared statistic
     is computed for each row. 
     
     If the variables are ordinal and a trend test should be applied 
     (e.g., in the two-class case, the Cochran-Armitage trend test), \code{method = trend.ebam}
     can be employed.
     
     It is also possible to employ an user-written function for computing an user-specified
     expression score. For details, see the vignette of \pkg{siggenes}.}
  \item{delta}{a numeric vector consisting of probabilities for which the number of differentially
    expressed genes and the FDR should be computed, where a gene is called differentially expressed
    if its posterior probability is larger than \eqn{\Delta}{Delta}.}
  \item{which.a0}{an integer between 1 and the length of \code{quan.a0} of \code{\link{find.a0}}. If
    \code{NULL}, the suggested choice of \code{find.a0} is used. Ignored if \code{x} is a matrix, data
    frame or ExpressionSet object.}
  \item{control}{further arguments for controlling the EBAM analysis. For these arguments,
    see \code{\link{ebamControl}}.}
  \item{gene.names}{a vector of length \code{nrow(x)} specifying the names of the variables. By default,
     the row names of the matrix / data frame comprised by \code{x} are used.}
  \item{\dots}{further arguments of the specific EBAM methods. If \code{method = z.ebam}, see 
    \code{\link{z.ebam}}. If \code{method = wilc.ebam}, see \code{\link{wilc.ebam}}. If
    \code{method = chisq.ebam}, see \code{\link{chisq.ebam}}.}
}
\value{
  An object of class EBAM.
}
\references{ 
   Efron, B., Tibshirani, R., Storey, J.D. and Tusher, V. (2001). Empirical Bayes Analysis
   of a Microarray Experiment. \emph{JASA}, 96, 1151-1160.
    
   Schwender, H., Krause, A., and Ickstadt, K. (2006). Identifying Interesting Genes with siggenes.
   \emph{RNews}, 6(5), 45-50.
   
   Storey, J.D. and Tibshirani, R. (2003). Statistical Significance for Genome-Wide
   Studies. \emph{Proceedings of the National Academy of Sciences}, 100, 9440-9445.
}


\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{\code{\link{EBAM-class}}, \code{\link{find.a0}}, \code{\link{z.ebam}},
   \code{\link{wilc.ebam}}, \code{\link{chisq.ebam}}}
\examples{\dontrun{
  # Load the data of Golub et al. (1999) contained in the package multtest.
  data(golub)
  
  # golub.cl contains the class labels.
  golub.cl
  
  # Perform an EBAM analysis for the two class unpaired case assuming
  # unequal variances. Specify the fudge factor a0 by the suggested
  # choice of find.a0
  find.out <- find.a0(golub, golub.cl, rand = 123)
  ebam.out <- ebam(find.out)
  ebam.out
    
  # Since a0 = 0 leads to the largest number of genes (i.e. the suggested
  # choice of a0), the following leads to the same results as the above
  # analysis (but only if the random number generator, i.e. rand, is set
  # to the same number).
  ebam.out2 <- ebam(golub, golub.cl, a0 = 0, fast = TRUE, rand = 123)
  ebam.out2

  # If fast is set to TRUE in ebam, a crude estimate of the number of
  # falsely called genes is used (see the help file for z.ebam). This
  # estimate is always employed in find.a0. 
  # The exact number is used in ebam when performing
  ebam.out3 <- ebam(golub, golub.cl, a0 = 0, rand = 123)
  ebam.out3  

  # Since this is the recommended way, we use ebam.out3 at the end of
  # the Examples section for further analyses.



  # Perform an EBAM analysis for the two class unpaired case assuming
  # equal group variances. Set a0 = 0, and use B = 50 permutations
  # of the class labels.
  ebam.out4 <- ebam(golub, golub.cl, a0 = 0, var.equal = TRUE, B = 50,
     rand = 123)
  ebam.out4
    
  # Perform an EBAM analysis for the two class unpaired cased assuming
  # unequal group variances. Use the median (i.e. the 50\% quantile)
  # of the standard deviations of the genes as fudge factor a0. And
  # obtain the number of genes and the FDR if a gene is called 
  # differentially when its posterior probability is larger than
  # 0.95.
  ebam.out5 <- ebam(golub, golub.cl, quan.a0 = 0.5, delta = 0.95,
     rand = 123)
  ebam.out5
    
  # For the third analysis, obtain the number of differentially
  # expressed genes and the FDR if a gene is called differentially
  # expressed if its posterior probability is larger than 0.8, 0.85,
  # 0.9, 0.95.
  print(ebam.out3, c(0.8, 0.85, 0.9, 0.95))
    
  # Generate a plot of the posterior probabilities for delta = 0.9.
  plot(ebam.out3, 0.9)
    
  # Obtain the list of genes called differentially expressed if their
  # posterior probability is larger than 0.99, and gene-specific 
  # statistics for these variables such as their z-value and their
  # local FDR.
  summary(ebam.out3, 0.99)
}}
\keyword{htest}
