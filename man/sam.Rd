\name{sam}
\alias{sam}
\title{Significance Analysis of Microarray}
\description{
  Performs a Significance Analysis of Microarrays (SAM). It is possible to perform
  one and two class analyses using either a modified t-statistic or a (standardized) 
  Wilcoxon rank statistic, and a multiclass analysis using a modified F-statistic. 
  Moreover, this function provides a SAM procedure for categorical data such as SNP data
  and the possibility to employ an user-written score function.
}
\usage{
  sam(data, cl, method = d.stat, control=samControl(),
      gene.names = dimnames(data)[[1]], ...)
}
\arguments{
  \item{data}{a matrix, a data frame, or an ExpressionSet object. Each row of \code{data}
    (or \code{exprs(data)}, respectively) must correspond to a variable (e.g., a gene), and
    each column to a sample (i.e.\ an observation).
    
    Can also be a list (if \code{method = chisq.stat} or 
     \code{method = trend.stat}). For details on how to specify data in this case, 
     see \code{\link{chisq.stat}}.}
  \item{cl}{a vector of length \code{ncol(data)} containing the class
     labels of the samples. In the two class paired case, \code{cl} can also 
     be a matrix with \code{ncol(data)} rows and 2 columns. If \code{data} is
     an ExpressionSet object, \code{cl} can also be a character string naming the column
     of \code{pData(data)} that contains the class labels of the samples. If \code{data}
     is a list, \code{cl} needs not to be specified. 
     
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
     
     In the multiclass case and if \code{method = chisq.stat}, \code{cl} should be a vector containing integers
     between 1 and \eqn{g}, where \eqn{g} is the number of groups. (In the case of \code{\link{chisq.stat}},
     \code{cl} needs not to be specified if \code{data} is a list of groupwise matrices.)
        
     For examples of how \code{cl} can be specified, see the manual of \pkg{siggenes}.}
  \item{method}{a character string or a name specifying the method/function that should be used
     in the computation of the expression scores \eqn{d}. 
     
     If \code{method = d.stat},
     a modified t-statistic or F-statistic, respectively, will be computed
     as proposed by Tusher et al. (2001). 
     
     If \code{method = wilc.stat}, a
     Wilcoxon rank sum statistic or Wilcoxon signed rank statistic will be used
     as expression score. 
     
     For an analysis of categorical data such as SNP data, 
     \code{method} can be set to \code{chisq.stat}. In this case Pearson's
     ChiSquare statistic is computed for each row. 
     
     If the variables are ordinal and a trend test should be applied 
     (e.g., in the two-class case, the Cochran-Armitage trend test), \code{method = trend.stat}
     can be employed.
     
     It is also possible to use
     an user-written function to compute the expression scores.
     For details, see \code{Details}.}
  \item{control}{further optional arguments for controlling the SAM analysis. For
     these arguments, see \code{\link{samControl}}.}
  \item{gene.names}{a character vector of length \code{nrow(data)} containing the
     names of the genes. By default the row names of \code{data} are used.}
  \item{\dots}{further arguments of the specific SAM methods. If \code{method = d.stat},
     see the help of \code{\link{d.stat}}. If \code{method = wilc.stat}, see the help
     of \code{\link{wilc.stat}}. If \code{method = chisq.stat}, see the help of
     \code{\link{chisq.stat}}.}
}
\details{
  \code{sam} provides SAM procedures for several types of analysis (one and two class analyses
  with either a modified t-statistic or a Wilcoxon rank statistic, a multiclass analysis
  with a modified F statistic, and an analysis of categorical data). It is, however, also 
  possible to write your own function for another type of analysis. The required arguments
  of this function must be \code{data} and \code{cl}. This function can also have other
  arguments. The output of this function must be a list containing the following objects:
  \describe{
     \item{\code{d}:}{a numeric vector consisting of the expression scores of the genes.}
     \item{\code{d.bar}:}{a numeric vector of the same length as \code{na.exclude(d)} specifying
        the expected expression scores under the null hypothesis.}
     \item{\code{p.value}:}{a numeric vector of the same length as \code{d} containing
        the raw, unadjusted p-values of the genes.}
     \item{\code{vec.false}:}{a numeric vector of the same length as \code{d} consisting of
        the one-sided numbers of falsely called genes, i.e. if \eqn{d > 0} the numbers
        of genes expected to be larger than \eqn{d} under the null hypothesis, and if
        \eqn{d<0}, the number of genes expected to be smaller than \eqn{d} under the
        null hypothesis.}
     \item{\code{s}:}{a numeric vector of the same length as \code{d} containing the standard deviations 
        of the genes. If no standard deviation can be calculated, set \code{s = numeric(0)}.}
     \item{\code{s0}:}{a numeric value specifying the fudge factor. If no fudge factor is calculated,
        set \code{s0 = numeric(0)}.}
     \item{\code{mat.samp}:}{a matrix with B rows and \code{ncol(data)} columns, where B is the number
        of permutations, containing the permutations used in the computation of the permuted
        d-values. If such a matrix is not computed, set \code{mat.samp = matrix(numeric(0))}.}
     \item{\code{msg}:}{a character string or vector containing information about, e.g., which type of analysis
        has been performed. \code{msg} is printed when the function \code{print} or 
        \code{summary}, respectively, is called. If no such message should be printed, set \code{msg = ""}.}
     \item{\code{fold}:}{a numeric vector of the same length as \code{d} consisting of the fold 
        changes of the genes. If no fold change has been computed, set \code{fold = numeric(0)}.}
  }
  If this function is, e.g., called \code{foo}, it can be used by setting \code{method = foo}
  in \code{sam}. More detailed information and an example will be contained in the siggenes
  manual.   
}
\value{
  An object of class SAM.
}

\references{
   Schwender, H., Krause, A., and Ickstadt, K. (2006). Identifying Interesting Genes with siggenes.
   \emph{RNews}, 6(5), 45-50.

   Schwender, H. (2004). Modifying Microarray Analysis Methods for 
   Categorical Data -- SAM and PAM for SNPs. To appear in: \emph{Proceedings
   of the the 28th Annual Conference of the GfKl}.

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.

}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{SAM-class}},\code{\link{d.stat}},\code{\link{wilc.stat}},
  \code{\link{chisq.stat}, \code{\link{samControl}}}
}
\examples{\dontrun{
  # Load the package multtest and the data of Golub et al. (1999)
  # contained in multtest.
  library(multtest)
  data(golub)
  
  # golub.cl contains the class labels.
  golub.cl

  # Perform a SAM analysis for the two class unpaired case assuming
  # unequal variances.
  sam.out <- sam(golub, golub.cl, B=100, rand=123)
  sam.out
  
  # Obtain the Delta plots for the default set of Deltas
  plot(sam.out)
  
  # Generate the Delta plots for Delta = 0.2, 0.4, 0.6, ..., 2
  plot(sam.out, seq(0.2, 0.4, 2))
  
  # Obtain the SAM plot for Delta = 2
  plot(sam.out, 2)
  
  # Get information about the genes called significant using 
  # Delta = 3.
  sam.sum3 <- summary(sam.out, 3, entrez=FALSE)
  
  # Obtain the rows of golub containing the genes called
  # differentially expressed
  sam.sum3@row.sig.genes
  
  # and their names
  golub.gnames[sam.sum3@row.sig.genes, 3] 

  # The matrix containing the d-values, q-values etc. of the
  # differentially expressed genes can be obtained by
  sam.sum3@mat.sig
  
  # Perform a SAM analysis using Wilcoxon rank sums
  sam(golub, golub.cl, method="wilc.stat", rand=123)
    

  # Now consider only the first ten columns of the Golub et al. (1999)
  # data set. For now, let's assume the first five columns were
  # before treatment measurements and the next five columns were
  # after treatment measurements, where column 1 and 6, column 2
  # and 7, ..., build a pair. In this case, the class labels
  # would be
  new.cl <- c(-(1:5), 1:5)
  new.cl
  
  # and the corresponding SAM analysis for the two-class paired
  # case would be
  sam(golub[,1:10], new.cl, B=100, rand=123)
  
  # Another way of specifying the class labels for the above paired
  # analysis is
  mat.cl <- matrix(c(rep(c(-1, 1), e=5), rep(1:5, 2)), 10)
  mat.cl
  
  # and the above SAM analysis can also be done by
  sam(golub[,1:10], mat.cl, B=100, rand=123)
}}
\keyword{htest}
