\name{sam.wilc}
\alias{sam.wilc}
\title{SAM Analysis Using Wilcoxon Rank Statistics}
\description{
  Performs a SAM (Significance Analysis of Microarrays) analysis using
  standardized Wilcoxon rank statistics. In the two class unpaired analysis, the 
  standardized Wilcoxon rank sum statistic is computed, 
  while in the one class analysis and in the the two class paired analysis, 
  the standardized Wilcoxon signed rank statistic is used as expression score.
}
\usage{
   sam.wilc(data, cl, delta = NULL, n.delta = 10, p0 = NA,
       lambda = seq(0, 0.95, 0.05), ncs.value = "max", ncs.weights = NULL,
       gene.names = dimnames(data)[[1]], q.version = 1, R.fold = 1,
       R.unlog = TRUE, na.replace = TRUE, na.method = "mean", approx50 = TRUE,
       check.ties = FALSE, rand = NA)
}
\arguments{
  \item{data}{a matrix, data frame or exprSet object. Each row of
    \code{data} (or \code{exprs(data)}, respectively) must correspond to a gene,
    and each column to a sample}
  \item{cl}{a numeric vector of length \code{ncol(data)} containing the class
     labels of the samples. In the two class paired case, \code{cl} can also 
     be a matrix with \code{ncol(data)} rows and 2 columns. If \code{data} is
     a exprSet object, \code{cl} can also be a character string. For details
     on how \code{cl} should be specified, see \code{?sam}}
  \item{delta}{a numeric vector specifying a set of values for the threshold 
     \eqn{\Delta}{Delta} that should be used. If \code{NULL}, \code{n.delta}
     \eqn{\Delta}{Delta} values will be computed automatically}
  \item{n.delta}{a numeric value specifying the number of \eqn{\Delta}{Delta} values
     that will be computed over the range of possible values of \eqn{\Delta}{Delta}
     if \code{delta} is not specified}
  \item{p0}{ a numeric value specifying the prior probability \eqn{\pi_0}{pi0} 
     that a gene is not differentially expressed. If \code{NA}, \code{p0} will
     be computed by the function \code{pi0.est}}
  \item{lambda}{a numeric vector or value specifying the \eqn{\lambda}{lambda}
     values used in the estimation of the prior probability. For details, see
     \code{?pi0.est}}
  \item{ncs.value}{a character string. Only used if \code{lambda} is a
     vector. Either \code{"max"} or \code{"paper"}. For details, see \code{?pi0.est}}
  \item{ncs.weights}{a numerical vector of the same length as \code{lambda}
     containing the weights used in the estimation of \eqn{\pi_0}{pi0}. By default
     no weights are used. For details, see \code{?pi0.est}}
  \item{gene.names}{a character vector of length \code{nrow(data)} containing the
     names of the genes. By default the row names of \code{data} are used}
  \item{q.version}{a numeric value indicating which version of the q-value should
     be computed. If \code{q.version=2}, the original version of the q-value, i.e.
     min\{pFDR\}, will be computed. If \code{q.version=1}, min\{FDR\} will be used
     in the calculation of the q-value. Otherwise, the q-value is not computed.
     For details, see \code{?qvalue.cal}}
  \item{R.fold}{a numeric value. If the fold change of a gene is smaller than or
     equal to \code{R.fold}, or larger than or equal to 1/\code{R.fold},respectively,
     then this gene will be excluded from the SAM analysis. The expression score 
     \eqn{d} of excluded genes is set to \code{NA}. By default, \code{R.fold}
     is set to 1 such that all genes are included in the SAM analysis. Setting 
     \code{R.fold} to 0 or a negative value will avoid the computation of the fold
     change. The fold change is only computed in the two-class cases}
  \item{R.unlog}{if \code{TRUE}, the anti-log of \code{data} will be used in the computation of the
     fold change. Otherwise, \code{data} is used. This transformation should be done
     if \code{data} is log2-tranformed (in a SAM analysis it is highly recommended
     to use log2-transformed expression data)}
  \item{na.replace}{if \code{TRUE}, missing values will be removed by the genewise/rowwise
     statistic specified by \code{na.method}. If a gene has less than 2 non-missing
     values, this gene will be excluded from further analysis. If \code{na.replace=FALSE},
     all genes with one or more missing values will be excluded from further analysis.
     The expression score \eqn{d} of excluded genes is set to \code{NA}}
  \item{na.method}{a character string naming the statistic with which missing values
     will be replaced if \code{na.replace=TRUE}. Must be either \code{"mean"} (default)
     or \code{median}}
  \item{approx50}{if \code{TRUE}, the null distribution will be approximated by
     the standard normal distribution. Otherwise, the exact null distribution is
     computed. This argument will automatically be set to \code{FALSE} if there
     are less than 50 samples in each of the groups}
  \item{check.ties}{if \code{TRUE}, a warning will be generated if there are ties or Zeros.
     This warning contains information about how many genes have ties or Zeros. Otherwise,
     this warning is not generated. Default is \code{FALSE} since checking for ties can
     take some time}
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state}
}
\details{
  Standardized versions of the Wilcoxon rank statistics are computed. This means that
  \eqn{W*=(W-W_{mean})/W_{sd}}{W*=(W-mean(W))/sd(W)} is used as expression 
  score \eqn{d}, where \eqn{W} is the usual Wilcoxon rank sum statistic or Wilcoxon
  signed rank statistic, respectively. 
  
  In the computation of these statistics, the ranks of ties are randomly
  assigned. In the computation of the Wilcoxon signed rank statistic, Zeros are randomly 
  set either to a very small positive or negative value.
  
  If there are less than 50 observations in each of the groups, the exact null distribution
  will be used. If there are more than 50 observations in at least one group, the null
  distribution will by default be approximated by the standard normal distribution. It is,
  however, still possible to compute the exact null distribution by  setting \code{approx50}
  to \code{FALSE}. 
}
\value{
   an object of class SAM
}
\note{ 
   SAM was deveoped by Tusher et al. (2001).
    
   !!! There is a patent pending for the SAM technology at Stanford University. !!!
}

\references{
   Schwender, H., Krause, A. and Ickstadt, K. (2003). Comparison of
   the Empirical Bayes and the Significance Analysis of Microarrays.
   \emph{Technical Report}, SFB 475, University of Dortmund, Germany.
   \url{http://www.sfb475.uni-dortmund.de/berichte/tr44-03.pdf}.

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.

}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{SAM-class}},\code{\link{sam}},\code{\link{sam.dstat}},
  \code{\link{sam.snp}}
}  
\examples{\dontrun{
  # Load the package multtest and the data of Golub et al. (1999)
  # contained in multtest.
  library(multtest)
  data(golub)
  
  # Perform a SAM analysis using Wilcoxon rank sum statistics.
  sam.wilc(golub,golub.cl,rand=123)
  
  # Alternative way of performing the same analysis
  sam(golub,golub.cl,method="wilc.stat",rand=123)
  
}}
\keyword{htest}
