\name{wilc.stat}
\alias{wilc.stat}

\title{SAM Analysis Using Wilcoxon Rank Statistics}
\description{
  Generates the required statistics for a Significance Analysis of Microarrays
  analysis using standardized Wilcoxon rank statistics. 
  
  Should not be called directly, but via sam(..., method = wilc.stat).
}
\usage{
   wilc.stat(x, y, gene.names = NULL, R.fold = 1, use.dm = FALSE,
       R.unlog = TRUE, na.replace = TRUE, na.method = "mean", 
       approx50 = TRUE, ties.method=c("min","random","max"), 
       use.row = FALSE, rand = NA)
}
\arguments{
  \item{x}{a matrix or a data frame. Each row of \code{x} must correspond to a variable (e.g., a gene),
    and each column to a sample (i.e.\ an observation).}
  \item{y}{a numeric vector of length \code{ncol(x)} containing the class
     labels of the samples. In the two class paired case, \code{y} can also 
     be a matrix with \code{ncol(x)} rows and 2 columns. For details
     on how \code{y} should be specified, see \code{?sam}.}
  \item{gene.names}{a character vector of length \code{nrow(x)} containing the
     names of the genes.}
  \item{R.fold}{a numeric value. If the fold change of a gene is smaller than or
     equal to \code{R.fold}, or larger than or equal to 1/\code{R.fold},respectively,
     then this gene will be excluded from the SAM analysis. The expression score 
     \eqn{d} of excluded genes is set to \code{NA}. By default, \code{R.fold}
     is set to 1 such that all genes are included in the SAM analysis. Setting 
     \code{R.fold} to 0 or a negative value will avoid the computation of the fold
     change. The fold change is only computed in the two-class unpaired case.}
  \item{use.dm}{if \code{TRUE}, the fold change is computed by 2 to the power of the difference between
     the mean log2 intensities of the two groups, i.e.\ 2 to the power of the numerator of the test statistic. 
     If \code{FALSE}, the fold change is determined
     by computing 2 to the power of \code{x} (if \code{R.unlog = TRUE}) and then calculating the ratio of the
     mean intensity in the group coded by 1 to the mean intensity in the group coded
     by 0. The latter is the default, as this definition of the fold change is used in
     Tusher et al.\ (2001).} 
  \item{R.unlog}{if \code{TRUE}, the anti-log of \code{x} will be used in the computation of the
     fold change. Otherwise, \code{x} is used. This transformation should be done
     if \code{x} is log2-tranformed. (In a SAM analysis, it is highly recommended
     to use log2-transformed expression data.) Ignored if \code{use.dm = TRUE}.}
  \item{na.replace}{if \code{TRUE}, missing values will be removed by the genewise/rowwise
     statistic specified by \code{na.method}. If a gene has less than 2 non-missing
     values, this gene will be excluded from further analysis. If \code{na.replace = FALSE},
     all genes with one or more missing values will be excluded from further analysis.
     The expression score \eqn{d} of excluded genes is set to \code{NA}.}
  \item{na.method}{a character string naming the statistic with which missing values
     will be replaced if \code{na.replace=TRUE}. Must be either \code{"mean"} (default)
     or \code{median}.}
  \item{approx50}{if \code{TRUE}, the null distribution will be approximated by
     the standard normal distribution. Otherwise, the exact null distribution is
     computed. This argument will automatically be set to \code{FALSE} if there
     are less than 50 samples in each of the groups.}
  \item{ties.method}{either \code{"min"} (default), \code{"random"}, or \code{"max"}. If
     \code{"random"}, the ranks of ties are randomly assigned. If \code{"min"} or \code{"max"},
     the ranks of ties are set to the minimum or maximum rank, respectively. For details,
     see the help of \code{\link{rank}}. If \code{use.row = TRUE}, \code{ties.method = "max"}
     will be used. For the handling of Zeros, see Details.}
  \item{use.row}{if \code{TRUE}, \code{\link{rowWilcoxon}} is used to compute the Wilcoxon
     rank statistics.}
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state.}
}
\details{
  Standardized versions of the Wilcoxon rank statistics are computed. This means that
  \eqn{W* = (W - W_{mean}) / W_{sd}}{W* = (W - mean(W)) / sd(W)} is used as expression 
  score \eqn{d}, where \eqn{W} is the usual Wilcoxon rank sum statistic or Wilcoxon
  signed rank statistic, respectively. 
  
  In the computation of these statistics, the ranks of ties are by default set to the
  minimum rank. In the computation of the Wilcoxon signed rank statistic, zeros are randomly 
  set either to a very small positive or negative value.
  
  If there are less than 50 observations in each of the groups, the exact null distribution
  will be used. If there are more than 50 observations in at least one group, the null
  distribution will by default be approximated by the standard normal distribution. It is,
  however, still possible to compute the exact null distribution by  setting \code{approx50}
  to \code{FALSE}. 
}
\value{
   A list containing statistics required by \code{sam}.
}

\references{
   Schwender, H., Krause, A. and Ickstadt, K. (2003). Comparison of
   the Empirical Bayes and the Significance Analysis of Microarrays.
   \emph{Technical Report}, SFB 475, University of Dortmund, Germany.

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.

}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{SAM-class}},\code{\link{sam}, \code{\link{wilc.ebam}}}
}  

\keyword{htest}
