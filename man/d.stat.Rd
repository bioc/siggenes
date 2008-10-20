\name{d.stat}
\alias{d.stat}
\alias{sam.dstat}
\title{SAM Analysis Using a Modified t-statistic}
\description{
  Computes the required statistics for a Significance Analysis of Microarrays (SAM)
  using either a (modified) t- or F-statistic.
  
  Should not be called directly, but via the function sam.
}
\usage{
  d.stat(data, cl, var.equal = FALSE, B = 100, med = FALSE, s0 = NA, 
      s.alpha = seq(0, 1, 0.05), include.zero = TRUE, n.subset = 10, 
      mat.samp = NULL, B.more = 0.1, B.max = 30000, gene.names = NULL,
      R.fold = 1, use.dm = TRUE, R.unlog = TRUE, na.replace = TRUE, 
      na.method = "mean", rand = NA)
}
\arguments{
  \item{data}{a matrix, data frame or \code{ExpressionSet} object. Each row of
    \code{data} (or \code{exprs(data)}, respectively) must correspond to a gene,
    and each column to a sample}
  \item{cl}{a numeric vector of length \code{ncol(data)} containing the class
     labels of the samples. In the two class paired case, \code{cl} can also 
     be a matrix with \code{ncol(data)} rows and 2 columns. If \code{data} is
     an \code{ExpressionSet} object, \code{cl} can also be a character string.
     For details on how \code{cl} should be specified, see \code{?sam}}
  \item{var.equal}{if \code{FALSE} (default), Welch's t-statistic will be computed.
     If \code{TRUE}, the pooled variance will be used in the computation of
     the t-statistic}
  \item{B}{numeric value indicating how many permutations should be used in
     the estimation of the null distribution}
  \item{med}{if \code{FALSE} (default), the mean number of falsely called genes
     will be computed. Otherwise, the median number is calculated}
  \item{s0}{a numeric value specifying the fudge factor. If \code{NA} (default),
     \code{s0} will be computed automatically}
  \item{s.alpha}{a numeric vector or value specifying the quantiles of the
     standard deviations of the genes used in the computation of \code{s0}. If
     \code{s.alpha} is a vector, the fudge factor is computed as proposed by
     Tusher et al. (2001). Otherwise, the quantile of the standard deviations
     specified by \code{s.alpha} is used as fudge factor}
  \item{include.zero}{if \code{TRUE}, \code{s0}=0 will also be a possible choice
     for the fudge factor. Hence, the usual t-statistic or F statistic, respectively,
     can also be a possible choice for the expression score \eqn{d}. If \code{FALSE},
     \code{s0=0} will not be a possible choice for the fudge factor. The latter
     follows Tusher et al. (2001) definition of the fudge factor in which only strictly
     positive values are considered}
  \item{n.subset}{a numeric value indicating how many permutations are considered
     simultaneously when computing the p-value and the number of falsely called
     genes. If \code{med = TRUE}, \code{n.subset} will be set to 1}
  \item{mat.samp}{a matrix having \code{ncol(data)} columns except for the two class
     paired case in which \code{mat.samp} has \code{ncol(data)}/2 columns. 
     Each row specifies one permutation of the group labels used in the computation
     of the expected expression scores \eqn{\bar{d}}{d.bar}. If not specified 
     (\code{mat.samp=NULL}), a matrix having \code{B} rows and \code{ncol(data)} is
     generated automatically and used in the computation of \eqn{\bar{d}}{d.bar}. In
     the two class unpaired case and the multiclass case, each row of \code{mat.samp} 
     must contain the same group labels as \code{cl}. In the one class and the two 
     class paired case, each row must contain -1's and 1's. In the one class case, 
     the expression values are multiplied by these -1's and 1's. In the two class paired
     case, each column corresponds to one observation pair whose difference is multiplied
     by either -1 or 1. For more details and examples, see the manual of siggenes}
  \item{B.more}{a numeric value. If the number of all possible permutations is smaller
     than or equal to (1+\code{B.more})*\code{B}, full permutation will be done. 
     Otherwise, \code{B} permutations are used. This avoids that \code{B} permutations
     will be used -- and not all permutations -- if the number of all possible permutations
     is just a little larger than \code{B}}
  \item{gene.names}{a character vector of length \code{nrow(data)} containing the
     names of the genes}
  \item{B.max}{a numeric value. If the number of all possible permutations is smaller
     than or equal to \code{B.max}, \code{B} randomly selected permutations will be used
     in the computation of the null distribution. Otherwise, \code{B} random draws
     of the group labels are used. In the latter way of permuting it is possible that
     some of the permutations are used more than once}
  \item{R.fold}{a numeric value. If the fold change of a gene is smaller than or
     equal to \code{R.fold}, or larger than or equal to 1/\code{R.fold},respectively,
     then this gene will be excluded from the SAM analysis. The expression score 
     \eqn{d} of excluded genes is set to \code{NA}. By default, \code{R.fold}
     is set to 1 such that all genes are included in the SAM analysis. Setting 
     \code{R.fold} to 0 or a negative value will avoid the computation of the fold
     change. The fold change is only computed in the two-class unpaired cases}
  \item{use.dm}{if \code{TRUE}, the fold change is computed by 2 to the power of the difference between
     the mean log2 intensities of the two groups, i.e.\ 2 to the power of the numerator of the test statistic. 
     If \code{FALSE}, the fold change is determined
     by computing 2 to the power of \code{data} (if \code{R.unlog = TRUE}) and then calculating the ratio of the
     mean intensity in the group coded by 1 to the mean intensity in the group coded
     by 0. The latter is the definition of the fold change used in Tusher et al.\ (2001)} 
  \item{R.unlog}{if \code{TRUE}, the anti-log of \code{data} will be used in the computation of the
     fold change. Otherwise, \code{data} is used. This transformation should be done
     when \code{data} is log2-tranformed (in a SAM analysis it is highly recommended
     to use log2-transformed expression data). Ignored if \code{use.dm = TRUE}}
  \item{na.replace}{if \code{TRUE}, missing values will be removed by the genewise/rowwise
     statistic specified by \code{na.method}. If a gene has less than 2 non-missing
     values, this gene will be excluded from further analysis. If \code{na.replace=FALSE},
     all genes with one or more missing values will be excluded from further analysis.
     The expression score \eqn{d} of excluded genes is set to \code{NA}}
  \item{na.method}{a character string naming the statistic with which missing values
     will be replaced if \code{na.replace=TRUE}. Must be either \code{"mean"} (default)
     or \code{median}}
  \item{rand}{numeric value. If specified, i.e. not \code{NA}, the random number generator
     will be set into a reproducible state}
}

\value{
   an object of class SAM
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
  \code{\link{SAM-class}},\code{\link{sam}}
}

\keyword{htest}
