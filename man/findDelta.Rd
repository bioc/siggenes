\name{findDelta}
\alias{findDelta}

\title{Finding the Threshold Delta}
\description{
  Computes the value of the threshold Delta for a given FDR or number of genes/variables
  in a SAM or EBAM analysis.
}
\usage{
   findDelta(object, fdr = NULL, genes = NULL, prec = 6, initial = NULL,
      verbose = FALSE)
}

\arguments{
  \item{object}{either a SAM or an EBAM object.}
  \item{fdr}{numeric value between 0 and 1 for which the threshold Delta and thus the
     number of genes/variables should be obtained. Only one of \code{fdr} and \code{genes}
     can be specified.}
  \item{genes}{integer specifying the number of genes/variables for which the threshold
     Delta and thus the estimated FDR should be obtained. Only one of \code{fdr} and
     \code{genes} can be specified.}
  \item{prec}{integer indicating the precision of the considered Delta values.}
  \item{initial}{a numeric vector of length two containing the minimum and the maximum value of
     Delta that is initially used in the search for Delta. Both values must be larger than 0.
     If \code{object} is an EBAM object, both values must also be smaller than or equal to 1.
     If not specified, the minimum is set to 0.1, and the maximum to either the maximum
     posterior (EBAM) or the maximum absolute distance between the observed and the corresponding
     expected values of the test statistic (SAM).}
  \item{verbose}{should more information about the search process be shown?}
}

\value{
  If a value of Delta is found for the exact value of \code{fdr} or \code{genes}, then a vector
  of length 3 consisting of Delta and the corresponding number of genes and the estimated FDR.
  If such a value is not found, then a matrix with two rows and three columns, where the two
  rows contain the number of genes/variables and the estimated FDR for the two considered values
  of Delta that provide the closest upper and lower bounds to the desired FDR (if \code{fdr} is
  specified) or number of genes/variables (if \code{genes} is specified.)
}


\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{sam}}, \code{\link{ebam}}}

\keyword{htest}
