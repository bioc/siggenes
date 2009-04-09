\name{limma2sam}
\alias{limma2sam}
\alias{limma2ebam}

\title{limma to SAM or EBAM}
\description{
 Transforms the output of an analysis with \pkg{limma} into a \code{SAM}
 or \code{EBAM} object, such that a SAM or EBAM analysis, respectively,
 can be performed using the test statistics provided by \pkg{limma}.
}

\usage{
limma2sam(fit, coef, moderate = TRUE, sam.control = samControl())

limma2ebam(fit, coef, moderate = TRUE, delta = 0.9, 
   ebam.control = ebamControl())
}

\arguments{
  \item{fit}{an object of class \code{MArrayLM}, i.e.\ the output of the functions
    \code{eBayes} and \code{lmFit} from the \pkg{limma} package.}
  \item{coef}{column number or name corresponding to the coefficient or contrast of
    interest. For details, see the argument \code{coef} of the function \code{topTable}
    in \pkg{limma}.}
  \item{moderate}{should the \pkg{limma} t-statistic be considered? If \code{FALSE}, the
    ordinary t-statistic is used in the trasnsformation to a \code{SAM} or \code{EBAM}
    object. If \code{TRUE}, it is expected that \code{fit} is the output of \code{eBayes}.
    Otherwise, \code{fit} can be the result of \code{lmFit} or \code{eBayes}.}
  \item{sam.control}{further arguments for the SAM analysis. See \code{\link{samControl}}
    for these arguments, which should only be changed if they are fully understood.}
  \item{delta}{the minimum posterior probability for a gene to be called differentially
    expressed (or more generally, for a variable to be called significant) in an EBAM
    analysis. For details, see \code{\link{ebam}}. Please note that the meaning of 
    \code{delta} differs substantially between \code{\link{sam}} and \code{\link{ebam}}}
  \item{ebam.control}{further arguments for an EBAM analysis. See \code{\link{ebamControl}}
    for these arguments, which should only be changed if their meaning is fully understood.}
}

\value{
  An object of class \code{SAM} or \code{EBAM}.
}
\references{
   Efron, B., Tibshirani, R., Storey, J.D. and Tusher, V. (2001). Empirical Bayes Analysis
   of a Microarray Experiment. \emph{JASA}, 96, 1151-1160.
   
   Smyth, G.K. (2004). Linear Models and Empirical Bayes Methods for Assessing Differential 
   Expression in Microarray Experiments. \emph{Statistical Applications in Genetics and
   Molecular Biology}, 3(1), Article 3.     

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance Analysis of Microarrays
   Applied to the Ionizing Radiation Response. \emph{PNAS}, 98, 5116-5121.

}
\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{sam}}, \code{\link{ebam}}, \code{\link{SAM-class}}, \code{\link{EBAM-class}},
   \code{\link{samControl}}, \code{\link{ebamControl}}}

\keyword{htest}
