\name{sam.delta}
\alias{sam.delta}
\title{Further SAM analysis}
\description{If the output of a previous analysis with \code{sam} was assigned to a vector, this function
    can be used to compute the number of significant genes and the FDR for values of the threshold Delta
    that was not used in the previous analysis.}

\usage{sam.delta(sam.out,delta)}

\arguments{
    \item{sam.out}{the object to which the output of a previous analysis with \code{sam} was assigned.}
    \item{delta}{vector of values for the threshold Delta.}
}

\value{a table of statistics (estimate for \eqn{p_0}{p0}, number of significant genes,
    number of falsely called genes and FDR) for the specified set of Deltas.
}

\note{
    This function can only be used for a SAM analysis -- not for a SAM.Wilc analysis.

    SAM was introduced by Tusher et al. (2001).

    !!! There is a patent pending for the SAM technology at Stanford University. !!!
}

\seealso{
    \code{\link{sam.plot}} \code{\link{sam}} \code{\link{sam.wilc}}
}


\references{
    Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
    applied to the ionizing radiation response, \emph{PNAS}, 98, 5116-5121.

    Schwender, H. (2003). Assessing the False Discovery Rate in a Statistical Analysis
    of Gene Expression Data, Chapter 5, \emph{Diploma thesis}, Department of Statistics,
    University of Dortmund, \url{http://de.geocities.com/holgerschw/thesis.pdf}.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}



\examples{
    library(multtest)
    # Load the data of Golub et al. (1999). data(golub) contains a 3051x38 gene expression
    # matrix called golub, a vector of length called golub.cl that consists of the 38 class labels,
    # and a matrix called golub.gnames whose third column contains the gene names.
    data(golub)
    
    # Performing a SAM Analysis of the Golub data. Setting rand=123, to make the results reproducible,
    # and setting med.fdr=FALSE, such that the mean number instead of the median number of falsely called
    # genes is computed. The output is assigned to an object for further analyses.
    
    sam.output<-sam(golub,golub.cl,med.fdr=FALSE,rand=123)
    
    # Now the number of significant genes and the FDR should be computed for other values of the threshold
    # Delta. Let's say Delta from 0.6 to 0.8 in steps of length 0.01.
    
    sam.delta(sam.output,seq(0.6,0.8,0.01))
}

\keyword{htest}
