\name{delta.plot}
\alias{delta.plot}
\title{Delta Plots}
\description{
  Generates both a plot of \eqn{\Delta}{Delta} vs. the FDR and a plot of 
  \eqn{\Delta}{Delta} vs. the number of identified genes in a SAM analysis.
}
\usage{
  delta.plot(object, delta = NULL, helplines = FALSE)
}
\arguments{
  \item{object}{a object of class SAM}
  \item{delta}{a vector of values for \eqn{\Delta}{Delta}. If \code{NULL},
      a default set of \eqn{\Delta}{Delta} values will be used }
  \item{helplines}{if \code{TRUE}, help lines will be drawn in the
      \eqn{\Delta}{Delta} plots}
}
\details{
  The \eqn{\Delta}{Delta} plots are a visualization of the table generated
  by \code{sam} that contains the estimated FDR and the number of identified
  genes for a set of \eqn{\Delta}{Delta} values. 
}
\value{
  Two plots in one graphsheet: The plot of \eqn{\Delta}{Delta} vs. FDR and
  the plot of \eqn{\Delta}{Delta} vs. the number of identified genes
}

\references{
   Tusher, V., Tibshirani, R., and Chu, G. (2001). Significance Analysis of
   Microarrays Applied to the Ionizing Radiation Response. \emph{PNAS}, 98,
   5116-5121.
}


\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{\code{\link{SAM-class}},\code{\link{sam}}}
\examples{\dontrun{
  # Load the package multtest and the data of Golub et al. (1999)
  # contained in multtest.
  library(multtest)
  data(golub)

  # Perform a SAM analysis.
  sam.out<-sam(golub,golub.cl,B=100,rand=123)
  
  # Generate the Delta plots for the default set of Deltas computed by sam.
  delta.plot(sam.out)
  
  # Another way of generating the same plot.
  plot(sam.out)

  # Generate the Delta plots for Delta = 0.2, 0.4, ..., 2.
  plot(sam.out,seq(0.2,2,0.2))
}}

\keyword{hplot}
