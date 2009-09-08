\name{md.plot}
\alias{md.plot}
\title{MD Plot}
\description{
  Generates an MD plot for a specified value of Delta.
  
  Contrary to a SAM plot in which the observed values of the test
  statistic \eqn{D} are plotted against the expected ones, the difference \eqn{M} 
  between the observed and the expected values are plotted against
  the observed values in an MD plot. 
}
\usage{
  md.plot(object, delta, pos.stats = 1, sig.col = 3, xlim = NULL, ylim = NULL, 
        main = NULL, xlab = NULL, ylab = NULL, xsym = NULL, ysym = NULL, 
        forceDelta = FALSE, includeZero = TRUE, lab = c(10, 10, 7), pch = NULL, 
        sig.cex = 1, ...)
}
\arguments{
  \item{object}{an object of class SAM.}
  \item{delta}{a numeric value specifying the value of \eqn{\Delta}{Delta}
     for which the SAM plot should be generated.}
  \item{pos.stats}{an integer between 0 and 2. If \code{pos.stats = 1}, general
     information as the number of significant genes and the estimated FDR for
     the specified value of \code{delta} will be plotted in the upper left corner
     of the plot. If \code{pos.stats = 2}, these information will be plotted in
     the lower right corner. If \code{pos.stats = 0}, no information will be plotted.}
  \item{sig.col}{a specification of the color of the significant genes. If \code{sig.col}
     has length 1, all the points corresponding to significant genes are marked in the
     color specified by \code{sig.col}. If \code{length(sig.col) == 2}, the down-regulated
     genes, i.e. the genes with negative expression score \eqn{d}, are marked in the
     color specified by \code{sig.col}[1], and the up-regulated genes, i.e. the genes
     with positive \eqn{d}, are marked in the color specified by \code{sig.col}[2]. 
     For a description of how colors are specified, see \code{\link{par}}.}
  \item{xlim}{a numeric vector of length 2 specifying the x limits (minimum and maximum)
     of the plot.}
  \item{ylim}{a numeric vector of length 2 specifying the y limits of the plot.}
  \item{main}{a character string naming the main title of the plot.}
  \item{xlab}{a character string naming the label of the x axis.}
  \item{ylab}{a character string naming the label of the y axis.}
  \item{xsym}{should the range of the plotted x-axis be symmetric about the origin? Ignored if \code{xlim}
     is specified. If \code{NULL}, \code{xsym} will be set to \code{TRUE}, if some of the observed
     values of the test statistic are negative. Otherwise, \code{xsym} will be set to \code{FALSE}.}
  \item{ysym}{should the range of the plotted y-axis be symmetric about the origin? Ignored if \code{ylim}
     is specified.If \code{NULL}, \code{ysym} will be set to \code{TRUE}, if some of the observed
     values of the test statistic are negative. Otherwise, \code{ysym} will be set to \code{FALSE}.}
  \item{forceDelta}{should the two horizontal lines at \code{delta} and -\code{delta} be
     within the plot region, no matter whether they are out of the range of the observed
     \eqn{d} values? Ignored if \code{ylim} is specified.} 
  \item{includeZero}{should \eqn{D = 0} and \eqn{M = 0} be included in the plot, although all observed
     values of \eqn{D} (or \eqn{M}) are larger than zero?}
  \item{lab}{a numeric vector of length 3 specifying the approximate number of tickmarks
     on the x axis and on the y axis and the label size.}
  \item{pch}{either an integer specifying a symbol or a single character to be used as the
     default in plotting points. For a description of how \code{pch} can be specified, see
     \code{\link{par}}.}
  \item{sig.cex}{a numerical value giving the amount by which the symbols of the significant
     genes should be scaled relative to the default.}
  \item{\dots}{further graphical parameters. See \code{\link{plot.default}} and 
     \code{\link{par}}.}
}
\value{
  A MD plot.
}



\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
  \code{\link{sam}}, \code{\link{sam.plot2}}
}
\examples{\dontrun{
  # Load the package multtest and the data of Golub et al. (1999)
  # contained in multtest.
  library(multtest)
  data(golub)
  
  # Perform a SAM analysis for the two class unpaired case assuming
  # unequal variances.
  sam.out <- sam(golub, golub.cl, B=100, rand=123)
  
  # Generate a SAM plot for Delta = 2
  plot(sam.out, 2)
  
  # As an alternative, the MD plot can be generated.
  md.plot(sam.out, 2)
  
}}
\keyword{hplot}
