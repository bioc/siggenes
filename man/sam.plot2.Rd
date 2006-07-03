\name{sam.plot2}
\alias{sam.plot2}
\alias{plotArguments}
\title{SAM Plot}
\description{
  Generates a SAM plot for a specified value of Delta.
}
\usage{
  sam.plot2(object, delta, pos.stats = NULL, sig.col = 3, xlim = NULL, 
        ylim = NULL, main = NULL, xlab = NULL, ylab = NULL, pty = "s", 
        lab = c(10, 10, 7), pch = NULL, sig.cex = 1, ...)
     
  plotArguments(pos.stats = NULL, sig.col = 3, xlim = NULL, ylim = NULL,
        main = NULL, xlab = NULL, ylab = NULL, pty = "s", lab = c(10, 10, 7),
        pch = NULL, sig.cex = 1)
}
\arguments{
  \item{object}{an object of class SAM}
  \item{delta}{a numeric value specifying the value of \eqn{\Delta}{Delta}
     for which the SAM plot should be generated}
  \item{pos.stats}{an integer between 0 and 2. If \code{pos.stats=1}, general
     information as the number of significant genes and the estimated FDR for
     the specified value of \code{delta} will be plotted in the upper left corner
     of the plot. If \code{pos.stats=2}, these information will be plotted in
     the lower right corner. If \code{pos.stats=0}, no information will be plotted.
     By default, \code{pos.stats=1} if the expression score \eqn{d} can be both 
     positive and negative, and \code{pos.stats=2} if \eqn{d} can only take positive
     values}
  \item{sig.col}{a specification of the color of the significant genes. If \code{sig.col}
     has length 1, all the points corresponding to significant genes are marked in the
     color specified by \code{sig.col}. If \code{length(sig.col)==2}, the down-regulated
     genes, i.e. the genes with negative expression score \eqn{d}, are marked in the
     color specified by \code{sig.col}[1], and the up-regulated genes, i.e. the genes
     with positive \eqn{d}, are marked in the color specified by \code{sig.col}[2]. 
     For a description of how colors are specified, see \code{?par}}
  \item{xlim}{a numeric vector of length 2 specifying the x limits (minimum and maximum)
     of the plot}
  \item{ylim}{a numeric vector of length 2 specifying the y limits of the plot}
  \item{main}{a character string naming the main title of the plot}
  \item{xlab}{a character string naming the label of the x axis}
  \item{ylab}{a character string naming the label of the y axis}
  \item{pty}{a character specifying the type of plot region to be used. \code{"s"} (default)
     generates a square plotting region, and \code{"m"} the maximal plotting region}
  \item{lab}{a numeric vector of length 3 specifying the approximate number of tickmarks
     on the x axis and on the y axis and the label size}
  \item{pch}{either an integer specifying a symbol or a single character to be used as the
     default in plotting points. For a description of how \code{pch} can be specified, see
     \code{?par}}
  \item{sig.cex}{a numerical value giving the amount by which the symbols of the significant
     genes should be scaled relative to the default}
  \item{\dots}{further graphical parameters. See \code{?plot.default} and \code{?par}}
}
\value{
  a SAM plot
}
\note{ 
   SAM was deveoped by Tusher et al. (2001).
    
   !!! There is a patent pending for the SAM technology at Stanford University. !!!
}

\references{
    Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
  \code{\link{SAM-class}},\code{\link{sam}},\code{\link{sam.dstat}},\code{\link{sam.wilc}},
  \code{\link{sam.snp}}
}
\examples{\dontrun{
  # Load the package multtest and the data of Golub et al. (1999)
  # contained in multtest.
  library(multtest)
  data(golub)
  
  # Perform a SAM analysis for the two class unpaired case assuming
  # unequal variances.
  sam.out<-sam(golub,golub.cl,B=100,rand=123)
  
  # Generate a SAM plot for Delta = 2
  sam.plot2(sam.out,2)
  
  # Alternatively way of generating the same SAM plot
  plot(sam.out,2)
  
}}
\keyword{hplot}
