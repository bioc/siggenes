\name{plotArguments}
\alias{plotArguments}
\title{Plot Arguments}
\description{
  Utility function for generating a plot of a SAM or an EBAM object in
  an html output.
}
\usage{
  plotArguments(pos.stats = NULL, sig.col = 3, xlim = NULL, ylim = NULL,
        main = NULL, xlab = NULL, ylab = NULL, pty = "s", lab = c(10, 10, 7),
        pch = NULL, sig.cex = 1, stats.cex = 0.8, y.intersp = 1.3)
}
\arguments{
  \item{pos.stats}{an integer between 0 and 2 for a SAM plot, and between 0 and
     4 for an EBAM plot. See \code{help.sam(plot)} or \code{help.ebam(plot)},
     respectively, for how \code{pos.stats} can be specified, and for its default}
  \item{sig.col}{a specification of the color of the significant genes. If \code{sig.col}
     has length 1, all the points corresponding to significant genes are marked in the
     color specified by \code{sig.col}. Only for a SAM plot: If \code{length(sig.col) == 2}, 
     the down-regulated
     genes, i.e. the genes with negative expression score \eqn{d}, are marked in the
     color specified by \code{sig.col}[1], and the up-regulated genes, i.e. the genes
     with positive \eqn{d}, are marked in the color specified by \code{sig.col}[2]. 
     For a description of how colors are specified, see \code{\link{par}}}
  \item{xlim}{a numeric vector of length 2 specifying the x limits (minimum and maximum)
     of the plot}
  \item{ylim}{a numeric vector of length 2 specifying the y limits of the plot}
  \item{main}{a character string naming the main title of the plot}
  \item{xlab}{a character string naming the label of the x axis}
  \item{ylab}{a character string naming the label of the y axis}
  \item{pty}{a character specifying the type of plot region to be used. \code{"s"} (default
     for a SAM plot) generates a square plotting region, and \code{"m"} (default for an
     EBAM plot) the maximal plotting region}
  \item{lab}{a numeric vector of length 3 specifying the approximate number of tickmarks
     on the x axis and on the y axis and the label size}
  \item{pch}{either an integer specifying a symbol or a single character to be used as the
     default in plotting points. For a description of how \code{pch} can be specified, see
     \code{\link{par}}}
  \item{sig.cex}{a numerical value giving the amount by which the symbols of the significant
     genes should be scaled relative to the default}
  \item{stats.cex}{the size of the statistics printed in the plot relative to the default size.
    Only available for an EBAM plot}
  \item{y.intersp}{a numeric value specifying the space between the rows in which the statistics
    are plotted. Only available for an EBAM plot}
}
\value{
  a list required by \code{sam2html} or \code{ebam2html} if \code{addPlot = TRUE}
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
  \code{\link{sam2html}},\code{\link{ebam2html}}
}
\keyword{hplot}
