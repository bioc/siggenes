\name{plotFindArguments}
\alias{plotFindArguments}
\title{Plot Arguments}
\description{
  Utility function for generating a plot of the posterior probabilities in an
  html file when searching for the optimal value of the fudge factor in an 
  EBAM analysis.
}
\usage{
  plotFindArguments(onlyTab = FALSE, logit = TRUE, pos.legend = NULL,
        legend.cex = 0.8, col = NULL, main = NULL, xlab = NULL, ylab = NULL,
        only.a0 = FALSE, lty = 1, lwd = 1, y.intersp = 1.1)
}
\arguments{
  \item{onlyTab}{if \code{TRUE}, then this plot is not generated and only the
     table of the number of differentially expressed genes and the estimated FDR
     for the different values of the fudge factor is shown.}
  \item{logit}{should the posterior probabilities be logit-transformed before
     they are plotted?}
  \item{pos.legend}{an integer between 0 and 4. See \code{help.finda0(plot)}
     for how \code{pos.legend} can be specified, and for its default.}
  \item{legend.cex}{the size of the text in the legend relative to the default size}
  \item{col}{a vector specifying the colors of the lines for the different values
     of the fudge factor. For a description of how colors can be specified, see
     \code{\link{par}}.}
  \item{main}{a character string naming the main title of the plot.}
  \item{xlab}{a character string naming the label of the x axis.}
  \item{ylab}{a character string naming the label of the y axis.}
  \item{only.a0}{if \code{TRUE}, only the values of \eqn{a_0}{a0} are shown in the
     legend. If \code{FALSE}, both the values of \eqn{a_0}{a0} and the corresponding
     number of differentially expressed genes are shown.}
  \item{lty}{a value or vector specifying the line type of the curves. For details,
     see \code{\link{par}}.}
  \item{lwd}{a numeric value specifying the width of the plotted lines. For details,
     see \code{\link{par}}.}
  \item{y.intersp}{a numeric value specifying the space between the rows of the legend.}
}
\value{
  A list required by \code{ebam2html} if \code{findA0} is specified.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
  \code{\link{ebam2html}}
}
\keyword{hplot}
