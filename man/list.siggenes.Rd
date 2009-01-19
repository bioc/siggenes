\name{list.siggenes}
\alias{list.siggenes}
\title{List of the significant genes}
\description{
  Lists the genes called differentially expressed by the SAM or the EBAM analysis for
  a specified value of the threshold \eqn{\Delta}{Delta}.
}
\usage{
  list.siggenes(object, delta, file = "", gene.names = NULL, order = TRUE, 
  text = NULL, append = FALSE)
}
\arguments{
  \item{object}{either a SAM- or an EBAM-object.}
  \item{delta}{a numeric value specifying the threshold \eqn{\Delta}{Delta}
      in the SAM or EBAM analysis. Note that the meaning of \eqn{\Delta}{Delta}
      differs between SAM and EBAM: In SAM, it is a strictly positive value, whereas
      in EBAM it is a probability.}
  \item{file}{a character string naming a file in which the output is stored.
     If \code{""}, the significant genes will be shown in the console.}
  \item{gene.names}{a character vector containing the names of the genes. Needs
     only to be specified, if the gene names were not specified in \code{sam} or
     \code{ebam}, respectively.}
  \item{order}{if \code{TRUE}, the gene names will be ordered by their 
     "significance".}
  \item{text}{a character string specifying the heading of the gene list. By
     default, the header specifies the type of analysis and the used value of
     \eqn{\Delta}{Delta}. To avoid a header, set \code{text = ""}.}
  \item{append}{If \code{TRUE}, the output will be appended to \code{file}. If
     \code{FALSE}, any existing file having the name \code{file} will be destroyed.}
}
\value{
  A list of significant genes either shown in the console or stored in a file.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{sam}}, \code{\link{ebam}}
}
\examples{\dontrun{
  # Load the package multtest and the data of Golub et al. (1999)
  # contained in \pkg{multtest}.
  library(multtest)
  data(golub)

  # Perform a SAM analysis.
  sam.out<-sam(golub, golub.cl, B=100, rand=123)
  
  # List the genes called significant by SAM using Delta = 3.1.
  list.siggenes(sam.out, 3.1, gene.names=golub.gnames[,2])
}}

\keyword{IO}
\keyword{file}
