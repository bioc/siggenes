\name{siggenes2excel}
\alias{sam2excel}
\alias{ebam2excel}
\alias{siggenes2excel}
\title{CSV file of a SAM or an EBAM object}
\description{
  Generates a csv file for either a SAM or an EBAM object for the use in Excel. This csv file can contain 
  general information as the number of differentially expressed genes and the estimated FDR,
  and gene-specific information on the differentially expressed genes.
}
\usage{
 sam2excel(object, delta, file, excel.version=1, n.digits = 3, what = "both", entrez = FALSE,
        chip = "", quote = FALSE)
        
 ebam2excel(object, delta, file, excel.version=1, n.digits = 4, what = "both", entrez = FALSE,
        chip = "", quote = FALSE)
}
\arguments{
  \item{object}{either a SAM or an EBAM object}
  \item{delta}{a numerical value specifying the Delta value}
  \item{file}{character string naming the file in which the output should be
        stored. Must have the suffix ".csv"}
  \item{excel.version}{either \code{1} or \code{2}. If \code{excel.version=1} (default) a csv file
        for the use in an Excel version with American standard settings (\code{sep=","} and 
        \code{dec="."}) will be generated. If \code{excel.version=2} a csv file for the European
        standard setting (\code{sep=";"} and \code{dec=","}) will be generated}
  \item{n.digits}{integer specifying the number of decimal places used in the output}
  \item{what}{either \code{"both"}, \code{"stats"} or \code{"genes"}. If \code{"stats"} general 
        information will be shown. If \code{"genes"} gene-specific information will be given. 
        If \code{"both"} both general and gene-specific information will be shown}
  \item{entrez}{logical indicating if both the Entrez links and the symbols of the genes will be added
        to the output}
  \item{chip}{character string naming the chip type used in this analysis. Must
        be specified as in the meta-data section of Bioconductor (e.g., \code{"hgu133a"}
        for the Affymetrix HG-U133A chip). Only needed if \code{ll=TRUE}. 
        If the argument \code{data} in \code{sam(data, cl, ...)} has been specified 
        by an \code{ExpressionSet} object \code{chip} need not to be specified}
  \item{quote}{logical indicating if character strings and factors should be surrounded by 
        double quotes. For details see \code{\link{write.table}}}
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{sam}}, \code{\link{sam2html}}, \code{\link{ebam}}, \code{\link{ebam2html}}
}

\keyword{IO}
\keyword{file}
