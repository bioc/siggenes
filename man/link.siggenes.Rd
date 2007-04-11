\name{link.siggenes}
\alias{link.siggenes}
\title{Links for a SAM or an EBAM object}
\description{
  Generates a html page with links to several public repositories for a list of genes
  called differentially expressed when using a specific Delta value in a SAM or an
  EBAM analysis.
}
\usage{
  link.siggenes(object, delta, filename, gene.names = NULL, addDataFrame = TRUE,
        ll = TRUE, refseq = TRUE, symbol = TRUE, omim = TRUE, ug =TRUE, 
        chipname = "", cdfname = NULL, n.digits = 3, title = NULL, bg.col = "white",
        text.col = "black", link.col = "blue", tableborder = 1, new.window = TRUE)
}
\arguments{
  \item{object}{a SAM or an EBAM object}
  \item{delta}{a numerical value specifying the Delta value}
  \item{filename}{character string naming the file in which the output should be
        stored. Must have the suffix ".html"}
  \item{gene.names}{a character vector of the same length as \code{object@d} containing
        the names of the genes. Must only be specified if it is not specified in \code{object},
        i.e. if it has not been specified in \code{sam(data, cl, ...)}}
  \item{addDataFrame}{logical indicating if gene-specific information on the differentially
        expressed genes should be added to the output}
  \item{ll}{logical indicating if Entrez/Locus Links should be added to the output}
  \item{refseq}{logical indicating if RefSeq links should be added to the output}
  \item{symbol}{logical indicating if the gene symbols should be added to the output}
  \item{omim}{logical indicating if OMIM links should be added to the output}
  \item{ug}{logical indicating if UniGene links should be added to the output}
  \item{chipname}{character string specifying the chip type used in the analysis. Must
        be specified as in the meta-data section of Bioconductor (e.g., \code{"hgu133a"}
        for the Affymetrix HG-U133A chip). Need not to be specified if \code{cdfname}
        is specified}
  \item{cdfname}{character string specifying the cdf name of the used chip. Must exactly follow
        the nomenclatur of the Affymetrix chips (e.g., \code{"HG-U133A"} for the
        Affymetrix HG-U133A chip). If specified, links to the Affymetrix webpage for the
        interesting genes will be added to the output}
  \item{n.digits}{integer specifying the number of decimal places used in the output}
  \item{title}{character string naming the title that should be used in the html page}
  \item{bg.col}{specification of the background color of the html page. See \code{?par} for
        how colors can be specified}
  \item{text.col}{specification of the color of the text used in the html page. See \code{?par} for
        how colors can be specified}
  \item{link.col}{specification of the color of the links used in the html file. See \code{?par}
        for how colors can be specified}
  \item{tableborder}{integer specifying the thickness of the border of the table}
  \item{new.window}{logical indicating if the links should be opened in a new window}
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{sam}}, \code{\link{ebam}}, \code{\link{link.genes}}, \code{\link{sam2html}},
   \code{\link{ebam2html}}
}

\keyword{IO}
\keyword{file}
