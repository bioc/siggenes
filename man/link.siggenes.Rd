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
        entrez = TRUE, refseq = TRUE, symbol = TRUE, omim = FALSE, ug = FALSE,
        fullname = FALSE, which.refseq = "NM", chipname = "", cdfname = NULL, 
        refsnp = NULL, max.associated = 2, n.digits = 3, title = NULL, 
        bg.col = "white", text.col = "black", link.col = "blue", tableborder = 1, 
        new.window = TRUE, load = TRUE)
}
\arguments{
  \item{object}{a SAM or an EBAM object.}
  \item{delta}{a numerical value specifying the Delta value.}
  \item{filename}{character string naming the file in which the output should be
        stored. Must have the suffix ".html".}
  \item{gene.names}{a character vector of the same length as \code{object@d} (or \code{object@z}) containing
        the names of the genes. Must only be specified if it is not specified in \code{object},
        i.e. if it has not been specified in \code{sam} (or \code{ebam}).}
  \item{addDataFrame}{logical indicating if gene-specific information on the differentially
        expressed genes should be added to the output.}
  \item{entrez}{logical indicating if Entrez links should be added to the output.}
  \item{refseq}{logical indicating if RefSeq links should be added to the output.}
  \item{symbol}{logical indicating if the gene symbols should be added to the output.}
  \item{omim}{logical indicating if OMIM links should be added to the output.}
  \item{ug}{logical indicating if UniGene links should be added to the output.}
  \item{fullname}{logical indicating whether the full gene names should be added to the output.}
  \item{which.refseq}{character string or vector naming the first two letters of the RefSeq
        links that should be displayed in the html file.}
  \item{chipname}{character string specifying the chip type used in the analysis. Must
        be specified as in the meta-data section of Bioconductor (e.g., \code{"hgu133a"}
        for the Affymetrix HG-U133A chip). Needs not to be specified if \code{cdfname}
        is specified. For Affymetrix SNP chips (starting with the 500k array set), 
        \code{chipname} can be specified by the metadata package name, i.e.\ either by
        \code{"pd.genomewidesnp.5"}, by \code{"pd.genomewidesnp.6"}, by \code{"pd.mapping250k.nsp"},
        or by \code{"pd.mapping250k.sty"}, to add links to the Affymetrix webpage of the SNPs
        to the html output.}
  \item{cdfname}{character string specifying the cdf name of the used chip. Must exactly follow
        the nomenclatur of the Affymetrix chips (e.g., \code{"HG-U133A"} for the
        Affymetrix HG-U133A chip). If specified, links to the Affymetrix webpage for the
        interesting genes will be added to the output. If SNP chips are considered, \code{chipname}
        instead of \code{cdfname} must be specified for obtaining these links.}
  \item{refsnp}{either a character vector or a data frame. If the former, \code{refsnp} containis
          the RefSNP IDs of the SNPs used in the SAM/EBAM analysis, where \code{names(refsnp)} specifies
          the names of these SNPs, i.e.\ their probe set IDs. If a data frame, then one column of \code{refsnp} must contain
          the RefSNP IDs of the SNPs, and the name of this column must be \code{RefSNP}. The other
          columns can contain additional annotations such as the chromosome or the physical position
          of each SNPs. The row names of \code{refsnp} must specify the SNPs, i.e.\ must be the
          probe set IDs of the SNPs. Using \code{buildSNPannotation} from the package \pkg{scrime}
          such a data frame can be
          generated automatically from the metadata package corresponding to the considered SNP chip.} 
  \item{max.associated}{integer specifying the maximum number of genes associated with the respective
          SNP displayed in the html output. If all entries should be shown, set \code{max.associated = 0}.
          This however might result in a very large html output. For details, see 
          \code{shortenGeneDescription} in the package \pkg{scrime}.}
  \item{n.digits}{integer specifying the number of decimal places used in the output.}
  \item{title}{character string naming the title that should be used in the html page.}
  \item{bg.col}{specification of the background color of the html page. See \code{?par} for
        how colors can be specified.}
  \item{text.col}{specification of the color of the text used in the html page. See \code{?par} for
        how colors can be specified.}
  \item{link.col}{specification of the color of the links used in the html file. See \code{?par}
        for how colors can be specified.}
  \item{tableborder}{integer specifying the thickness of the border of the table.}
  \item{new.window}{logical indicating if the links should be opened in a new window.}
  \item{load}{logical value indicating whether to attempt to load the required annotation data package
        if it is not already loaded. For details, see the man page of \code{lookUp} in the package
        \pkg{annotate}.}
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{sam}}, \code{\link{ebam}}, \code{\link{link.genes}}, \code{\link{sam2html}},
   \code{\link{ebam2html}}
}

\keyword{IO}
\keyword{file}
