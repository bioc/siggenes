\name{siggenes2html}
\alias{sam2html}
\alias{ebam2html}
\alias{siggenes2html}
\title{HTML page for a SAM or an EBAM object}
\description{
  Generates a html page for a SAM or an EBAM object. This html page can contain general information
  as the number of differentially expressed genes and the estimated FDR, the SAM or EBAM plot, and
  gene-specific information on the differentially expressed genes.
}
\usage{
 ebam2html(object, delta, filename, addStats = TRUE, addPlot = TRUE, 
        addGenes = TRUE, findA0 = NULL, varName = NULL, entrez = TRUE, 
        refseq = TRUE, symbol = TRUE, omim = FALSE, ug = FALSE, 
        chipname = "", cdfname = NULL, which.refseq = "NM", refsnp = NULL, 
        n.digits = 3, bg.col = "white", text.col = "black", link.col = "blue", 
        plotArgs = plotArguments(), plotFindArgs = plotFindArguments(), 
        bg.plot.adjust = FALSE, plotname = NULL, plotborder = 0, 
        tableborder = 1, new.window = TRUE, load = TRUE, ...)

 sam2html(object, delta, filename, addStats = TRUE, addPlot = TRUE, 
        addGenes = TRUE, varName = NULL, entrez = TRUE, refseq = TRUE, 
        symbol = TRUE, omim = FALSE, ug = FALSE, chipname = "", 
        cdfname = NULL, which.refseq = "NM", refsnp = NULL, n.digits = 3, 
        bg.col = "white", text.col = "black", link.col = "blue", 
        plotArgs = plotArguments(), bg.plot.adjust = FALSE, 
        plotname = NULL, plotborder = 0, tableborder = 1, 
        new.window = TRUE, load = TRUE, ...)
}
\arguments{
  \item{object}{a SAM or an EBAM object}
  \item{delta}{a numerical value specifying the Delta value}
  \item{filename}{character string naming the file in which the output should be
        stored. Must have the suffix ".html"}
  \item{addStats}{logical indicating if general information as the number of differentially
        expressed genes and the estimated FDR should be added to the html page}
  \item{addPlot}{logical indicating if the SAM/EBAM plot should be added to the html page}
  \item{addGenes}{logical indicating if gene-specific information on the differentially
        expressed genes should be added to the html page}
  \item{findA0}{an object of class FindA0. If specified, the numbers of differentially expressed
        genes and the estimated FDRs for the different possible values of the fudge factor
        and the corresponding plot of the logit-transformed posterior probabilities are
        included in the html file}
  \item{varName}{character string indicating how the variables should be named. If \code{NULL},
        the variables will be referred to as SNPs in the output if 
        \code{method = cat.stat}, and as Genes otherwise}
  \item{entrez}{logical indicating if Entrez links should be added to the output. Ignored if
        \code{addGenes = FALSE}}
  \item{refseq}{logical indicating if RefSeq links should be added to the output. Ignored
        if \code{addGenes = FALSE}}
  \item{symbol}{logical indicating if the gene symbols should be added to the output.
        Ignored if \code{addGenes = FALSE}}
  \item{omim}{logical indicating if OMIM links should be added to the output. Ignored
        if \code{addGenes = FALSE}}
  \item{ug}{logical indicating if UniGene links should be added to the output. Ignored if
        \code{addGenes = FALSE}}
  \item{chipname}{character string specifying the chip type used in the analysis. Must
        be specified as in the meta-data section of Bioconductor (e.g., \code{"hgu133a"}
        for the Affymetrix HG-U133A chip). Need not to be specified if \code{cdfname}
        is specified. Ignored if \code{addGenes = FALSE}}
  \item{cdfname}{character string specifying the cdf name of the used chip. Must exactly follow
        the nomenclatur of the Affymetrix chips (e.g., \code{"HG-U133A"} for the
        Affymetrix HG-U133A chip). If specified, links to the Affymetrix webpage for the
        interesting genes will be added to the output. Ignored if \code{addGenes = FALSE}}
  \item{which.refseq}{character string or vector naming the first two letters of the RefSeq
        links that should be displayed in the html file}
  \item{refsnp}{character vector containing the RefSNP IDs of the SNPs used in the SAM/EBAM
        analysis, where \code{names(refsnp)} specifies the names of these SNPs} 
  \item{n.digits}{integer specifying the number of decimal places used in the output}
  \item{bg.col}{specification of the background color of the html page. See \code{\link{par}} for
        how colors can be specified}
  \item{text.col}{specification of the color of the text used in the html page. See 
        \code{\link{par}} for how colors can be specified}
  \item{link.col}{specification of the color of the links used in the html file. 
        See \code{\link{par}} for how colors can be specified}
  \item{plotArgs}{further arguments for generating the SAM/EBAM plot. These are the arguments used
        by the SAM/EBAM specific \code{plot} method. See the help of \code{\link{plotArguments}} 
        for these arguments. Ignored if \code{addPlot = FALSE}}
  \item{plotFindArgs}{further arguments for generating the (logit-transformed) posterior
        probabilities for the different values of the fudge factor. Ignored if \code{findA0 = NULL}.
        See the help of \code{\link{plotFindArguments}} for these arguments}
  \item{bg.plot.adjust}{logical indicating if the background color of the SAM plot should be
        the same as the background color of the html page. If \code{FALSE} (default) the 
        background of the plot is white. Ignored if \code{addPlot = FALSE}}
  \item{plotname}{character string naming the file in which the SAM/EBAM plot is stored. This file
        is needed when the SAM/EBAM plot should be added to the html page. If not specified the SAM/EBAM
        plot will be stored as png file in the same folder as the html page. Ignored if \code{addPlot = FALSE}}
  \item{plotborder}{integer specifying the thickness of the border around the plot. By default,
        \code{plotborder = 0}, i.e.\ no border is drawn around the plot. Ignored if 
        \code{addPlot = FALSE}} 
  \item{tableborder}{integer specifying the thickness of the border of the table. Ignored if
        \code{addGenes = FALSE}}
  \item{new.window}{logical indicating if the links should be opened in a new window}
  \item{load}{logical value indicating whether to attempt to load the required annotation data package
        if it is not already loaded. For details, see the man page of \code{lookUp} in the package
        \pkg{annotate}}
  \item{...}{further graphical arguments for the SAM/EBAM plot. See \code{\link{plot.default}} and 
        \code{\link{par}}. Ignored if \code{addPlot = FALSE}}
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{SAM-class}}, \code{\link{sam}}, \code{\link{EBAM-class}}, \code{\link{ebam}},
   \code{\link{link.genes}}, \code{\link{link.siggenes}}, \code{\link{plotArguments}},
   \code{\link{plotFindArguments}}
}

\keyword{IO}
\keyword{file}
