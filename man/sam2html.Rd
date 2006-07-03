\name{sam2html}
\alias{sam2html}
\title{HTML page for a SAM object}
\description{
  Generates a html page for a SAM object. This html page can contain general information
  as the number of differentially expressed genes and the estimated FDR, the SAM plot and
  gene-specific information on the differentially expressed genes.
}
\usage{
 sam2html(sam.out, delta, filename, addStats = TRUE, addPlot = TRUE, addGenes = TRUE,
        ll = TRUE, refseq = TRUE, symbol = TRUE, omim = TRUE, ug = TRUE,
        chipname = "", cdfname = NULL, n.digits = 3, bg.col = "white", text.col = "black",
        link.col = "blue", plotArgs = plotArguments(), bg.plot.adjust = FALSE, 
        plotname = NULL, plotborder = 0, tableborder = 1, new.window = TRUE, ...)
}
\arguments{
  \item{sam.out}{a SAM object}
  \item{delta}{a numerical value specifying the Delta value}
  \item{filename}{character string naming the file in which the output should be
        stored. Must have the suffix ".html"}
  \item{addStats}{logical indicating if general information as the number of differentially
        expressed genes and the estimated FDR should be added to the html page}
  \item{addPlot}{logical indicating if the SAM plot should be added to the html page}
  \item{addGenes}{logical indicating if gene-specific information on the differentially
        expressed genes should be added to the html page}
  \item{ll}{logical indicating if Entrez/Locus Links should be added to the output. Ignored if
        \code{addGenes=FALSE}}
  \item{refseq}{logical indicating if RefSeq links should be added to the output. Ignored
        if \code{addGenes=FALSE}}
  \item{symbol}{logical indicating if the gene symbols should be added to the output.
        Ignored if \code{addGenes=FALSE}}
  \item{omim}{logical indicating if OMIM links should be added to the output. Ignored
        if \code{addGenes=FALSE}}
  \item{ug}{logical indicating if UniGene links should be added to the output. Ignored if
        \code{addGenes=FALSE}}
  \item{chipname}{character string specifying the chip type used in the analysis. Must
        be specified as in the meta-data section of Bioconductor (e.g., \code{"hgu133a"}
        for the Affymetrix HG-U133A chip). Need not to be specified if \code{cdfname}
        is specified. Ignored if \code{addGenes=FALSE}}
  \item{cdfname}{character string specifying the cdf name of the used chip. Must exactly follow
        the nomenclatur of the Affymetrix chips (e.g., \code{"HG-U133A"} for the
        Affymetrix HG-U133A chip). If specified, links to the Affymetrix webpage for the
        interesting genes will be added to the output. Ignored if \code{addGenes=FALSE}}
  \item{n.digits}{integer specifying the number of decimal places used in the output}
  \item{bg.col}{specification of the background color of the html page. See \code{?par} for
        how colors can be specified}
  \item{text.col}{specification of the color of the text used in the html page. See \code{?par} for
        how colors can be specified}
  \item{link.col}{specification of the color of the links used in the html file. See \code{?par}
        for how colors can be specified}
  \item{plotArgs}{further arguments for generating the SAM plot. These are the arguments used
        by \code{sam.plot2} and \code{plot}, respectively. See \code{?plotArgs} for these
        arguments. Ignored if \code{addPlot=FALSE}}
  \item{bg.plot.adjust}{logical indicating if the background color of the SAM plot should be
        the same as the background color of the html page. If \code{FALSE} (default) the 
        background of the plot is white. Ignored if \code{addPlot=FALSE}}
  \item{plotname}{character string naming the file in which the SAM plot is stored. This file
        is needed when the SAM plot should be added to the html page. If not specified the SAM
        plot will be stored as png file in the same folder as the html page. Ignored if \code{addPlot=FALSE}}
  \item{plotborder}{integer specifying the thickness of the border around the plot. By default,
        \code{plotborder=0}, i.e.\ no border is drawn around the plot. Ignored if 
        \code{addPlot=FALSE}} 
  \item{tableborder}{integer specifying the thickness of the border of the table. Ignored if
        \code{addGenes=FALSE}}
  \item{new.window}{logical indicating if the links should be opened in a new window}
  \item{...}{further graphical arguments for the SAM plot. See \code{?plot.default} and \code{?par}.
        Ignored if \code{addPlot=FALSE}}
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{SAM-class}}, \code{\link{sam}}, \code{\link{link.genes}}, \code{\link{link.siggenes}},
   \code{\link{sam2excel}}, \code{\link{plotArguments}}
}

\keyword{IO}
\keyword{file}
