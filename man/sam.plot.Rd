\name{sam.plot}
\alias{sam.plot}
\title{SAM Analysis for a specific threshold}
\description{
    This functions performs a SAM or a SAM-Wilc analysis,
    respectively, for a specific threshold Delta, and (optionally) stores information
    on the differentially expressed genes in a file.
}
\usage{
    sam.plot(sam.out,delta,q.values=TRUE,R.fold=TRUE,R.unlog=TRUE,
        na.rm=FALSE,file.out=NA,gene.names=NULL)
}
\arguments{
    \item{sam.out}{the object to which the output of a previous analysis with
        \code{sam} or \code{sam.wilc} was stored.}
    \item{delta}{the value of the threshold \eqn{\Delta} for which the analysis should be performed.}
    \item{q.values}{if \code{TRUE} (default), the q-value for each gene will be computed.}
    \item{R.fold}{if \code{TRUE} (default), the fold change for each differentially 
        expressed gene will be computed. If the previous SAM analysis was an one class analysis
        \code{R.fold} is automatically set to \code{FALSE}.}
    \item{R.unlog}{if \code{TRUE}, 2\^data will be used in the computation of the
        R.fold. This is recommended if \code{data} consists of log2 transformed gene expression
        data.}
    \item{na.rm}{if \code{FALSE} (default), the fold change of genes with at least one
        missing value will be set to \code{NA}. If \code{TRUE}, missing values will be
        replaced by the genewise mean.}
    \item{file.out}{if specified (i.e. not \code{NA}), general information like the number
        of significant genes and the estimated FDR and gene-specific information is
        stored in the corresponding file.}
    \item{gene.names}{a vector containing the names of the genes.}
}
\value{
    a SAM Plot for the specified \eqn{\Delta}{Delta} and (optionally) an output file
    containing general information like the number of significant genes and the FDR and
    information about the differentially expressed genes like their names, q-values and
    fold changes.
    
    \item{vec.fdr}{a list containing the number of \eqn{p_0}{p0}, significant genes, the FDR
        etc. for the specified \eqn{\Delta}{Delta}.}
    \item{sam.output}{a table containing gene-specific information about the differentially
        expressed genes like their IDs (i.e. the rows of the data set that contain the 
        expression data of these genes), their expression score, q-values, the fold changes
        etc.}
    \item{row.sig.genes}{vector that consists of the rows of the data set that contain the
        expression data of the differentially expressed genes.}
}

\note{
    SAM was developed by Tusher et al. (2001).
    
    !!! There is a patent pending for the SAM technology at Stanford University. !!!
}

\seealso{
    \code{\link{sam}} \code{\link{sam.wilc}} \code{\link{sam.delta}}
}
\references{
    Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
    applied to the ionizing radiation response, \emph{PNAS}, 98, 5116-5121.

    Storey, J.D., and Tibshirani, R. (2003). Statistical significance for genome-wide
    experiments, \emph{Technical Report}, Department of Statistics, Stanford University.

    Schwender, H. (2003). Assessing the false discovery rate in a statistical analysis of
    gene expression data, Chapters 5 and 6, \emph{Diploma thesis}, Department of Statistics,
    University of Dortmund, \url{http://de.geocities.com/holgerschw/thesis.pdf}.
}

\author{Holger Schwender \email{holger.schw@gmx.de}}

\examples{\dontrun{
    library(multtest)
    # Load the data of Golub et al. (1999). data(golub) contains a 3051x38 gene expression
    # matrix called golub, a vector of length called golub.cl that consists of the 38 class labels,
    # and a matrix called golub.gnames whose third column contains the gene names.
    data(golub)
    
    # Performing a SAM Analysis of the Golub data. Setting rand=123, to make the results reproducible,
    # and setting med.fdr=FALSE, such that the mean number instead of the median number of falsely called
    # genes is computed. The output is assigned to an object for further analyses.
    
    sam.output<-sam(golub,golub.cl,med.fdr=FALSE,rand=123)
    
    # Now sam.plot can be used to obtain the SAM plot and information about the differentially expressed
    # genes for a specified value of delta.
    
    sam.14<-sam.plot(sam.output,1.4,gene.names=golub.gnames[,3])
    
    # sam.plot can also be used in a SAM-Wilc analysis for a specific Delta. 
    
    # For further analyses, the row numbers of the differentially expressed genes can be obtained by
    
    sam14$row.sig.genes
}}


\keyword{htest}
