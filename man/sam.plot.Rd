\name{sam.plot}
\alias{sam.plot}
\title{SAM Analysis for a specific threshold}
\description{
    This functions performs a SAM or a SAM-Wilc analysis,
    respectively, for a specific threshold Delta, and (optionally) stores information
    on the differentially expressed genes in a file.
}
\usage{
    sam.plot(sam.out,delta,data,q.values=TRUE,R.fold=TRUE,na.rm=FALSE,pty.square=TRUE,file.out=NA,
        col.accession=NA,col.gene.name=NA,use.numbers=sam.out$use.numbers,rand=sam.out$rand)
}
\arguments{
    \item{sam.out}{the object to which the output of a previous analysis with
        \code{sam} or \code{sam.wilc} was stored.}
    \item{delta}{the value of the threshold \eqn{\Delta} for which the analysis should be performed.}
    \item{data}{the used data set. If the fold change should be computed, this data set
        can be a transformed version of the data set used in the previous analysis
        with \code{sam} or \code{sam.wilc}.}
    \item{q.values}{if \code{TRUE} (default), the q-value for each gene will be computed.}
    \item{R.fold}{if \code{TRUE} (default), the fold change for each differentially 
        expressed gene will be computed. If the previous SAM analysis was an one class analysis
        \code{R.fold} is automatically set to \code{FALSE}.}
    \item{na.rm}{if \code{FALSE} (default), the fold change of genes with at least one
         missing value will be set to \code{NA}. If \code{TRUE}, missing values will
         be replaced by the genewise mean.}
    \item{pty.square}{if \code{TRUE} (default), a square SAM Plot will be generated with
         x and y axes having the same range.}
    \item{file.out}{if specified (i.e. not \code{NA}), general information like the number
        of significant genes and the estimated FDR and gene-specific information is
        stored in the corresponding file.}
    \item{col.accession}{the column of \code{data} containing the accession numbers of
        the genes. If specified, the accession numbers of the significant genes
        will be added to the output.}
    \item{col.gene.name}{the column of \code{data} that contains the names of the genes.
        If specified, the names of the significant genes will be added to the output.}
    \item{use.numbers}{if \code{TRUE}, the number of observations that correspond to a point
        in the SAM Plot will be used as symbol for this point. Will only be used, if
        \code{sam.out} contains the results of \code{sam.wilc}.}
    \item{rand}{if specified, the random number generator is set in a reproducible state.
        By default, the same \code{set.seed} is used as in the previous analysis with
        either \code{sam} or \code{sam.wilc}}
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
    \code{\link{sam}} \code{\link{sam.wilc} \code{\link{sam.lambda}}}
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

\keyword{htest}
