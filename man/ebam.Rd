\name{ebam}
\alias{ebam}
\title{Empirical Bayes Analysis of Microarrays}
\description{
    Performs an Empirical Bayes Analysis of Microarrays for a specified value of the
    fudge factor \eqn{a_0}{a0}. Modified versions of the t statistics are used.
}
\usage{
    ebam(a0.out,a0=NA,p0=NA,delta=NA,local.bin=.1,gene.names=NULL,q.values=TRUE,
    R.fold=TRUE,R.unlog=TRUE,na.rm=FALSE,file.out=NA)
}

\arguments{
    \item{a0.out}{the object to which the output of a previous analysis with \code{find.a0}
        was assigned.}
    \item{a0}{the fudge factor. If \code{NA}, the value suggested by \code{find.a0}
        will be used.}
    \item{p0}{prior probability that a gene is differentially expressed. If not specified
        (i.e. \code{NA}), it will automatically be computed.}
    \item{delta}{a gene will be called differentially expressed, if its posterior
        probability of being differentially expressed is large than or equal to
        \code{delta}. By default, the same \code{delta} is used as in \code{find.a0}.}
    \item{local.bin}{specifies the interval used in the estimation of the local FDR for the
        expression score \eqn{z}. By default, this interval is \eqn{[z-0.1,z+0.1]}.}
    \item{gene.names}{a vector containing the names of the genes}
    \item{q.values}{if \code{TRUE} (default), the q-value for each gene will be computed.}
    \item{R.fold}{if \code{TRUE} (default), the fold change for each differentially
        expressed gene will be computed.}
    \item{R.unlog}{if \code{TRUE}, 2\^data will be used in the computation of
    the R.fold. This is recommend if \code{data} contains the log2 transformed gene
    expression levels.}
    \item{na.rm}{if \code{FALSE} (default), the fold change of genes with at least one
        missing value will be set to \code{NA}. If \code{TRUE}, missing values will be
        replaced by the genewise mean.}
    \item{file.out}{if specified, general information like the number of significant 
        genes and the estimated FDR and gene-specific information like the expression
        scores, the q-values, the R fold etc. of the differentially expressed genes
        are stored in this file.}
}


\value{a plot of the expression scores against their posterior probability of
    being differentially expressed, and (optional) a file containing general information
    like the estimated FDR and the number of differentially expressed genes and 
    gene-specific information about the differentially expressed genes like their names,
    their expression scores, q values and their fold changes.

    \item{FDR}{vector containing the estimated \eqn{p_0}{p0}, the number of significant genes,
        the number of falsely called genes and the estimated FDR.}
    \item{ebam.out}{table containing gene-specific information about the differentially
        expressed genes.}
    \item{row.sig.genes}{vector consisting of the row numbers that belong to the differentially
        expressed genes.}
    \item{...}
}

\note{The number of false positives are computed by \eqn{p_0}{p0} times the number of falsely
    called genes.}

\seealso{
    \code{\link{find.a0}}   \code{\link{ebam.wilc}}
}

\references{
    Efron, B., Tibshirani, R., Storey, J.D., and Tusher, V. (2001). Empirical Bayes Analysis
    of a Microarray Experiment, \emph{JASA}, 96, 1151-1160.}

    Storey, J.D., and Tibshirani, R. (2003). Statistical significance for genome-wide
    experiments, \emph{Technical Report}, Department of Statistics, Stanford University.

    Schwender, H. (2003). Assessing the false discovery rate in a statistical analysis of
    gene expression data, Chapter 7, \emph{Diploma thesis}, Department of Statistics,
    University of Dortmund, \url{http://de.geocities.com/holgerschw/thesis.pdf}.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}


\examples{
    library(multtest)
    # Load the data of Golub et al. (1999). data(golub) contains a 3051x38 gene expression
    # matrix called golub, a vector of length called golub.cl that consists of the 38 class labels,
    # and a matrix called golub.gnames whose third column contains the gene names.
    data(golub)
    
    # The optimal value for the fudge factor a0 is computed, where possible values of the a0 are
    # 0 and the 0, 0.05 and 0.1 quantile of the standard deviations of the genes. Setting rand=123
    # makes the results reproducible.
    
    find.out<-find.a0(golub,golub.cl,alpha=c(0,0.05,0.1),rand=123)
    
    # Now that we have find the optimal value of a0, an empirical Bayes analysis can be performed.
    
    ebam.out<-ebam(find.out,gene.names=golub.gnames[,3])
    
    # For further analyses the row numbers of the differentially expressed genes are obtained by
    
    ebam.out$row.sig.genes
}

\keyword{htest}
