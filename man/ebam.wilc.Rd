\name{ebam.wilc}
\alias{ebam.wilc}
\title{Empirical Bayes Analysis using Wilcoxon Rank Sums}
\description{Performs an Empirical Bayes Analysis of Microarrays by using Wilcoxon Rank Sums
as expression scores for the genes.}
\usage{
    ebam.wilc(data,cl,delta=.9,p0=NA,ties.rand=TRUE,zero.rand=TRUE,gene.names=NULL,
        R.fold=TRUE,R.unlog=TRUE,file.out=NA,na.rm=FALSE,rand=NA)
}

\arguments{
    \item{data}{the data set that should be analyzed. Every row of this data set must
         correspond to a gene, and each column to a biological sample.}
    \item{cl}{a vector containing the class labels of the samples. In the two class unpaired case,
    the label of a sample is either 0 (e.g, control group) or 1 (e.g., case group).
    In the two class paired case, the labels are the integers between 1 and \eqn{n/2}
    (e.g., after treatment group) and between -1 and \eqn{-n/2} (e.g., before treatment
    group), where \eqn{n} is the length of \code{cl} and \eqn{k} is paired with \eqn{-k}.}
    \item{delta}{a gene will be called significant, if its posterior probability of
        being differentially expressed is larger than or equal to \code{delta}.}
    \item{p0}{prior probability that a gene is differentially expressed. If not specified,
        it will automatically be computed.}
    \item{ties.rand}{if \code{TRUE} (default), non-integer expression scores will be randomly
        assigned to the next lower or upper integer. Otherwise, they are assigned to
        the integer that is closer to the mean.}
    \item{zero.rand}{if \code{TRUE} (default), the sign of each Zero in the computation of
        the Wilcoxon signed rank sums will be randomly assigned. If \code{FALSE}, the
        sign of the Zeros will be set to '--'.}
    \item{gene.names}{a vector containing the names of the genes.}
    \item{R.fold}{if \code{TRUE} (default), the fold change for each differentially
        expressed gene will be computed.}
    \item{R.unlog}{if \code{TRUE}, the anti-log of \code{data} will be used in the computation of the
        R.fold. This is recommended if \code{data} consists of log2 transformed gene expression
        data.}
    \item{file.out}{if specified, general information like the number of significant 
        genes and the estimated FDR and gene-specific information like the expression
        scores, the q-values, the R fold etc. of the differentially expressed genes
        are stored in this file.}
    \item{na.rm}{if \code{FALSE} (default), the fold change of genes with at least one
        missing value will be set to \code{NA}. If \code{TRUE}, missing values will be
        replaced by the genewise mean.}
    \item{rand}{if specified, the random number generator will be set in a reproducible state.}
}

\value{a plot of the expression scores vs.\ their posterior probability of being differentially
    expressed, and (optionally) a file containing general information like the FDR and the
    number of differentially expressed genes and gene-specific information on the differentially
    expressed genes like their names, their q-values and their fold change.

    \item{nsig}{number of significant genes.}
    \item{fdr}{estimated FDR.}
    \item{ebam.output}{table containing gene-specific information on the differentially
        expressed genes.}
    \item{row.sig.genes}{vector containing of the row numbers that belong to the differentially
        expressed genes.}
    \item{...}
}

\seealso{
    \code{\link{ebam}}
}

\references{

    Efron, B., Storey, J.D., Tibshirani, R.\ (2001). Microarrays, empirical Bayes methods, and
    the false discovery rate, \emph{Technical Report}, Department of Statistics, Stanford
    University.

    Storey, J.D., and Tibshirani, R. (2003). Statistical significance for genome-wide
    experiments, \emph{Technical Report}, Department of Statistics, Stanford University.

    Schwender, H. (2003). Assessing the false discovery rate in a statistical analysis of
    gene expression data, Chapter 8, \emph{Diploma thesis}, Department of Statistics,
    University of Dortmund, \url{http://de.geocities.com/holgerschw/thesis.pdf}.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\examples{\dontrun{
    library(multtest)
    # Load the data of Golub et al. (1999). data(golub) contains
    # a 3051x38 gene expression matrix called golub, a vector of
    # length called golub.cl that consists of the 38 class labels,
    # and a matrix called golub.gnames whose third column contains
    # the gene names.
    data(golub)
    
    # An EBAM-Wilc analysis of the Golub data is performed by
    
    ebam.wilc.out<-ebam.wilc(golub,golub.cl,gene.names=golub.gnames[,3],rand=123)
    
    # For further analyses, the row numbers of the differentially expressed
    # genes are obtained by
    
    ebam.wilc.out$row.sig.genes
}}

\keyword{htest}
