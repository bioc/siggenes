\name{ebam.wilc}
\alias{ebam.wilc}
\title{Empirical Bayes Analysis using Wilcoxon Rank Sums}
\description{Performs an Empirical Bayes Analysis of Microarrays by using Wilcoxon Rank Sums
as expression scores for the genes.}
\usage{
    ebam.wilc(data,x,y,paired=FALSE,delta=.9,p0=NA,stable.p0=TRUE,use.offset=TRUE,use.weights=TRUE,
    ties.rand=TRUE,zero.rand=TRUE,ns.df=5,col.accession=NA,col.gene.name=NA,R.fold=TRUE,
    R.dataset=data,file.out=NA,rand=NA,na.rm=FALSE)
}

\arguments{
    \item{data}{the data set that should be analyzed. Every row of this data set must
         correspond to a gene.}
    \item{x}{vector of the columns of \code{data} that correspond to the treatment group. In the
        paired case, \eqn{(x[i],y[i])} build a pair. If, e.g., the first \eqn{n_1}{n1} columns of 
        \code{data} build the treatment group, \code{x=1:\eqn{n_1}{n1}}.}
    \item{y}{vector of the columns of \code{data} that correspond to the control group. In the
        paired case, \eqn{(x[i], y[i])} are an observation pair.}
    \item{paired}{paired (\code{TRUE}) or unpaired (\code{FALSE}) data. Default is \code{FALSE}.}
    \item{delta}{a gene will be called significant, if its posterior probability of
        being differentially expressed is larger than or equal to \code{delta}.}
    \item{p0}{prior probability that a gene is differentially expressed. If not specified,
        it will automatically be computed.}
    \item{stable.p0}{if \code{TRUE} (default), \eqn{p_0}{p0} will be computed by the algorithm
        of Storey and Tibshirani (2003). If \code{FALSE}, the (unstable) estimate will
        be computed that ensures that the posterior probability of being differentially
        expressed is always nonnegative.}
    \item{use.offset}{if \code{TRUE} (default), an offset will be used in the Poisson regression
        for the estimation of the density of the expression scores of all genes.}
    \item{use.weights}{if \code{TRUE} (default), weights are used in the natural cubic spline
        fit for the estimation of \eqn{p_0}{p0}.}
    \item{ties.rand}{if \code{TRUE} (default), non-integer expression scores will be randomly
        assigned to the next lower or upper integer. Otherwise, they are assigned to
        the integer that is closer to the mean.}
    \item{zero.rand}{if \code{TRUE} (default), the sign of each Zero in the computation of
        the Wilcoxon signed rank sums will be randomly assigned. If \code{FALSE}, the
        sign of the Zeros will be set to '--'.}
    \item{ns.df}{the number of degrees of freedom used in the Poisson regression for the 
        estimation of the mixture density of the expression scores of all genes.}
    \item{col.accession}{the column of \code{data} containing the accession numbers of
        the genes. If specified, the accession numbers of the significant genes
        will be added to the output.}
    \item{col.gene.name}{the column of \code{data} that contains the names of the genes.
        If specified, the names of the significant genes will be added to the output.}
    \item{R.fold}{if \code{TRUE} (default), the fold change for each differentially
        expressed gene will be computed.}
    \item{R.dataset}{the data set used in the computation of the fold change. This data
        set can be a transformed version of \code{data}.}
    \item{file.out}{if specified, general information like the number of significant 
        genes and the estimated FDR and gene-specific information like the expression
        scores, the q-values, the R fold etc. of the differentially expressed genes
        are stored in this file.}
    \item{rand}{if specified, the random number generator will be set in a reproducible state.}
    \item{na.rm}{if \code{FALSE} (default), the fold change of genes with at least one
        missing value will be set to \code{NA}. If \code{TRUE}, missing values will be
        replaced by the genewise mean.}
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

\keyword{htest}
