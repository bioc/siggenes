\name{sam.wilc}
\alias{sam.wilc}
\title{SAM Analysis using Wilcoxon Rank Sums}
\description{Performs a Significance Analysis of Microarrays for a set of integer thresholds Delta. Instead
of the modified t statistics, it uses Wilcoxon Rank Sums.}
\usage{
    sam.wilc(data,x,y,paired=FALSE,na.rm=FALSE,zero.rand=TRUE,rand=NA,use.weights=TRUE,
    delta=1:max(abs(W.diff)),graphic=TRUE,pty.square=TRUE,
    thres=round(quantile(2:max(abs(W.diff)),(0:3)/3)),use.numbers=TRUE,
    helplines=TRUE)
}
\arguments{
    \item{data}{the data set that should be analyzed. Every row of this data set must
        correspond to a gene.}
    \item{x}{vector of the columns of the data set that correspond to the treatment group.
        In the paired case, \eqn{(x[i],y[i])} build a pair. If, e.g., the first \eqn{n_1}{n1} columns
        of \code{data} contain the gene expression values of the treatment group, \code{x=1:n1}.}
    \item{y}{vector of the columns of the data set that correspond to the control group. In the
        paired case, \eqn{(x[i], y[i])} are an observation pair.}
    \item{paired}{paired (\code{TRUE}) or unpaired (\code{FALSE}) data. Default is \code{FALSE}.}
    \item{na.rm}{if \code{FALSE} (default), the expression scores W of genes with one or
        more missing values will be set to \code{NA}. If \code{TRUE}, the missing values
        will be replaced by the genewise mean of the non-missing values.}
    \item{zero.rand}{if \code{TRUE} (default), the sign of each Zero in the calculation of
        the Wilcoxon signed rank score will be randomly assigned. If \code{FALSE},
        the sign of the Zeroes will be set to '--'.}
    \item{rand}{if specified (i.e. not \code{NA}), the random number generator will be
        put in a reproducible state.}
    \item{use.weights}{if \code{T} (default), then the data points are weighed by \eqn{1-\lambda}
        in the fit of a natural cubic spline.}
    \item{delta}{a vector of integer for which the SAM-Wilc analysis should be performed.}
    \item{graphic}{if \code{TRUE} (default), both the SAM plot and the plots of Delta vs. 
        FDR and Delta vs. number of significant genes are generated. To avoid this
        plotting, set \code{graphic=FALSE}.}
    \item{pty.square}{if \code{TRUE} (default), a square SAM plot will be generated with
        x and y axes having the same range.}
    \item{thres}{a vector of integer values for \eqn{\Delta} for which two lines parallel
        to the 45-degree line are generated.}
    \item{use.numbers}{if \code{TRUE} (default), the symbol for each point in the SAM Plot
        of a SAM-Wilc analysis will be the number of observations that correspond to
        this point.}
    \item{helplines}{if \code{TRUE} (default), help lines will be generated in both the
        Delta vs.\ FDR and the Delta vs.\ number of significant genes
        plot.}
}

\value{a table of statistics (estimate for \eqn{p_0}{p0}, number of significant genes,
    number of falsely called genes and FDR) for the specified set of Deltas, a
    SAM Plot, a Delta vs.\ FDR plot, and a plot of Delta vs.\ the number of significant genes.
}

\note{For further analyses with \code{sam.plot}, the results of \code{sam.wilc} must be assigned
    to an object.
    
    SAM was developed by Tusher et al. (2001).
    
    !!! There is a patent pending for the SAM technology at Stanford University. !!!
}

\seealso{
    \code{\link{sam.plot}} \code{\link{sam}}
}

\references{
    Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
    applied to the ionizing radiation response, \emph{PNAS}, 98, 5116-5121.
    
    Schwender, H. (2003). Assessing the False Discovery Rate in a Statistical Analysis
    of Gene Expression Data, Chapter 6, \emph{Diploma thesis}, Department of Statistics,
    University of Dortmund, \url{http://de.geocities.com/holgerschw/thesis.pdf}.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\keyword{htest}
