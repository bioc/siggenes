\name{sam.wilc}
\alias{sam.wilc}
\title{SAM Analysis using Wilcoxon Rank Sums}
\description{Performs a Significance Analysis of Microarrays for a set of integer thresholds Delta. Instead
of the modified t statistics, it uses Wilcoxon Rank Sums.}
\usage{
    sam.wilc(data,cl,delta=1:max(abs(W.diff)),na.rm=FALSE,zero.rand=TRUE,
        rand=NA,graphic=TRUE,thres=round(quantile(2:max(abs(W.diff)),(0:3)/3)),
        use.numbers=TRUE)
}
\arguments{
    \item{data}{the data set that should be analyzed. Every row of this data set must
        correspond to a gene.}
    \item{cl}{a vector containing the class labels of the samples. In the two class unpaired case,
        the label of a sample is either 0 (e.g., control group) or 1 (e.g., case group). 
        In the two class paired case, the labels are the integers between 1 and \eqn{n/2}
        (e.g., before treatment group) and between -1 and \eqn{-n/2} (e.g., after treatment
        group), where \eqn{n} is the length of \code{cl} and \eqn{k} is paired with \eqn{-k}.}
    \item{delta}{a vector of integer for which the SAM-Wilc analysis should be performed.}
    \item{na.rm}{if \code{FALSE} (default), the expression scores W of genes with one or
        more missing values will be set to \code{NA}. If \code{TRUE}, the missing values
        will be replaced by the genewise mean of the non-missing values.}
    \item{zero.rand}{if \code{TRUE} (default), the sign of each Zero in the calculation of
        the Wilcoxon signed rank score will be randomly assigned. If \code{FALSE},
        the sign of the Zeroes will be set to '--'.}
    \item{rand}{if specified (i.e. not \code{NA}), the random number generator will be
        put in a reproducible state.}
    \item{graphic}{if \code{TRUE} (default), both the SAM plot and the plots of Delta vs. 
        FDR and Delta vs. number of significant genes are generated. To avoid this
        plotting, set \code{graphic=FALSE}.}
    \item{thres}{a vector of integer values for \eqn{\Delta} for which two lines parallel
        to the 45-degree line are generated.}
    \item{use.numbers}{if \code{TRUE} (default), the symbol for each point in the SAM Plot
        of a SAM-Wilc analysis will be the number of observations that correspond to
        this point.}
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

\examples{
    library(multtest)
    # Load the data of Golub et al. (1999). data(golub) contains a 3051x38 gene expression
    # matrix called golub, a vector of length called golub.cl that consists of the 38 class labels,
    # and a matrix called golub.gnames whose third column contains the gene names.
    data(golub)
    
    # Performing a SAM-Wilc Analysis of the Golub data. Setting rand=123 makes the results reproducible.
    # The output is assigned to an object for further analyses.
    
    sam.output<-sam.wilc(golub,golub.cl,rand=123)
}

\keyword{htest}
