\name{sam}
\alias{sam}
\title{Significance Analysis of Microarrays}
\description{Performs a Significance Analysis of Microarrays (SAM) for
    a set of positive thresholds. It is possible to do either an one class or a two class
    SAM analysis.}
\usage{sam(data,cl,B=100,balanced=FALSE,mat.samp=NULL,delta=(1:10)/5,med.fdr=TRUE,
        s0=NA,alpha.s0=seq(0,1,.05),include.s0=TRUE,p0=NA,lambda.p0=1,vec.lambda.p0=(0:95)/100,
        na.rm=FALSE,graphic.fdr=TRUE,thres.fdr=seq(0.5,2,0.5),ngenes=NA,iteration=3,
        initial.delta=c(.1,seq(.2,2,.2),4),rand=NA)}
\arguments{
    \item{data}{the data set that should be analyzed. Every row of this data set
        must correspond to a gene.}
    \item{cl}{a vector containing the class labels of the samples. In the two class unpaired case,
        the label of a sample is either 0 (e.g., control group) or 1 (e.g., case group). 
        In the two class paired case, the labels are the integers between 1 and \eqn{n/2}
        (e.g., before treatment group) and between -1 and \eqn{-n/2} (e.g., after treatment
        group), where \eqn{n} is the length of \code{cl} and \eqn{k} is paired with \eqn{-k}.
        For one group data, the label for each sample should be 1.}
    \item{B}{number of permutations used in the calculation of the null density.
        Default is \code{B=100}.}
    \item{balanced}{if \code{TRUE}, balanced permutations will be used. Default is \code{FALSE}.}
    \item{mat.samp}{a permutation matrix. If specified, this matrix will be used,
        even if \code{rand} and \code{B} are specified.}
    \item{delta}{a vector of values for the threshold \eqn{\Delta}{Delta} for which the SAM
         analysis is performed.}
    \item{med.fdr}{if \code{TRUE} (default), the median number, otherwise the expected
         number, of falsely called genes will be computed.}
    \item{s0}{the fudge factor. If \code{NA} (default), the fudge factor \eqn{s_0}{s0} will be computed
         automatically.}
    \item{alpha.s0}{the possible values of the fudge factor \eqn{s_0}{s0} in terms of quantiles of the
        standard deviations of the genes.}
    \item{include.s0}{if \code{TRUE} (default), \eqn{s_0=0}{s0=0} is a possible choice for the
        fudge factor.}
    \item{p0}{the probability that a gene is not differentially expressed. If not
        specified (default), it will be computed.}
    \item{lambda.p0}{number between 0 and 1 that is used to estimate \eqn{p_0}{p0}. 
        If set to \code{1} (default), the automatic \eqn{p_0}{p0} selection using 
        the natural cubic spline fit is used.}
    \item{vec.lambda.p0}{vector of values for \eqn{\lambda} used in the automatic
         computation of \eqn{p_0}{p0}.}
    \item{na.rm}{if \code{FALSE} (default), the expression scores \eqn{d} of genes with one or more
         missing values will be set to \code{NA}. If \code{TRUE}, the missing
         values will be replaced by the genewise mean of the non-missing values.}
    \item{graphic.fdr}{if \code{TRUE} (default), both the SAM plot and the plots of Delta vs.
        FDR and Delta vs. number of significant genes will be generated.}
    \item{thres.fdr}{for each value contained in \code{thres.fdr}, two lines parallel
         to the 45-degree line are generated in the SAM plot.}
    \item{ngenes}{a number or proportion of genes for which the FDR is estmated.}
    \item{iteration}{the number of iterations used in the estimation of the FDR
        for a given number or proportion of genes.}
    \item{initial.delta}{a set of initial guesses for \eqn{\Delta}{Delta} in the computation
         of the FDR for a given number or proportion of genes.}
    \item{rand}{if specified, the random number generator will be put in a 
        reproducible state.}
}
\value{a table of statistics (estimate of \eqn{p_0}{p0}, number of significant genes,
    number of falsely called genes and FDR) for the specified set of Deltas, a
    SAM Plot, a Delta vs.\ FDR plot, and a plot of Delta vs.\ the number of significant genes.}

\note{For further analyses with \code{sam.plot}, the results of \code{sam} must be assigned
    to an object.
    
    SAM was deveoped by Tusher et al. (2001).
    
    !!! There is a patent pending for the SAM technology at Stanford University. !!! 
}

\section{Warning}{In the one class case, the null distribution will only be computed correctly,
    if the expression values are log ratios. So in the one class case only log ratios should be
    used. (There will be no checking, if the expression values are
    really log ratios.)}

\seealso{
    \code{\link{sam.plot}} \code{\link{sam.wilc}} \code{\link{sam.delta}}
}

\references{
    Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
    applied to the ionizing radiation response, \emph{PNAS}, 98, 5116-5121.

    Storey, J.D. (2002). A direct approach to the false discovery rate, \emph{Journal of
    the Royal Statistical Society, Series B}, 64, 479-498.

    Storey, J.D., and Tibshirani, R. (2003). Statistical significance for genome-wide
    experiments, \emph{Technical Report}, Department of Statistics, Stanford University.

    Schwender, H. (2003). Assessing the false discovery rate in a statistical analysis of
    gene expression data, Chapter 5, \emph{Diploma thesis}, Department of Statistics,
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
}}


\keyword{htest}
