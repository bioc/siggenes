\name{find.a0}
\alias{find.a0}
\title{Computation of the Fudge Factor}
\description{
        Provides the required information for obtaining the optimal choice of the fudge factor
    in the Empirical Bayes Analysis of Microarrays that uses the modified t statistics.}
\usage{
        find.a0(data,cl,B=100,balanced=FALSE,mat.samp=NULL,delta=0.9,alpha=(0:9)/10,
        include.0=TRUE,p0=NA,plot.legend=TRUE,na.rm=FALSE,rand=TRUE)
}

\arguments{
    \item{data}{the data set that should be analyzed. Every row of this data set must
            correspond to a gene, and each column to a biological sample.}
    \item{cl}{a vector containing the class labels of the samples. In the two class unpaired case,
        the label of a sample is either 0 (e.g., control group) or 1 (e.g., case group). 
        In the two class paired case, the labels are the integers between 1 and \eqn{n/2}
        (e.g., before treatment group) and between -1 and \eqn{-n/2} (e.g., after treatment
        group), where \eqn{n} is the length of \code{cl} and \eqn{k} is paired with \eqn{-k}.
        For one group data, the label for each sample should be 1.}
    \item{B}{number of permutations used in the calculation of the null density.}
    \item{balanced}{if \code{TRUE}, only balanced permutations will be used. Default is
            \code{FALSE}.}
    \item{mat.samp}{a permutation matrix. If specified, this matrix will be used, even if
            \code{rand} and \code{B} are specified.}
    \item{delta}{a gene will be called differentially expressed, if its posterior
        probability of being differentially expressed is large than or equal to
            \code{delta}.}    
    \item{alpha}{a vector of possible values for the fudge factor \eqn{a_0}{a0} in terms
        of quantiles of the standard deviations of the genes.}
    \item{include.0}{if \code{TRUE} (default), \eqn{a_0=0}{a0=0} will also be a possible
        choice for the fudge factor.}
        \item{p0}{the prior probability that a gene is differentially expressed. If not
         specified, it will automatically be computed.}
    \item{plot.legend}{if \code{TRUE} (default), a legend will be added to the plot of the
            expression scores vs.\ their logit-transformed posterior probability.}
    \item{na.rm}{if \code{FALSE} (default), the expression score of genes with one or more
            missing values will be set to \code{NA}. If \code{TRUE}, the missing values
            will be replaced by the genewise mean of the non-missing values.}
    \item{rand}{if specified, the random number generator will be put in a reproducible
            state.}
}

\note{The results of \code{find.a0} must be assigned to an object for the further analysis
    with \code{ebam}.}
\value{a list of the numbers of genes called differentially expressed by the EBAM analysis
        for several choices of \eqn{a_0}{a0}, and the plot of the expression scores vs.\ their
        corresponding logit-transformed posterior probability of being significant.
    \item{sig.a0}{vector containing the number of differentially expressed genes for the
            specified set of values for \eqn{a_0}{a0}.}
    \item{a0}{the optimal choice of the fudge factor using the criterion of Efron et al. (2001)
            that the \eqn{a_0}{a0} should be used which leads to the most differentially
            expressed genes.}
}
\seealso{
    \code{\link{ebam}}  \code{\link{ebam.wilc}}
}

\references{
        Efron, B., Tibshirani, R., Storey, J.D., and Tusher, V. (2001). Empirical Bayes
    Analysis of a Microarray Experiment, \emph{JASA}, 96, 1151-1160.

    Storey, J.D., and Tibshirani, R. (2003). Statistical significance for genome-wide
        experiments, \emph{Technical Report}, Department of Statistics, Stanford University.

    Schwender, H. (2003). Assessing the false discovery rate in a statistical analysis of
    gene expression data, Chapter 7, \emph{Diploma thesis}, Department of Statistics,
    University of Dortmund, \url{http://de.geocities.com/holgerschw/thesis.pdf}.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\examples{\dontrun{
    library(multtest)
    # Load the data of Golub et al. (1999). data(golub) contains a 3051x38 gene expression
    # matrix called golub, a vector of length called golub.cl that consists of the 38 class labels,
    # and a matrix called golub.gnames whose third column contains the gene names.
    data(golub)
    
    # Now the optimal value for the fudge factor a0 is computed, where possible values of the a0 are
    # 0 and the 0, 0.05 and 0.1 quantile of the standard deviations of the genes. Setting rand=123
    # makes the results reproducible.
    
    find.out<-find.a0(golub,golub.cl,alpha=c(0,0.05,0.1),rand=123)
}}



\keyword{htest}
