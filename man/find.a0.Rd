\name{find.a0}
\alias{find.a0}
\title{Computation of the Fudge Factor}
\description{
    Provides the required information for obtaining the optimal choice of the fudge
    factor in the Empirical Bayes Analysis of Microarrays that uses the modified
    t statistics.
}
\usage{
    find.a0(data,x,y,paired=FALSE,mat.samp=NULL,B=100,balanced=FALSE,na.rm=FALSE,delta=0.9,
    alpha=(0:9)/10,include.0=TRUE,p0=NA,stable=TRUE,number.int=139,rand=NA,plot.legend=TRUE)
}

\arguments{
    \item{data}{the data set that should be analyzed. Every row of this data set must
        correspond to a gene.}
    \item{x}{vector of the columns of \code{data} that correspond to the treatment group.
        In the paired case, \eqn{(x[i],y[i])} build a pair. If, e.g., the first \eqn{n_1}{n1}
        columns contain the gene expression values of the treatment group,\code{x=1:\eqn{n_1}{n1}}.}
    \item{y}{vector of the columns of \code{data} that correspond to the control group. In the
        paired case, \eqn{(x[i], y[i])} are an observation pair.}
    \item{paired}{paired (\code{TRUE}) or unpaired (\code{FALSE}) data. Default is \code{FALSE}.}
    \item{mat.samp}{a permutation matrix. If specified, this matrix will be used, even if
        \code{rand} and \code{B} are specified.}
    \item{B}{number of permutations used in the calculation of the null density.}
    \item{balanced}{if \code{TRUE}, only balanced permutations will be used. Default is
        \code{FALSE}.}
    \item{na.rm}{if \code{FALSE} (default), the expression score of genes with one or more
        missing values will be set to \code{NA}. If \code{TRUE}, the missing values
        will be replaced by the genewise mean of the non-missing values.}
    \item{delta}{a gene will be called differentially expressed, if its posterior
        probability of being differentially expressed is large than or equal to
        \code{delta}.}
    \item{alpha}{a vector of possible values for the fudge factor \eqn{a_0}{a0} in terms of quantiles of the
        standard deviations of the genes.}
    \item{include.0}{if \code{TRUE} (default), \eqn{a_0=0}{a0=0} will also be a possible choice
        of the fudge factor.}
    \item{p0}{the prior probability that a gene is differentially expressed. If not specified,
        it will automatically be computed.}
    \item{stable}{if \code{TRUE} (default), \eqn{p_0}{p0} will be computed by the algorithm of
        Storey and Tibshirani (2003). If \code{FALSE}, the (unstable) estimate that
        ensures that the posterior probability of being differentially expressed is
        always non-negative is computed.}
    \item{number.int}{number of equally spaced intervals between the minimum and and the maximum
        of the expression scores \eqn{z} that are used in the logistic regression for estimating
        the ratio of the null density to the mixture density.}
    \item{rand}{if specified, the random number generator will be put in a reproducible
        state.}
    \item{plot.legend}{if \code{TRUE} (default), a legend will be added to the plot of the
        expression scores vs.\ their logit-transformed posterior probability.}
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
    Efron, B., Tibshirani, R., Storey, J.D., and Tusher, V. (2001). Empirical Bayes Analysis
    of a Microarray Experiment, \emph{JASA}, 96, 1151-1160.

    Storey, J.D., and Tibshirani, R. (2003). Statistical significance for genome-wide
    experiments, \emph{Technical Report}, Department of Statistics, Stanford University.

    Schwender, H. (2003). Assessing the false discovery rate in a statistical analysis of
    gene expression data, Chapter 7, \emph{Diploma thesis}, Department of Statistics,
    University of Dortmund, \url{http://de.geocities.com/holgerschw/thesis.pdf}.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\keyword{htest}
