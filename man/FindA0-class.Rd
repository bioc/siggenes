\name{FindA0class}
\docType{class}
\alias{FindA0-class}
\alias{FindA0}
\alias{show,FindA0-method}
\alias{print,FindA0-method}
\alias{plot,FindA0-method}

\title{Class FindA0}

\description{This is a class representation for the specification of the fudge factor
  in an EBAM analysis as proposed by Efron et al. (2001).
}

\section{Objects from the Class}{
    Objects can be created using the function \code{find.a0}.
}

\section{Slots}{
  \describe{
     \item{\code{mat.z}:}{Object of class \code{"matrix"} containing the expression
        scores of the genes for each of the possible values for the fudge factor,
        where each row corresponds to a gene, and each column to one of the values
        for the fudge factor \eqn{a_0}{a0}.}
     \item{\code{mat.posterior}:}{Object of class \code{"matrix"} consisting of
        the posterior probabilities of the genes for each of the possible values 
        for the fudge factor, where each row of \code{mat.posterior} corresponds to
        a gene, and each column to one of the values for \eqn{a_0}{a0}. The 
        probabilities in \code{mat.posterior} are computed using the monotonically
        transformed test scores (see the Details section of \code{\link{find.a0}}).}
     \item{\code{mat.center}:}{Object of class \code{"matrix"} representing the
        centers of the \code{nrow(mat.center)} intervals used in the logistic
        regression with repeated observations for estimating \eqn{f/f_0}{f/f0}
        for each of the \code{ncol(mat.center)} possible values for the fudge 
        factor.}
     \item{\code{mat.success}:}{Object of class \code{"matrix"} consisting of 
        the numbers of observed test scores in the \code{nrow(mat.success)} intervals
        used in the logistic regression with repeated observations for each
        of the \code{ncol(mat.success)} possible values for the fudge factor.}
     \item{\code{mat.failure}:}{Object of class \code{"matrix"} containing the
        numbers of permuted test scores in the \code{nrow(mat.failure)} intervals
        used in the logistic regression with repeated observations for each
        of the \code{ncol(mat.failure)} possible values for the fudge factor.}
     \item{\code{z.norm}:}{Object of class \code{"numeric"} comprising the
        values of the \code{nrow(mat.z)} quantiles of the standard normal 
        distribution (if any \code{mat.z<0}) or an F-distribution (if all
        \code{mat.z >= 0}).}
     \item{\code{p0}:}{Object of class \code{"numeric"} specifying the prior
        probability that a gene is not differentially expressed.}
     \item{\code{mat.a0}:}{Object of class \code{"data.frame"} comprising
        the number of differentially expressed genes and the estimated FDR
        for the possible choices of the fudge factor specified by \code{vec.a0}.}
     \item{\code{mat.samp}:}{Object of class \code{"matrix"} consisting of the
        \code{nrow{mat.samp}} permutations of the class labels.}
     \item{\code{vec.a0}:}{Object of class \code{"numeric"} representing the
        possible values of the fudge factor \eqn{a_0}{a0}.}
     \item{\code{suggested}:}{Object of class \code{"numeric"} revealing the
        suggested choice for the fudge factor, i.e. the value of \code{vec.a0}
        that leads to the largest number of differentially expressed genes.}
     \item{\code{delta}:}{Object of class \code{"numeric"} specifying the
        minimum posterior probability that a gene must have to be called 
        differentially expressed.}
     \item{\code{df.ratio}:}{Object of class \code{"numeric"} representing the
        degrees of freedom of the natural cubic spline used in the logistic
        regression with repeated observations.}
    \item{\code{msg}:}{Object of class \code{"character"} containing information
       about, e.g., the type of analysis. \code{msg} is printed when 
       \code{print} is called.}
    \item{\code{chip}:}{Object of class \code{"character"} naming the microarray
       used in the analysis. If no information about the chip is available,
        \code{chip} will be set to \code{""}.}
}}

\section{Methods}{
  \describe{
     \item{plot}{\code{signature(object = "FindA0")}: Generates a plot of the
        (logit-transformed) posterior probabilities of the genes for a specified 
        value of \eqn{\Delta}{Delta} and a set of possible values for the fudge
        factor. For details, see \code{help.finda0(plot)}. For
        the arguments, see \code{args.finda0(plot)}.}
     \item{print}{\code{signature(object = "FindA0)}: Prints the number of 
        differentially expressed genes and the estimated
        FDR for each of the possible values of the fudge factor specified by 
        \code{vec.a0}. For details, see \code{help.finda0(print)}. 
        For arguments, see \code{args.finda0(print)}.} 
     \item{show}{\code{signature(object = "FindA0")}: Shows the output of an
        analysis with \code{find.a0}.}
}}

\references{ 
   Efron, B., Tibshirani, R., Storey, J.D. and Tusher, V. (2001). Empirical Bayes Analysis
   of a Microarray Experiment, \emph{JASA}, 96, 1151-1160.
    
   Schwender, H., Krause, A. and Ickstadt, K. (2003). Comparison of
   the Empirical Bayes and the Significance Analysis of Microarrays.
   \emph{Technical Report}, SFB 475, University of Dortmund, Germany.
}


\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{\code{\link{find.a0}}, \code{\link{ebam}}, \code{\link{EBAM-class}}}

\keyword{classes}


      
