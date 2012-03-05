\name{EBAM-class}
\docType{class}
\alias{EBAM-class}
\alias{EBAM}
\alias{show,EBAM-method}
\alias{print,EBAM-method}
\alias{plot,EBAM-method}
\alias{summary,EBAM-method}
\alias{plot,EBAM,ANY-method}

\title{Class EBAM}

\description{This is a class representation for the Empirical Bayes Analysis of
  Microarrays (EBAM) proposed by Efron et al. (2001).
}

\section{Objects from the Class}{
    Objects can be created using the function \code{ebam}.
}

\section{Slots}{
  \describe{
     \item{\code{z}:}{Object of class \code{"numeric"} representing the 
        expression scores of the genes.}
     \item{\code{posterior}:}{Object of class \code{"numeric"} representing
        the posterior probabilities of the genes.}
     \item{\code{p0}:}{Object of class \code{"numeric"} specifying the
        prior probability that a gene is not differentially expressed.}
     \item{\code{local}:}{Object of class \code{"numeric"} consisting of the
        local FDR estimates for the genes.}
     \item{\code{mat.fdr}:}{Object of class \code{"matrix"} containing general
        statistics such as the number of differentially expressed genes and 
        the estimated FDR for the specified values of \code{delta}.}
     \item{\code{a0}:}{Object of class \code{"numeric"} specifying the used
        value of the fudge factor. If not computed, \code{a0} will be
        set to \code{numeric(0)}.}
    \item{\code{mat.samp}:}{Object of class \code{"matrix"} containing
       the permuted group labels used in the estimation of the null
       distribution. Each row represents one permutation, each column
       one observation (pair). If no permutation procedure has been used,
       \code{mat.samp} will be set to \code{matrix(numeric(0))}.}
     \item{\code{vec.pos}:}{Object of class \code{"numeric"} consisting of
       the number of positive permuted test scores that are absolutely
       larger than the test score of a particular gene for each gene.
       If not computed \code{vec.pos} is set to \code{numeric(0)}.}
     \item{\code{vec.neg}:}{Object of class \code{"numeric"} consisting of
       the number of negative permuted test scores that are absolutely
       larger than the test score of a particular gene for each gene.
       If not computed \code{vec.neg} is set to \code{numeric(0)}.}
    \item{\code{msg}:}{Object of class \code{"character"} containing information
       about, e.g., the type of analysis. \code{msg} is printed when the functions
       \code{print} and \code{summary} are called.}
    \item{\code{chip}:}{Object of class \code{"character"} naming the microarray
       used in the analysis. If no information about the chip is available,
        \code{chip} will be set to \code{""}.}
}}

\section{Methods}{
  \describe{
     \item{plot}{\code{signature(object = "EBAM")}: Generates a plot of the
        posterior probabilities of the genes for a specified value of 
        \eqn{\Delta}{Delta}. For details, see \code{help.ebam(plot)}. For
        the arguments, see \code{args.ebam(plot)}.}
     \item{print}{\code{signature(object = "EBAM")}: Prints general information
        such as the number of differentially expressed genes and the estimated
        FDR for several values of \eqn{\Delta}{Delta}. For details, see
        \code{help.ebam(print)}. Arguments can be listed by \code{args.ebam(print)}.} 
     \item{show}{\code{signature(object = "EBAM")}: Shows the output of an
        EBAM analysis.}
     \item{summary}{\code{signature(object = "EBAM")}: Summarizes the results
        of an EBAM analysis for a specified value of \eqn{\Delta}{Delta}.
        For details, see \code{help.ebam(summary)}. For the arguments,
        see \code{args.ebam(summary)}.}
}}

\references{ 
   Efron, B., Tibshirani, R., Storey, J.D. and Tusher, V. (2001). Empirical Bayes Analysis
   of a Microarray Experiment, \emph{JASA}, 96, 1151-1160.
    
   Schwender, H., Krause, A. and Ickstadt, K. (2003). Comparison of
   the Empirical Bayes and the Significance Analysis of Microarrays.
   \emph{Technical Report}, SFB 475, University of Dortmund, Germany.
}


\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{\code{\link{ebam}}, \code{\link{find.a0}}, \code{\link{FindA0-class}}}

\examples{\dontrun{
  # Load the data of Golub et al. (1999) contained in the package multtest.
  data(golub)
  
  # golub.cl contains the class labels.
  golub.cl
  
  # Perform an EBAM analysis for the two class unpaired case assuming
  # unequal variances. Specify the fudge factor a0 by the suggested
  # choice of find.a0
  find.out <- find.a0(golub, golub.cl, rand = 123)
  ebam.out <- ebam(find.out)
  ebam.out
    
  # Obtain the number of differentially
  # expressed genes and the FDR if a gene is called differentially
  # expressed if its posterior probability is larger than 0.8, 0.85,
  # 0.9, 0.95.
  print(ebam.out, c(0.8, 0.85, 0.9, 0.95))
    
  # Generate a plot of the posterior probabilities for delta = 0.9.
  plot(ebam.out, 0.9)
    
  # Obtain the list of genes called differentially expressed if their
  # posterior probability is larger than 0.99, and gene-specific 
  # statistics for these variables such as their z-value and their
  # local FDR.
  summary(ebam.out, 0.9)
  
}}

\keyword{classes}
