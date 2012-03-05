\name{SAM-class}
\docType{class}
\alias{SAM-class}
\alias{SAM}
\alias{identify,SAM-method}
\alias{plot,SAM-method}
\alias{print,SAM-method}
\alias{show,SAM-method}
\alias{summary,SAM-method}
\alias{plot,SAM,ANY-method}

\title{Class SAM}
\description{This is a class representation for several versions of the
   SAM (Significance Analysis of Microarrays) procedure proposed by Tusher
   et al. (2001).}
\section{Objects from the Class}{
   Objects can be created using the functions \code{sam}, \code{sam.dstat}, 
   \code{sam.wilc} and \code{sam.snp}.
}
\section{Slots}{
  \describe{
    \item{\code{d}:}{Object of class \code{"numeric"} representing the
       expression scores of the genes.}
    \item{\code{d.bar}:}{Object of class \code{"numeric"} representing
       the expected expression scores under the null hypothesis.}
    \item{\code{vec.false}:}{Object of class \code{"numeric"} containing
       the one-sided expected number of falsely called genes.}
    \item{\code{p.value}:}{Object of class \code{"numeric"} consisting of
       the p-values of the genes.}
    \item{\code{s}:}{Object of class \code{"numeric"} representing the
       standard deviations of the genes. If the standard deviations are
       not computed, \code{s} will be set to \code{numeric(0)}. }
    \item{\code{s0}:}{Object of class \code{"numeric"} representing the
       value of the fudge factor. If not computed, \code{s0} will be
       set to \code{numeric(0)}.}
    \item{\code{mat.samp}:}{Object of class \code{"matrix"} containing
       the permuted group labels used in the estimation of the null
       distribution. Each row represents one permutation, each column
       one observation (pair). If no permutation procedure has been used,
       \code{mat.samp} will be set to \code{matrix(numeric(0))}.}
    \item{\code{p0}:}{Object of class \code{"numeric"} representing the
       prior probability that a gene is not differentially expressed.}
    \item{\code{mat.fdr}:}{Object of class \code{"matrix"} containing general
       information as the number of significant genes and the estimated FDR
       for several values of \eqn{\Delta}{Delta}. Each row represents one
       value of \eqn{\Delta}{Delta}, each of the 9 columns one statistic.}
    \item{\code{q.value}:}{Object of class \code{"numeric"} consisting of
       the q-values of the genes. If not computed, \code{q.value} will be
       set to \code{numeric(0)}.}
    \item{\code{fold}:}{Object of class \code{"numeric"} representing the
       fold changes of the genes. If not computed, \code{fold} will be 
       set to \code{numeric(0)}.}
    \item{\code{msg}:}{Object of class \code{"character"} containing information
       about, e.g., the type of analysis. \code{msg} is printed when the functions
       \code{print} and \code{summary}, respectively, are called.}
    \item{\code{chip}:}{Object of class \code{"character"} naming the microarray
       used in the analysis. If no information about the chip is available,
        \code{chip} will be set to \code{""}.}
  }
}
\section{Methods}{
  \describe{
    \item{identify}{\code{signature(x = "SAM")}: After generating a SAM plot,
       \code{identify} can be used to obtain information about the genes by
       clicking on the symbols in the SAM plot. For details, see 
       \code{help.sam(identify)}. Arguments are listed by \code{args.sam(identify)}.}
    \item{plot}{\code{signature(x = "SAM")}: Generates a SAM plot or the Delta
       plots. If the specified \code{delta} in \code{plot(object,delta)} is
       a numeric value, a SAM plot will be generated. If \code{delta} is either
       not specified or a numeric vector, the Delta plots will be generated.
       For details, see \code{?sam.plot2}, \code{?delta.plot} or 
       \code{help.sam(plot)},respectively. Arguments are listed by \code{args.sam(plot)}.}
    \item{print}{\code{signature(x = "SAM")}: Prints general information such as 
       the number of significant genes and the estimated FDR for a set of 
       \eqn{\Delta}{Delta}. For details, see \code{help.sam(print)}. Arguments are
       listed by \code{args.sam(print)}.}
    \item{show}{\code{signature(object = "SAM")}: Shows the output of the SAM
       analysis.}
    \item{summary}{\code{signature(object = "SAM")}: Summarizes the results of
        a SAM analysis. If \code{delta} in \code{summary(object,delta)} is not
        specified or a numeric vector, the information shown by print and some
        additional information will be shown. If \code{delta} is a numeric
        vector, the general information for the specific \eqn{\Delta}{Delta} is
        shown and additionally gene-specific information about the genes called 
        significant using this value of \eqn{\Delta}{Delta}. The output of summary
    is an object of class sumSAM which has the slots \code{row.sig.genes},
    \code{mat.fdr}, \code{mat.sig} and \code{list.args}. For details, 
    see \code{help.sam(summary)}. All arguments are listed by \code{args.sam(summary)}.}
  }
}
\note{ 
   SAM was developed by Tusher et al. (2001).
    
   !!! There is a patent pending for the SAM technology at Stanford University. !!!
}

\references{
   Schwender, H., Krause, A. and Ickstadt, K. (2003). Comparison of
   the Empirical Bayes and the Significance Analysis of Microarrays.
   \emph{Technical Report}, SFB 475, University of Dortmund, Germany.
   \url{http://www.sfb475.uni-dortmund.de/berichte/tr44-03.pdf}.

   Schwender, H. (2004). Modifying Microarray Analysis Methods for 
   Categorical Data -- SAM and PAM for SNPs. To appear in: \emph{Proceedings
   of the the 28th Annual Conference of the GfKl}.

   Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance analysis of microarrays
   applied to the ionizing radiation response. \emph{PNAS}, 98, 5116-5121.

}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{sam}},\code{\link{args.sam}},\code{\link{sam.plot2}},
  \code{\link{delta.plot}}
}
\examples{\dontrun{
  # Load the package multtest and the data of Golub et al. (1999)
  # contained in multtest.
  library(multtest)
  data(golub)
  
  # Perform a SAM analysis for the two class unpaired case assuming
  # unequal variances.
  sam.out <- sam(golub, golub.cl, B=100, rand=123)
  sam.out
  
  # Alternative ways to show the output of sam.
  show(sam.out)
  print(sam.out)
  
  # Obtain a little bit more information.
  summary(sam.out)
  
  # Print the results of the SAM analysis for other values of Delta.
  print(sam.out, seq(.2, 2, .2))
  
  # Again, the same with additional information.
  summary(sam.out, seq(.2, 2, .2))
    
  # Obtain the Delta plots for the default set of Deltas.
  plot(sam.out)
  
  # Generate the Delta plots for Delta = 0.2, 0.4, 0.6, ..., 2.
  plot(sam.out, seq(0.2, 0.4, 2))
  
  # Obtain the SAM plot for Delta = 2.
  plot(sam.out, 2)
  
  # Get information about the genes called significant using 
  # Delta = 3.
  sam.sum3 <- summary(sam.out, 3)
  sam.sum3
  
  # Obtain the rows of the Golub et al. (1999) data set containing
  # the genes called differentially expressed
  sam.sum3@row.sig.genes
  
  # and their names
  golub.gnames[sam.sum3@row.sig.genes, 3] 

  # The matrix containing the d-values, q-values etc. of the
  # differentially expressed genes can be obtained by
  sam.sum3@mat.sig
}}
\keyword{classes}
