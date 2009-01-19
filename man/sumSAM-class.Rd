\name{sumSAM-class}
\docType{class}
\alias{sumSAM-class}
\alias{print,sumSAM-method}
\alias{show,sumSAM-method}
\alias{sumEBAM-class}
\alias{print,sumEBAM-method}
\alias{show,sumEBAM-method}

\title{Classes sumSAM and sumEBAM}
\description{These classes are just used for a nicer output of the summary of
   an object of class SAM or EBAM, respectively.}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("sumSAM", ...)},
    or by using the function \code{summary(object)} when object is a 
    SAM-class object.  
    
    Objects can be created by calls of the form \code{new("sumEBAM", ...)},
    or by using the function \code{summary(object)} when object is an 
    EBAM-class object.      
}
\section{Slots}{
  \describe{
    \item{\code{row.sig.genes}:}{Object of class \code{"numeric"} consisting
       of the row numbers of the significant genes in the data matrix.}
    \item{\code{mat.fdr}:}{Object of class \code{"matrix"} containing general
       information as the number of differentially expressed genes and the
       estimated FDR for either one or several values of Delta.}
    \item{\code{mat.sig}:}{Object of class \code{"data.frame"} containing gene-specific
       statistics as the d-values (or z-values) and the q-values or (the local FDR)
       of the differentially expressed genes.}
    \item{\code{list.args}:}{Object of class \code{"list"} consisting of some of
       the specified arguments of summary needed for internal use.}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "sumSAM")}: Prints the output of the SAM-specific
       method summary.}
    \item{show}{\code{signature(object = "sumSAM")}: Shows the output of the summary
       of a SAM analysis.}
       
    \item{print}{\code{signature(x = "sumEBAM")}: Prints the output of the EBAM-specific
       method summary.}
    \item{show}{\code{signature(object = "sumEBAM")}: Shows the output of the summary
       of a EBAM analysis.}
  }
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{SAM-class}}, \code{\link{EBAM-class}}
}


\keyword{classes}
