\name{sumSAM-class}
\docType{class}
\alias{sumSAM-class}
\alias{print,sumSAM-method}

\title{Class sumSAM}
\description{This class is just used for a nicer output of the summary of
   an SAM-class object.}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("sumSAM", ...)}
    or by using the function \code{summary(object)} when object is an 
    SAM-class object. 
}
\section{Slots}{
  \describe{
    \item{\code{row.sig.genes}:}{Object of class \code{"numeric"} consisting
       of the row numbers of the significant genes in the data matrix}
    \item{\code{mat.fdr}:}{Object of class \code{"matrix"} containing general
       information as the number of differentially expressed genes and the
       estimated FDR for either one or several values of Delta.}
    \item{\code{mat.sig}:}{Object of class \code{"data.frame"} containing gene-specific
       statistics as the d-values and the q-values of the differentially expressed
       genes.}
    \item{\code{list.args}:}{Object of class \code{"list"} consisting of some of
       the specified arguments of summary needed for internal use.}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "sumSAM")}: Prints the output of the SAM-specific
       method summary.}
  }
}
\note{ 
   SAM was developed by Tusher et al. (2001).
    
   !!! There is a patent pending for the SAM technology at Stanford University. !!!
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{
  \code{\link{SAM-class}}, \code{\link{summary,SAM-method}}
}


\keyword{classes}
