\name{nclass.wand}
\alias{nclass.wand}

\title{Number of cells in a histogram}
\description{
  Computes the number of cells in a histogram using the method of Wand (1994).
}
\usage{
   nclass.wand(x, level = 1)
}

\arguments{
  \item{x}{numeric vector of observations}
  \item{level}{integer specifying the number of levels of functional estimation used in
     the estimation. For details, see the help page of \code{dpih} from the package
     \pkg{KernSmooth}}
}

\details{
  \code{nclass.wand} calls \code{dpih}, and then computes the number of cells
  corresponding to the optimal bin width returned by \code{dpih}
}

\value{
  a numeric value specifying the number of cells for the histogram of \code{x}
}

\references{  
  Wand, M.P. (1997). Data-based choice of histogram bin width. 
  \emph{American Statistician}, 51, 59--64. 
}

\seealso{\code{\link{denspr}}}

\keyword{utilities}
