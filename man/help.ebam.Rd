\name{help.ebam}
\alias{help.ebam}
\alias{args.ebam}
\title{Help files or argument list for EBAM-specific methods }
\description{Displays the help page or the argument list, respectively,
    for a EBAM-specific method.
}
\usage{
   help.ebam(method)
   args.ebam(method)
}
\arguments{
  \item{method}{a name or a character string specifying the method for which 
    the arguments or the help page, respectively, should be shown. Currently 
    available are \code{print}, \code{plot}, and \code{summary}}
}
\value{
  The arguments of the specified \code{method} are displayed or a html page
  containing the help for the specified \code{method} is opened, respectively.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{EBAM-class}}, \code{\link{ebam}}
}
\examples{\dontrun{
   # Displays the arguments of the function summary
   args.ebam(summary)
   
   # Opens the help page in the browser
   help.ebam(summary)
}}
\keyword{documentation}
