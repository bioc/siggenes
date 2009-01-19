\name{help.sam}
\alias{help.sam}
\alias{args.sam}
\title{Help files or argument list for SAM-specific methods }
\description{Displays the help page or the argument list, respectively,
    for a SAM-specific method.
}
\usage{
   help.sam(method)
   args.sam(method)
}
\arguments{
  \item{method}{a name or a character string specifying the method for which 
    the arguments or the help page, respectively, should be shown. Currently 
    available are \code{print}, \code{plot}, \code{summary} and \code{identify}.}
}
\value{
  The arguments of the specified \code{method} are displayed or a html page
  containing the help for the specified \code{method} is opened, respectively.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{SAM-class}},\code{\link{sam}}
}
\examples{\dontrun{
   # Displays the arguments of the function summary
   args.sam(summary)
   
   # Opens the help page in the browser
   help.sam(summary)
}}
\keyword{documentation}
