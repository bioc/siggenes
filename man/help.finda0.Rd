\name{help.finda0}
\alias{help.finda0}
\alias{args.finda0}
\title{Help files or argument list for FindA0-specific methods }
\description{Displays the help page or the argument list, respectively,
    for a FindA0-specific method.
}
\usage{
   help.finda0(method)
   args.finda0(method)
}
\arguments{
  \item{method}{a name or a character string specifying the method for which 
    the arguments or the help page, respectively, should be shown. Currently 
    available are \code{print} and \code{plot}.}
}
\value{
  The arguments of the specified \code{method} are displayed or a html page
  containing the help for the specified \code{method} is opened, respectively.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{FindA0-class}}, \code{\link{find.a0}}
}
\examples{\dontrun{
   # Displays the arguments of the function summary
   args.finda0(summary)
   
   # Opens the help page in the browser
   help.finda0(summary)
}}
\keyword{documentation}
