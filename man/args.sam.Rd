\name{args.sam}
\alias{args.sam}
\title{Argument list for SAM-specific methods }
\description{Displays the argument names of the methods
   for a object of class SAM.
}
\usage{
   args.sam(method)
}
\arguments{
  \item{method}{a character string specifying the method whose arguments
      should be shown. \code{name} must be either \code{"print"}, \code{"plot"},
      \code{"summary"} or \code{"identify"}}
}
\value{
  The arguments of the specified \code{method} are displayed
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{SAM-class}},\code{\link{sam}},\code{\link{sam.dstat}},
   \code{\link{sam.wilc}},\code{\link{sam.snp}}
}
\examples{
   # Displays the arguments of the function summary
   args.sam("summary")
}
\keyword{documentation}
