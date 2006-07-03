\name{na.replace.dist}
\alias{na.replace.dist}
\title{Handling of Missing Values in Categorical Data}
\description{
  Replaces missing values by the distribution of the levels.
}
\usage{
  na.replace.dist(x)
}
\arguments{
  \item{x}{a vector}
}
\details{
  Each missing value is replaced by a random draw from the levels
  of \code{x}, where the probability that a specific level is
  drawn is given by the number of non-missing values of \code{x}
  having this level divided by the number of
  non-missing values of \code{x}. 
}
\value{
  the vector \code{x} with missing values replaced by the described
  method
}
\author{Holger Schwender, \email{holger.schw@gmx}}
\seealso{
   \code{\link{na.handling}},\code{\link{na.replace.dist}}
}
\keyword{utilities}
