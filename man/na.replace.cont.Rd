\name{na.replace.cont}
\alias{na.replace.cont}
\title{Replacing of Missing Values in Continuous Data}
\description{
  Replaces missing values by either the rowwise mean or median.
}
\usage{
  na.replace.cont(X, na.method = "mean")
}
\arguments{
  \item{X}{a matrix}
  \item{na.method}{character string specifying the statistic
     with which the missing values should be replaced. Must be either
     \code{"mean"} or \code{"median"}}
}
\value{
   The matrix \code{X} with missing values replaced by either the
   rowwise mean or median.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{na.handling}},\code{\link{na.replace.dist}}
}
\keyword{utilities}
