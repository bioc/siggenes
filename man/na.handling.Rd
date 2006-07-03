\name{na.handling}
\alias{na.handling}
\title{Handling of Missing Values}
\description{
  Removes missing values from a data set either by replacing them or
  by removing the rows that contain missing values.
}
\usage{
  na.handling(X, na.replace = TRUE, na.method = "mean")
}
\arguments{
  \item{X}{a matrix}
  \item{na.replace}{If \code{TRUE}, missing values are replaced by
     the rowwise statistic specified by \code{na.method}. If 
     \code{FALSE}, all rows with one or more missing values are
     removed}
  \item{na.method}{character string specifying the statistic
     with which the missing values should be replaced. Must be either
     \code{"mean"} or \code{"median"}}
}
\details{
   Since missing values are replaced by either the rowwise mean or median,
   \code{na.handling} can only handle continuous data. There is, however, no
   check if the data is continuous or categorical/discrete.
   
   If there are no or only one non-missing value in a row, \code{na.handling}
   will remove this row from the data set even if \code{na.replace=TRUE}.
}
\value{
  \item{X}{matrix without missing values. The number of rows of \code{X} can
     differ from the number of rows of the input matrix \code{X}, since all
     rows with either less than one non-missing value (if \code{na.replace=TRUE})
     or with at least one missing value (if \code{na.replace=FALSE}) are
     removed}
  \item{NA.genes}{the rows of the input matrix \code{X} that are removed and
     thus are excluded from the output matrix \code{X}}
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}
\seealso{
   \code{\link{na.replace.cont}},\code{\link{na.replace.dist}}
}
\examples{
   mat<-matrix(1:20,5,4)
   mat[1,2]<-NA
   mat
   
   # Replace the missing value by the mean of the first row.
   na.handling(mat)
   
   # Remove the first row.
   na.handling(mat,na.replace=FALSE)
   
   # Replace the missing value in the first row by the mean of
   # this row, and remove the second row containing only NAs.
   mat[2,]<-NA
   mat
   na.handling(mat)
}
\keyword{utilities}
