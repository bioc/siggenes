\name{ebamControl}
\alias{ebamControl}
\alias{find.a0Control}

\title{Further EBAM Arguments}
\description{
  Specifies most of the optional arguments of \code{ebam} and \code{find.a0}.
}
\usage{   
ebamControl(p0 = NA, p0.estimation = c("splines", "interval", "adhoc"), 
   lambda = NULL, ncs.value = "max", use.weights = FALSE)
   
find.a0Control(p0.estimation = c("splines", "adhoc", "interval"), 
   lambda = NULL, ncs.value = "max", use.weights = FALSE,
   n.chunk = 5, n.interval = 139, df.ratio = NULL)
}

\arguments{
  \item{p0}{a numeric value specifying the prior probability \eqn{p_0}{p0} that a gene is not
    differentially expressed. If \code{NA}, \code{p0} will be estimated automatically.}
  \item{p0.estimation}{either \code{"splines"} (default), \code{"interval"}, or \code{"adhoc"}. 
    If \code{"splines"}, the spline based method of Storey and Tibshirani (2003) is used to estimate
    \eqn{p_0}{p0}. If \code{"adhoc"} (\code{"interval"}), the adhoc (interval based) method 
    proposed by Efron et al.\ (2001) is used to estimate \eqn{p_0}{p0}.}
  \item{lambda}{a numeric vector or value specifying the \eqn{\lambda}{lambda} values used in
    the estimation of \eqn{p_0}{p0}. If \code{NULL}, \code{lambda} is set to \code{seq(0, 0.95, 0.05)}
    if \code{p0.estimation = "splines"}, and to \code{0.5} if \code{p0.estimation = "interval"}.
    Ignored if \code{p0.estimation = "adhoc"}. For details, see \code{\link{pi0.est}}.}
  \item{ncs.value}{a character string. Only used if \code{p0.estimation = "splines"} and
    \code{lambda} is a vector. Either \code{"max"} or \code{"paper"}. For details, see
    \code{\link{pi0.est}}.}
  \item{use.weights}{should weights be used in the spline based estimation of \eqn{p_0}{p0}? If
    \code{TRUE}, 1 - \code{lambda} is used as weights. For details, see \code{\link{pi0.est}}.}
  \item{n.chunk}{an integer specifying in how many subsets the \code{B} permutations
     should be split when computing the permuted test scores.}
  \item{n.interval}{the number of intervals used in the logistic regression with
     repeated observations for estimating the ratio \eqn{f_0/f}{f0/f}.}
  \item{df.ratio}{integer specifying the degrees of freedom of the natural cubic
     spline used in the logistic regression with repeated observations.}    
}
\details{
  These parameters should only be changed if they are fully understood.
}
\value{
  A list containing the values of the parameters that are used in \code{\link{ebam}} or \code{\link{find.a0}},
  respectively.
}
\references{   
   Efron, B., Tibshirani, R., Storey, J.D. and Tusher, V. (2001). Empirical Bayes Analysis
   of a Microarray Experiment. \emph{JASA}, 96, 1151-1160.
   
   Storey, J.D. and Tibshirani, R. (2003). Statistical Significance for Genome-Wide
   Studies. \emph{Proceedings of the National Academy of Sciences}, 100, 9440-9445.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{limma2ebam}}, \code{\link{ebam}}, \code{\link{find.a0}}}

\keyword{utilities}

