\name{impmv}
\alias{impmv}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Impute Missing Values}
\description{
  Imputes missing values in each row of an array based on the
  non-missing values in the row.
}
\usage{
impmv(data, IFUN = mean, min.nonmiss = 2, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{ A numeric matrix }
  \item{IFUN}{ Function that takes the nonmissing values in a row as input and
    returns the imputed values as output }
  \item{min.nonmiss}{ The minimum number of non-missing values needed to
    impute missing values }
  \item{\dots}{ Additional arguments to \code{IFUN} }
}
\details{
  \code{IFUN} can either return a single value (in which case each
  missing value is replaced by the same value), or a vector of length
  equal to the number of missing values in the row (see Examples below).
  The \code{min.nonmiss} argument can be used to apply different
  imputations for rows with different numbers of missing values (by
  making multiple calls to \code{impmv} with different values of
  \code{min.nonmiss}) or to remove rows with too few non-missing values
  (eg by applying \code{na.omit} after calling \code{impmv}).
}
\value{
  The input matrix \code{data} with dimension and dimnames preserved,
  but with missing values in rows with more than \code{min.nonmiss}
  non-missing values replaced by imputed values.
}
\examples{
x <- matrix(c(1:10,NA,NA),4,dimnames=list(format(1:4),c('a','b','c')))
impmv(x)
#   a b  c
# 1 1 5  9
# 2 2 6 10
# 3 3 7  5
# 4 4 8  6
na.omit(impmv(x,min.nonmiss=3))
#   a b  c
# 1 1 5  9
# 2 2 6 10
# attr(,"na.action")
# 3 4 
# 3 4 
# attr(,"na.action")attr(,"class")
# [1] "omit"
irn <- function(x,nc) rnorm(nc-length(x),mean(x),sqrt(var(x)))
impmv(x,IFUN=irn,nc=ncol(x)) 
#   a b         c
# 1 1 5  9.000000
# 2 2 6 10.000000
# 3 3 7  4.933022
# 4 4 8  5.039667
impmv(t(x),IFUN=irn,nc=nrow(x)) 
#   1  2       3        4
# a 1  2 3.00000  4.00000
# b 5  6 7.00000  8.00000
# c 9 10 9.76208 10.07452
}
\keyword{manip }% at least one, from doc/KEYWORDS
