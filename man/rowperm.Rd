\name{rowperm}
\alias{rowperm}
\title{
Generates random permutations within rows of a matrix
}
\description{
Given a matrix x, returns a copy of x with a separate random
permutations applied within each row.
}
\usage{
rowperm(x)
}
\arguments{
\item{x}{
a numeric matrix
A numeric matrix with each row containing the elements from the
corresponding row of x permuted in a random order.
}}
\examples{
 x <- matrix(1:12,3)
 x
#     [,1] [,2] [,3] [,4] 
#[1,]    1    4    7   10
#[2,]    2    5    8   11
#[3,]    3    6    9   12
  rowperm(x)
#     [,1] [,2] [,3] [,4] 
#[1,]   10    7    1    4
#[2,]    8   11    2    5
#[3,]    6   12    9    3

}
\keyword{utilities}
% Converted by Sd2Rd version 1.19.
