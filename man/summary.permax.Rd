\name{summary.permax}
\alias{summary.permax}
\alias{print.summary.permax}
\title{
Summarizes the output of permax
}
\description{
Finds and prints the most significant genes in the output of \code{permax}
}
\usage{
summary.permax(object, data, nl=25, nr=25, dig=3, ...)
}
\arguments{
\item{object}{
A dataframe of class permax (created by \code{permax})
}
\item{data}{
The data matrix used as input to permax.  If given, the rows of d
corresponding to the most significant genes will also be printed.
}
\item{nl}{
The nl most significant genes in the lower tail will be printed
}
\item{nr}{
The nr most significant genes in the upper tail will be printed
}
\item{dig}{
Number of significant digits in the printed output (roughly)
Returns a list with components \code{lower} giving the \code{row.names}
of object 
corresponding to the nl most significant genes in the upper tail, and
\code{upper} giving the \code{nr} most significant genes in the upper
tail.}
\item{...} Supplied for compatibility but not used.
}  }
\details{
If d is given, it must be a data frame with \code{row.names(d)} corresponding
to \code{row.names(object)}, or a matrix with \code{dimnames(d)[[1]]}
corresponding to \code{row.names(object)}.
The purpose of including d is primarily to print the
rows of the original data corresponding to the most significant
statistics.  
}
\seealso{
permax
}
\examples{
# An example is given in the permax help file
}
\keyword{print}
% Converted by Sd2Rd version 1.19.
