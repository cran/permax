\name{summary.permax}
\alias{summary.permax}
\alias{print.summary.permax}
\title{
Summarizes the output of permax
}
\description{
Finds and prints the most significant genes in the output of
\code{\link{permax}} 
}
\usage{
summary.permax(object, data, nl=25, nr=25, ...)

print.summary.permax(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{object}{
A data frame of class \code{permax} (created by \code{\link{permax}})
}
\item{data}{
The data matrix used as input to \code{permax}.  If given, the rows of
\code{data} corresponding to the most significant genes will also be printed.
}
\item{nl}{
The \code{nl} most significant genes in the lower tail will be printed
}
\item{nr}{
The \code{nr} most significant genes in the upper tail will be printed
}
\item{x}{
An object of class \code{summary.permax}
}
\item{digits}{
Number of significant digits in the printed output (roughly)
}
\item{...}{
  Supplied for compatibility, but not used.
}
}
\value{
  \code{summary.permax} returns a list with components
  \item{nl,nr}{the corresponding input values}
  \item{objectl}{the \code{nl} most significant lower tail rows of \code{object}}
  \item{objectr}{the \code{nr} most significant upper tail rows of \code{object}}
  \item{dataL}{the rows of \code{data} corresponding to \code{$objectl}
    (if \code{data} specified)}
  \item{dataR}{the rows of \code{data} corresponding to \code{$objectr}
    (if \code{data} specified)}
  \item{lower}{giving the \code{row.names} of \code{objectl}} 
  \item{upper}{giving the \code{row.names} of \code{objectr}} 
  \item{dist}{giving the \code{dist} attribute of \code{object}}
  \item{call}{giving the \code{call} attribute of \code{object}}
  
\code{print.summary.permax} prints the \code{$call}, \code{$objectl},
and \code{$objectr}, components, and if non-NULL, the \code{$dist},
\code{$dataL}, and \code{$dataR} components.
}
\details{
If \code{data} is given, it must be a data frame with
\code{row.names(data)} corresponding to \code{row.names(object)}, or a
matrix with \code{dimnames(data)[[1]]} corresponding to
\code{row.names(object)}.  The purpose of including \code{data} is
primarily to print the rows of the original data corresponding to the
most significant statistics.  
}
\seealso{
\code{\link{permax}}
}
\examples{
# An example is given in the permax() help file
}
\keyword{print}
% Converted by Sd2Rd version 1.19.
