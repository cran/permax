\name{plot.permax}
\alias{plot.permax}
\title{
Image plot of the most significant genes (attributes) from a permax
analysis 
}
\description{
  Given the output of \code{\link{permax}}, and the array of expression
  levels, creates a color image plot of the expression levels of the most
  significant genes 
}
\usage{
plot.permax(x, data, nl=25, nr=25, logs=TRUE, ig1=NULL, ig2=NULL,
            clmn.lab=dimnames(data)[[2]], row.lab=dimnames(data)[[1]],
            clmn.off=NULL, row.off=NULL,\dots)
}
\arguments{
\item{x}{
A permax object (output from \code{\link{permax}} or \code{\link{permcor}})
}
\item{data}{
Matrix or data frame of expression levels used as input to \code{permax}
}
\item{nl}{
The \code{nl} most significant genes in the lower tail will be plotted
}
\item{nr}{
The \code{nr} most significant genes in the upper tail will be plotted
}
\item{logs}{
If \code{logs=TRUE}, then log values are used.
}
\item{ig1}{
The columns of \code{data} for cases in group 1 (see Details)
}
\item{ig2}{
The columns in group 2.  By default, all the columns not in group 1.
}
\item{clmn.lab}{
  Labels for the columns in the image plot.
}
\item{row.lab}{
    Labels for the rows in the image plot.  (Note: the default is assumed to
    be evaluated after \code{data} has been subset to the rows to be plotted.)
}
\item{clmn.off}{
  Offset for printing the column labels (<0 to put labels
  outside the plot; larger negative numbers move the labels further from
  the axis).
}
\item{row.off}{
  Offset for printing the row labels (<0 to put labels outside
  the plot; larger negative numbers move the labels further from
  the axis). 
}
\item{...}{
Additional arguments to \code{\link{image}} (see \code{\link{par}})
}}
\section{Side Effects}{
An image plot is created on the current graphics device
}
\details{

Values within a row of \code{data} are centered and normalized to have 
variance 1.  If \code{ig1} is not given, then the values are centered to have
mean 0.  If \code{ig1} is given, the values are centered so the means of
the columns in \code{ig1} and \code{ig2} are equal in magnitude and
opposite in direction (usually \code{ig1} and \code{ig2} should match
the values used in the \code{\link{permax}} call).  The plot is thus
useful for comparing within rows, but differences in colors between rows
have no meaning.

The plot will give the most significant lower tail genes in the top
portion (most significant at the top), and the most significant upper
tail genes in the bottom portion (most significant at the bottom). 

This function just selects out the appropriate rows of \code{data}, and
calls \code{\link{plot.expr}()}.  \code{row.names(data)} or
\code{dimnames(data)[[1]]} must correspond to
\code{row.names(x)} for the selection to work properly.

}
\seealso{
\code{\link{permax}}, \code{\link{plot.expr}}, \code{\link{image}}
}
\examples{
   set.seed(1292)
   ngenes <- 1000
   m1 <- rnorm(ngenes,4,1)
   m2 <- rnorm(ngenes,4,1)
   exp1 <- cbind(matrix(exp(rnorm(ngenes*5,m1,1)),nrow=ngenes),
               matrix(exp(rnorm(ngenes*10,m2,1)),nrow=ngenes))
   exp1[exp1<20] <- 20
   sub <- exp1>20 & exp1<150
   exp1[sub] <- ifelse(runif(length(sub[sub]))<.5,20,exp1[sub])
   dimnames(exp1) <- list(paste('x',format(1:ngenes,justify='l'),sep=''),
                     paste('sample',format(1:ncol(exp1),justify='l'),sep=''))
   dimnames(exp1) <- list(paste('x',1:ngenes,sep=''),
                     paste('sample',1:ncol(exp1),sep=''))
   exp1 <- round(exp1)

   uu <- permax(exp1,1:5)
   plot(uu,exp1,ig1=1:5)
}
\keyword{hplot}

% Converted by Sd2Rd version 1.19.
