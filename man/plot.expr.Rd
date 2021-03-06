\name{plot.expr}
\alias{plot.expr}
\title{
Color image plot of gene expression levels
}
\description{
Represents values in the rows of a matrix as colored rectangles in an
image plot
}
\usage{
plot.expr(x, logs=TRUE, ig1=NULL, ig2=NULL, clmn.lab=dimnames(x)[[2]], 
    row.lab=dimnames(x)[[1]], clmn.off=NULL, row.off=NULL,\dots)
}
\arguments{
\item{x}{
matrix or \code{data.frame} containing the values to be plotted.
}
\item{logs}{
If \code{logs=TRUE}, then log values are used.
}
\item{ig1}{
The columns of \code{x} for cases in group 1 (see Details)
}
\item{ig2}{
The columns in group 2.  By default, all the columns not in group 1.
}
\item{clmn.lab}{
Labels for the columns in the array.
}
\item{row.lab}{
Labels for the rows in the array.
}
\item{clmn.off}{
Offset for printing the column labels (<0 to put labels
outside the plot).
}
\item{row.off}{
Offset for printing the row labels (<0 to put labels
outside the plot).
}
\item{...}{
Additional arguments to \code{image} and \code{text} (see \code{par})
none
}}
\section{Side Effects}{
An image plot is created on the current graphics device
}
\details{
Values within a row are centered and normalized to have variance 1.  If
\code{ig1} is not given, then the values are centered to have mean 0.
If \code{ig1} 
is given, the values are centered so the means of the columns in
\code{ig1} 
and \code{ig2} are equal in magnitude and opposite in direction.  The
plot is 
thus useful for comparing within rows, but differences in colors
between rows have no meaning.


A graphics device supporting image plots must be initialized prior to
calling this function.  Under Splus 3.4 for unix, the following command
(without the line breaks) initializes the X window motif plot window to
use 30 colors from blue (lowest levels) to yellow (highest levels) for
the image plots (in this scheme a value half way between the lowest and
highest values would be a medium intensity gray).


  motif("-xrm 'sgraphMotif.colorSchemes : background : black; lines :
  yellow cyan magenta green MediumBlue red; text : white yellow cyan
  magenta green MediumBlue red; images : blue 30 yellow'")


}
\seealso{
\code{plot.permax}
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

  plot.expr(exp1[1:20,])
}
\keyword{hplot}

% Converted by Sd2Rd version 1.19.
