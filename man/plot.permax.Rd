\name{plot.permax}
\alias{plot.permax}
\title{
Image plot of the most significant genes (attributes) from a permax()
analysis 
}
\description{
Given the output of permax(), and the array of expression levels,
creates a color image plot of the expression levels of the most
significant genes 
}
\usage{
plot.permax(Z, data, nl=25, nr=25, logs=T, ig1=NULL, ig2=NULL, \dots)
}
\arguments{
\item{Z}{
A permax object (output from permax())
}
\item{data}{
Matrix or data frame of expression levels used as input to permax()
}
\item{nl}{
The nl most significant genes in the lower tail will be plotted
}
\item{nr}{
The nr most significant genes in the upper tail will be plotted
}
\item{logs}{
If logs=T, then log values are used.
}
\item{ig1}{
The columns of data for cases in group 1 (see Details)
}
\item{ig2}{
The columns in group 2.  By default, all the columns not in group 1.
}
\item{...}{
Additional arguments to image() and text() (see par())
none
}}
\section{Side Effects}{
An image plot is created on the current graphics device
}
\details{
Values within a row of data are centered and normalized to have
variance 1.  If ig1 is not given, then the values are centered to have
mean 0.  If ig1 is given, the values are centered so the means of the
columns in ig1 and ig2 are equal in magnitude and opposite in direction
(usually ig1 and ig2 should match the values used in the permax()
call).  The plot is thus useful for comparing within rows, but
differences in colors between rows have no meaning.


The plot will give the most significant lower tail genes in the top
portion (most significant at the top), and the most significant upper
tail genes in the bottom portion (most significant at the bottom). 


This function just selects out the appropriate rows of data, and calls
plot.expr().  row.names(data) or dimnames(data)[[1]] must correspond to
row.names(Z) for the selection to work properly.


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
permax, plot.expr
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
  plot(uu,exp1,ig1=1:5,cex=.7)
}
\keyword{hplot}

% Converted by Sd2Rd version 1.19.
