\name{permcor}
\alias{permcor}
\title{
permutation tests for correlations in high dimensional data
}
\description{
For high dimensional vectors of observations, computes the correlation
coefficient for each attribute with a specified vector of values, and
assesses significance using the permutation distribution of the maximum
and minimum over all attributes.
}
\usage{
permcor(data, phen, nperm=1000, logs=TRUE, ranks=FALSE, min.np=1, WHseed=NULL)
}
\arguments{
\item{data}{
Data matrix or data frame.  Each case is a column, and each row is an
attribute (the opposite of the standard configuration).
}
\item{phen}{
A vector of values (the ideal phenotype pattern).  The correlations of
each row of data with phen will be computed.
}
\item{nperm}{
The number of random permutations to use in computing the p-values.
The default is 1000.  If nperm is < 0, the entire permutation
distribution will be used, which is only feasible if the sample size is
fairly small 
}
\item{logs}{
If logs=TRUE (the default), then logs of the values in data are used in
computing correlations (the actual values of phen are used, though).
}
\item{ranks}{
If ranks=T, then within row ranks are used in place of the values in
data in the correlations.  The actual values of phen are still used.  
Default is ranks=FALSE.
}
\item{min.np}{
data will be subset to only rows with at least min.np values larger
than min(data).
}
\item{WHseed}{
Initial random number seed (a vector of 3 integers).  If missing, an
initial seed is generated from the runif() function.  Not needed if all
permutations are calculated.
}
}
\value{
Output is a data.frame of class c('permcor','permax'), with columns
  stat: the Pearson correlation coeffcients for each row of data
  pind: individual permutation p-values (2-sided)
  p2: 2-sided p-value using the distribution of the max overall rows
  p.lower: 1-sided p-value for lower levels in group 1
  p.upper: 1-sided p-value for higher levels in group 1
  nml: \# permutations where this row was the most significant for p.lower
  nmr: \# permutations where this row was the most sig for p.upper
  np: \# values > min(data) in each row


Also, if nperm>0, then output includes attributes 'seed.start' giving
the initial random number seed, and 'seed.end' giving the value of the
seed at the end.  These can be accessed with the attributes() and
attr() functions.
}

\details{
For DNA array data, this function is designed to identify the genes
with the largest positive and negative correlations with the phenotype
in phen.  Upper and lower p-values (p.upper, p.lower) are computed by
comparing each correlation to the permutation distribution of the maximum
and minimum (largest negative) correlations over all genes.  The \code{pind}
component of the output gives the p-value for the permutation
distribution of each individual gene.


If phen is a vector of 1's for the columns in group 1 and 0's for the
other columns, then the p-values from permcor() should be the same as
from permax() (to within simulation precision if random permutations
are used).  permax() is substantially more efficient in this setting.


The functions summary.permax() and plot.permax() can be used with the
output of permcor().


It is strongly recommended that different seeds be used for different
runs, and ideally the final seed from one run, attr(output,'seed.end'),
would be used as the initial seed in the next run.
}
\seealso{
permax, summary.permax, plot.permax.
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

#see the permax help file for the definition of exp1
  u8 <- permcor(exp1,1:15)
 summary(u8,nr=4,nl=4)
 u10 <- permcor(exp1[,c(1:3,5:8)],c(1,1,1,0,0,0,0),nperm=0)


 summary(u10,nl=4,nr=4)

}
\keyword{htest}

% Converted by Sd2Rd version 1.19.
