\name{permsep}
\alias{permsep}
\title{
Permutation analysis for complete separation
}
\description{
Given two groups of samples of high dimensional attribute vectors (eg
DNA array expression levels), determines the number of attributes which
completely separate the two groups (all values in one group strictly
larger than in the other), and a permutation p-value for this quantity.
}
\usage{
permsep(data, ig1, nperm=0, ig2, WHseed=NULL)
}
\arguments{
\item{data}{
Data matrix or data frame.  Each case is a column, and each row is an
attribute (the opposite of the standard configuration).
}
\item{ig1}{
The columns of data corresponding to group 1
}
\item{nperm}{
The number of random permutations to use in computing the p-values.
The default is to use the entire permutation distribution, which is
only feasible if the sample sizes are fairly small
}
\item{ig2}{
The columns of data corresponding to group 2.  The default is to
include all columns not in ig1 in group 2.  When both ig1 and ig2 are
given, columns not in either are excluded.
}
\item{WHseed}{
Initial random number seed (a vector of 3 integers).  If missing, an
initial seed is generated from the runif() function.  Not needed if all
permutations are calculated.


Prints a vector giving the \# genes with complete separation (all in one group 
larger than all in the other, the proportion of permutations with
this many or more genes with complete separation (p-value) ('permutation'
actually means a distinct rearrangement of columns into 2 groups), the 
average number of genes per permutation with complete separation, 
and the proportion of permutations with any genes with complete
separation.  


The value returned is a list with components \code{ics} = a vector indicating
(with 1) which rows of data have complete separation, and \code{dtcs} = a
vector containing the printed output.


Also, if \code{nperm}>0, then the output includes attributes \code{seed.start}
giving the initial random number seed, and \code{seed.end} giving the value
of the seed at the end.  These can be accessed with the
\code{\link{attributes}} 
and \code{\link{attr}} functions.


}}
\details{
For each gene there will be 0, 1 or 2 rearrangements with complete
separation, depending on the number of unique values and the sizes of
the two groups.  Adding these numbers over genes and dividing by the
number of rearrangements gives the average number per permutation.  The
value returned averages only over the rearrangements actually used,
though.


It is strongly recommended that different seeds be used for different
runs, and ideally the final seed from one run, attr(output,'seed.end'),
would be used as the initial seed in the next run.


}
\seealso{
permax
}
\examples{
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

 uuu <- permsep(exp1,1:5)

}
\keyword{htest}
% Converted by Sd2Rd version 1.19.
