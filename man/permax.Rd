\name{permax}
\alias{permax}
\title{
2-sample permutation t-tests for high dimensional data
}
\description{
For high dimensional vectors of observations, computes t statistics for
each attribute, and assesses significance using the permutation
distribution of the maximum and minimum over all attributes.
}
\usage{
permax(data, ig1, nperm=0, logs=T, ranks=F, min.np=1, ig2, WHseed=NULL)
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
\item{logs}{
If logs=T (the default), then logs of the values in data are used in
the statistics.
}
\item{ranks}{
If ranks=T, then within row ranks are used in place of the values in
data in the t statistics.  This is equivalent to using the Wilcoxon
statistic.  Default is ranks=F
}
\item{min.np}{
data will be subset to only rows with at least min.np values larger
than min(data) in the columns in ig1 and ig2 
}
\item{ig2}{
The columns of data corresponding to group 2.  The default is to
include all columns not in ig1 in group 2.  When both ig1 and ig2 are
given, columns not in either are excluded from the tests.
}
\item{WHseed}{
Initial random number seed (a vector of 3 integers).  If missing, an
initial seed is generated from the runif() function.  Not needed if all
permutations are calculated.
Output is a data.frame of class 'permax', with columns
  stat: the standardized test statistics for each row
  pind: individual permutation p-values (2-sided)
  p2: 2-sided p-value using the distribution of the max overall rows
  p.lower: 1-sided p-value for lower levels in group 1
  p.upper: 1-sided p-value for higher levels in group 1
  nml: \# permutations where this row was the most significant for p.lower
  nmr: \# permutations where this row was the most sig for p.upper
  m1, m2: means of groups 1 and 2 (means of logs if logs=T)
  s1, s2: std deviations of groups 1 and 2 (of logs if logs=T)
  np1,np2: \# values > min(data) in groups 1 and 2
  mdiff: difference of means (if logs=T the difference of geometric means)
  mrat: ratio of means (if logs=T ratio of geometric means)


Also, if nperm>0, then output includes attributes 'seed.start' giving
the initial random number seed, and 'seed.end' giving the value of the
seed at the end.  These can be accessed with the attributes() and
attr() functions.


}}
\details{
For DNA array data, this function is designed to identify the genes
which best discriminate between two tissue types.  2-sample t
statistics are computed for each gene using logs (default), raw values,
or ranks.  Upper and lower p-values (p.upper, p.lower) are computed by
comparing each statistic to the permutation distribution of the maximum
and minimum (largest negative) statistic over all genes.  The \code{pind}
component of the output gives the p-value for the permutation
distribution of each individual gene.


It is strongly recommended that different seeds be used for different
runs, and ideally the final seed from one run, attr(output,'seed.end'),
would be used as the initial seed in the next run.
}
\seealso{
summary.permax, plot.permax, permcor, permsep.
}
\examples{
#generate make believe data
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
  summary(uu,nl=5,nr=5) # 5 most extreme in each direction
}
\keyword{htest}

% Converted by Sd2Rd version 1.19.
