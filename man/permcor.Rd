\name{permcor}
\alias{permcor}
\title{
permutation tests for correlations in high dimensional data
}
\description{
For high dimensional vectors of observations, computes the correlation
coefficient for each attribute with a specified vector of values, and
assesses significance using the permutation distribution of the maximum
and minimum correlation over all attributes.
}
\usage{
permcor(data, phen, nperm=1000, logs=TRUE, ranks=FALSE, min.np=1,
        WHseed=NULL, cluster=NULL, stratify=FALSE, weights=NULL,
        permute.cluster=FALSE)
}
\arguments{
\item{data}{
Data matrix or data frame.  Each case is a column, and each row is an
attribute (the opposite of the standard configuration).
}
\item{phen}{
A vector of values (the ideal phenotype pattern or covariate).  The
correlations of each row of \code{data} with \code{phen} will be computed.
}
\item{nperm}{
The number of random permutations to use in computing the p-values.
The default is 1000.  If \code{nperm} is <= 0, the entire permutation
distribution will be used, which is only feasible if the sample size is
fairly small.
}
\item{logs}{
If \code{logs=TRUE} (the default), then logs of the values in
\code{data} are used in computing correlations (the actual values of
\code{phen} are used, though). 
}
\item{ranks}{
If \code{ranks=TRUE}, then within row ranks are used in place of the values
in \code{data} in the correlations.  The actual values of \code{phen}
are still used.  Default is \code{ranks=FALSE}.
}
\item{min.np}{
\code{data} will be subset to only rows with at least \code{min.np}
values larger than \code{min(data)}.
}
\item{WHseed}{
Initial random number seed (a vector of 3 integers).  If \code{NULL}, an
initial seed is generated from the \code{\link{runif}()} function.  Not
used if \code{nperm <= 0}.
}
\item{cluster}{
  A vector of cluster membership indicators for the columns of
  \code{data}.  Distinct values of \code{cluster} define distinct
  clusters or blocks.  
}
\item{stratify}{
If \code{TRUE}, correlations will be estimated within each cluster, and
combined into an overall estimate using a weighted average.  Only
applicable if \code{cluster} is given.
}
\item{weights}{Weights for combining statistics over clusters.  If
  \code{NULL}, defaults to
  equal weight for all clusters.  Only applicable if \code{stratify=TRUE}.  
}
\item{permute.cluster}{If \code{FALSE}, values of \code{phen} are
  permuted among columns within
  clusters.  If \code{TRUE}, values are permuted among clusters.  The latter
  is only allowed if all columns in each cluster have the same value of
  \code{phen}.  Only
  applicable if \code{cluster} is given.  
}
}
\value{
Output is a data frame of class \code{c('permcor','permax')}, with columns
  \item{$stat}{the estimated correlation coefficients for each row of \code{data}}
  \item{$pind.lower}{1-sided individual p-values for high negative correlations}
  \item{$pind.upper}{1-sided individual p-values for high positive correlations}
  \item{$p.lower}{1-sided simultaneous p-values for high negative correlations}
  \item{$p.upper}{1-sided simultaneous p-values for high positive correlations}
  \item{$np}{\# values > \code{min(data)} in each row}


The output also includes an attribute \code{call} giving the call to
\code{permcor}, and if \code{nperm>0}, the output includes attributes
\code{seed.start} giving the initial random number seed, and
\code{seed.end} giving the value of the
seed at the end.  These can be accessed with the
\code{\link{attributes}()} and \code{\link{attr}()} functions.
}


\details{

For DNA array data, this function is designed to identify the genes with
the largest positive and negative correlations with the phenotype in
\code{phen}.  As of version 2.2, upper and lower simultaneous p-values
(\code{p.upper}, \code{p.lower}) are computed using the free step-down
algorithm of Westfall and Young (1993, Resampling-Based Multiple
Testing, Section 2.6), applied to the correlations, and thus may be
slightly smaller than those obtained with earlier versions, which
compared each correlation to the permutation distribution of the maximum
and minimum (largest negative) correlations over all genes.  These
p-values are corrected for simultaneous inference (and control the
familywise error rate), and can be used to test the global null of no
correlation with \code{phen} for any of the attributes in \code{data}.
The \code{pind.lower} and \code{pind.upper} components of the output
gives the p-value for the permutation distribution of each individual gene.

If \code{phen} is a vector of 1's for the columns in group 1 and 0's for the
other columns, then the p-values from \code{permcor()} should be the same as
from \code{\link{permax}()} (to within simulation precision if random
permutations are used).  \code{\link{permax}()} is more efficient in
this setting. 

If \code{cluster=NULL}, then the columns of \code{data} are assumed to
be independent.  When \code{cluster} is specified, there are 3 options.  If
\code{stratify=FALSE} and \code{permute.cluster=FALSE}, then
correlations are computed as in the independence case, but columns are
only permuted within clusters, since under the null it is assumed that
columns from the same cluster are exchangeable, while columns from
different clusters may not be.  If \code{stratify=TRUE} and
\code{permute.cluster=FALSE}, then the analysis is also stratified on
\code{cluster}.  In this case correlations are computed separately
within the strata and combined with a weighted average over strata.  The
weights are given by the \code{weights} argument (\code{weights} will be
normalized to sum to 1).  This might be more appropriate if cluster
effects are very large.  If \code{permute.cluster=TRUE}, then it is
assumed that phenotype is defined at the cluster level (that is, all
members of each cluster have the same value of \code{phen}).  In this case
correlations are computed the same as in the independence case, but
phenotype is permuted among the whole clusters, rather then among
individual columns, to preserve the within cluster dependencies.  The
combination \code{stratify=TRUE} and \code{permute.cluster=TRUE} is not
allowed.

When \code{permute.cluster=TRUE} and cluster sizes vary, the
distribution of \code{phen} will vary over the permutations, since the
number of columns with a given values of \code{phen} will depend on
which cluster that value is assigned to.

When \code{ranks=TRUE} and \code{stratify=TRUE}, the values in \code{data}
are ranked separately within each stratum.

The functions \code{\link{summary.permax}()} and
\code{\link{plot.permax}()} can be used with the 
output of \code{permcor()}.

With random permutations, it is strongly recommended that different
seeds be used for different runs, and ideally the final seed from one run,
\code{\link{attr}(output,'seed.end')}, would be used as the initial seed in the
next run. 
}
\seealso{
\code{\link{permax}}, \code{\link{summary.permax}}, \code{\link{plot.permax}}.
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

   u8 <- permcor(exp1,1:15)
   summary(u8,nr=4,nl=4)
   u10 <- permcor(exp1[,c(1:3,5:8)],c(1,1,1,0,0,0,0),nperm=0)
   summary(u10,nl=4,nr=4)
# columns 1,2,6:9 from cluster 1, columns 3,4,10:13 from cluster 2,
# columns 5,14,15 from cluster 3
   clustind <- c(1,1,2,2,3,1,1,1,1,2,2,2,2,3,3)
   u8 <- permcor(exp1,1:15,cluster=clustind,WHseed=attr(u8,'seed.end'))
   summary(u8,nr=4,nl=4)
# correlations estimated within clusters; average weighted by cluster sizes
   u8 <- permcor(exp1,1:15,cluster=clustind,stratify=TRUE,
                 weights=table(clustind),WHseed=attr(u8,'seed.end'))
   summary(u8,nr=4,nl=4)
# permuting clusters: phenotype is a characteristic of the cluster
   ucorx <- permcor(exp1[,3:10],c(1,2,2,3,3,4,4,4),nperm=0,
      cluster=c(1,2,2,3,3,4,4,4),permute.cluster=TRUE)
   summary(ucorx,nl=5,nr=5)
}
\keyword{htest}

% Converted by Sd2Rd version 1.19.
