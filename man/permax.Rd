\name{permax}
\alias{permax}
\title{
2-sample permutation tests for high dimensional data
}
\description{
For high dimensional vectors of observations, computes statistics for
each attribute, and assesses significance using the permutation
distribution of the maximum and minimum statistics over all attributes.
}
\usage{
permax(data, ig1, nperm=0, logs=TRUE, ranks=FALSE, min.np=1, ig2=NULL,
       WHseed=NULL, cluster=NULL, stratify=FALSE, weights=NULL,
       nl=50, nr=50, permute.cluster=FALSE, signed.rank=FALSE)
}
\arguments{
\item{data}{
Data matrix or data frame.  Each case is a column, and each row is an
attribute (the opposite of the standard configuration).
}
\item{ig1}{
The column numbers in \code{data} corresponding to group 1
}
\item{nperm}{
The number of random combinations to use in computing the p-values.
If \code{nperm <= 0}, then all possible combinations are enumerated, which is
only feasible if the sample sizes are fairly small
}
\item{logs}{
If \code{logs=TRUE} (the default), then logs of the values in
\code{data} are used in the statistics.
}
\item{ranks}{
If \code{ranks=TRUE}, then within row ranks are used in place of the values 
in \code{data} in the statistics.  This is equivalent to using the Wilcoxon
statistic.  Default is \code{ranks=FALSE}.
}
\item{min.np}{
\code{data} will be subset to only rows with at least \code{min.np}
values larger than \code{min(data)} in the columns in \code{ig1} and
\code{ig2} 
}
\item{ig2}{
The column numbers in \code{data} corresponding to group 2.  If
\code{NULL}, all columns not in \code{ig1} are assumed to be in group 2.
When both \code{ig1} and \code{ig2} are given, columns not in either are
excluded from the tests. 
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
If \code{TRUE}, statistics will be computed with each cluster, and
combined into an overall statistic using a weighted average.  Only
applicable if \code{cluster} is given.
}
\item{weights}{Weights for combining statistics over clusters.  If
  \code{NULL}, defaults to equal weight for all clusters.  Only
  applicable if \code{stratify=TRUE}.  
}
\item{nl}{A lower tail critical value is determined using the
  \code{nl}th largest value of the test statistics.  The number of genes
  with statistics <= this critical value is determined for each
  combination. 
}
\item{nr}{An upper tail critical value is determined using the
  \code{n-nr+1}st largest value of the test statistics, where \code{n}
  is the number of genes.  The number of genes
  with statistics >= this critical value is determined for each
  combination. 
}
\item{permute.cluster}{If \code{FALSE}, columns are permuted between
  groups within
  clusters.  If \code{TRUE}, clusters are permuted between groups.  The latter
  is only allowed if all columns in each cluster are in the same group.  Only
  applicable if \code{cluster} is given.  
}
\item{signed.rank}{If \code{TRUE}, the Wilcoxon signed rank statistic is
  computed.  This option can only be used if cluster is given, all
  cluster sizes = 2 and each cluster contains both groups (that is, only
  if there is a paired data structure).  Also, requires
  \code{stratify=FALSE}.  
}
}
\value{
Output is a data frame of class \code{permax}, with columns
  \item{$stat}{the test statistics for each row}
  \item{$pind.lower}{1-sided individual p-values for lower levels in group 1}
  \item{$pind.upper}{1-sided individual p-values for higher levels in group 1}
  \item{$p.lower}{1-sided simultaneous p-values for lower levels in group 1}
  \item{$p.upper}{1-sided simultaneous p-values for higher levels in group 1}
  \item{$m1, m2}{means of groups 1 and 2 (means of logs if \code{logs=TRUE})}
  \item{$s1, s2}{std deviations of groups 1 and 2 (of logs if \code{logs=TRUE})}
  \item{$np1,np2}{\# values \code{> min(data)} in groups 1 and 2}
  \item{$mdiff}{difference of means (if \code{logs=TRUE} the difference of geometric means)}
  \item{$mrat}{ratio of means (if \code{logs=TRUE}, ratio of geometric means)}

  NOTE: the \code{m1}, \code{m2}, \code{s1}, \code{s2}, \code{mdiff},
  and \code{mrat}, components of the output are computed ignoring
  \code{cluster}, even when \code{stratify=TRUE}.

The output also includes an attribute \code{call} giving the call to
\code{permax}, and an attribute \code{dist}, which is a vector of values
described below.  If \code{nperm>0}, then the output includes
attributes \code{seed.start} giving the initial random number seed, and
\code{seed.end} giving the value of the seed at the end.  These can be
accessed with the \code{\link{attributes}()} and \code{\link{attr}()}
functions.

The \code{dist} attribute gives a summary of the null permutation
distribution of the number of `positives'.  The \code{nl} and \code{nr}
arguments are used to determine critical values as the \code{nl} and
\code{n-nr+1} largest values of the test statistic computed from the
original data.  The number of genes <= the lower tail critical value and
the number >= the upper tail critical value are determined for each
permutation.  These are the `positives' in each direction.  The
\code{prop.nl} and \code{prop.nr} components give the proportion of
permutations with >= \code{nl} lower positives and >=\code{nr} upper
positives, the \code{prop.1l} and \code{prop.1r} components give the
proportion of permutations with any lower positives and any upper
positives, and the \code{ave.l} and \code{ave.r} components give the
average number of lower and upper positives per permutation.  The
\code{ave.l} and \code{ave.r} components thus give an estimate of the
expected number of false positives if the global null hypothesis of no
differences is true.

}
\details{
For DNA array data, this function is designed to identify the genes
which best discriminate between two tissue types.  2-sample t statistics
are computed for each gene using logs (default), raw values, or ranks.
As of version 2.2, upper and lower simultaneous p-values (\code{p.upper},
\code{p.lower}) are computed using the free step-down algorithm of
Westfall and Young (1993, Resampling-Based Multiple Testing, Section
2.6), applied to the test statistics, and thus may be slightly smaller 
than those obtained with earlier versions, which compared each statistic
to the permutation distribution of the maximum and minimum (largest
negative) statistics over all genes.  These p-values are corrected for
simultaneous inference (and control the familywise error rate), and can
be used to test the global null of no difference for any of the
attributes in \code{data}.  The \code{pind.lower} and
\code{pind.upper} components of the output give the p-values for the
permutation distribution of each individual gene (with no correction for
multiple comparisons).

If \code{cluster=NULL}, then the columns of \code{data} are assumed to
be independent.  When \code{cluster} is specified, there are 3 options.  If
\code{stratify=FALSE} and \code{permute.cluster=FALSE}, then statistics
are computed as in the independence case, but columns are only permuted
within clusters, since under the null it is assumed that columns from
the same cluster are exchangeable, while columns from different clusters
may not be.  If \code{stratify=TRUE} and \code{permute.cluster=FALSE},
then statistics are also stratified on \code{cluster}.  In this case,
statistics are computed within the strata and combined with a weighted
average over strata.  The weights are given by the \code{weights}
argument.  This might be more appropriate if clusters are unbalanced or
if cluster effects are very large.  If \code{permute.cluster=TRUE}, then
it is assumed that group or treatment is defined at the cluster level
(that is, all members of each cluster are in the same group).  In this
case statistics are computed the same as in the independence case, but
only whole clusters are permuted to preserve the within cluster
dependencies.  If the data come from a split plot experiment, then the
\code{stratify=TRUE} option would be appropriate for the subplot
effects, while the \code{permute.cluster=TRUE} option would examine the
whole plot effects.  The combination \code{stratify=TRUE} and
\code{permute.cluster=TRUE} is not allowed.

When \code{ranks=FALSE} and \code{stratify=TRUE}, the data values are first
standardized to have mean 0 and variance 1 within each stratum for each
row.  The statistics computed are the difference of group means of the
standardized values within each stratum.  When \code{ranks=TRUE} and
\code{stratify=TRUE}, the values are ranked separately within each
stratum, and the statistics computed are the difference of group means
of the ranks within each stratum.  

The test statistic returned is the standard t statistic if \code{ranks=FALSE}
and \code{stratify=FALSE}.  It is the sum of the mean corrected ranks in
group 1 if \code{ranks=TRUE} and \code{stratify=FALSE}.  It is the
weighted average of the within stratum differences of means of the
standardized values if \code{ranks=FALSE} and \code{stratify=TRUE}, and
the weighted average of the within stratum differences of the
ranks if \code{ranks=TRUE} and \code{stratify=TRUE}.  See the document
Permax Computations for additional details (a .pdf copy should be in the
\code{doc} subdirectory of the \code{permax} library).

In the special case where each cluster consists of a pair of objects,
one from each group, the Wilcoxon signed-rank statistic can be used by
specifying \code{signed.rank=TRUE}.  In this case, within pair
differences of expression levels are taken, and the absolute values are
ranked.  The statistic computed is then the sum of the signed ranks for
all pairs, with pairs with a difference of 0 given a rank of 0.

In this function, the permutation distribution is generated by selecting
combinations of columns to be in group 1, rather than by generating full
permutations.

With random combinations, it is strongly recommended that different
seeds be used for different runs, and ideally the final seed from one run,
\code{\link{attr}(output,'seed.end')}, would be used as the initial seed
in the next run. 
}
\seealso{
\code{\link{summary.permax}}, \code{\link{plot.permax}}, \code{\link{permcor}}, \code{\link{permsep}}.
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
# columns 1,2,6:9 from cluster 1, columns 3,4,10:13 from cluster 2,
# columns 5,14,15 from cluster 3
   clustind <- c(1,1,2,2,3,1,1,1,1,2,2,2,2,3,3)
   uuc <- permax(exp1,1:5,nperm=5000,cluster=clustind)
   summary(uuc,nl=5,nr=5) 
# stratified on cluster, using ranks; equal weights
   uus <- permax(exp1,1:5,nperm=5000,cluster=clustind,stratify=TRUE,
            ranks=TRUE,WHseed=attr(uuc,'seed.end'))
   summary(uus,nl=5,nr=5) 
# stratified on cluster, using ranks; average weighted by cluster sizes
   uus <- permax(exp1,1:5,nperm=5000,cluster=clustind,stratify=TRUE,
            ranks=TRUE,weights=table(clustind),WHseed=attr(uus,'seed.end'))
   summary(uus,nl=5,nr=5)
# permute clusters
   clust2 <- c(1,1,2,2,2,3,3,4,4,5,5,5,6,6,6)
   uu5 <- permax(exp1,1:5,nperm=5000,cluster=clust2,permute.cluster=TRUE,
            ranks=TRUE,WHseed=attr(uus,'seed.end'))
   summary(uu5,nl=5,nr=5)
# paired data
   clust <- c(1,2,3,4,1,2,3,4)
   dat <- rbind(c(1:8),c(8,1:7),c(5,1:7),c(8,7,1:6))
   dimnames(dat) <- list(1:4,1:8)
   permax(dat,1:4,cluster=clust,signed.rank=TRUE)
}
\keyword{htest}

% Converted by Sd2Rd version 1.19.
