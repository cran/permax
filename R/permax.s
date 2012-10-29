# By R Gray, June 24, 2000, DFCI
# Copyright (C) 2000 Robert Gray
# Distributed under the GNU Public License (see the file COPYING)

permax <- function(data,ig1,nperm=0,logs=TRUE, ranks=FALSE, min.np=1,
                   ig2, WHseed=NULL) {
### data=data matrix; markers in rows, samples in columns, gene codes used for
### matching should be in dimnames(data)[[1]]
### ig1=columns of data in group 1
### nperm <= 0 => compute full permutation distribution
### nperm >0 => compute nperm random samples
### if logs=T then summary statistics are computed from logs of data, and logs
###   are used in the t statistics (if rank=F)
### if ranks=T then ranks are used in the t statistics (giving the Wilcoxon
###   test)
### min.np : data will be subset to only rows with at least min.np
###   values > min(data) in the columns in ig1 and ig2
### ig2 : column numbers for group 2.  If missing, all columns not in group
###   1 are assumed to be in group 2. min.np applies only to columns in ig1
###   and ig2
### WHseed = Initial random number seed (vector of 3 integers).  If missing,
###   generated from the runif() function.  Not needed if all permutations
###   are calculated.
### Output is a data.frame with columns
###   stat: the standardized test statistics
###   pind: individual permutation p-values (2-sided)
###   p2: 2-sided p-value using dist of max overall rows
###   p.lower: 1-sided p-value for lower levels in group 1
###   p.upper: 1-sided p-value for higher levels in group 1
###   nml: # permutations where this row was the most significant for p.lower
###   nmr: # permutations where this row was the most sig for p.upper
###   m1, m2: means of groups 1 and 2 (means of logs if logs=T)
###   s1, s2: std deviations of groups 1 and 2 (of logs if logs=T)
###   np1,np2: # pos (actually number > minimum value in data) in grps 1 & 2
###   mdiff: difference of means (if logs=T the diff of geometric means)
###   mrat: ratio of means (if logs=T ratio of geometric means)
  data <- as.matrix(data)
  if (logs) {
	tmp <- data<=0
	if(any(tmp))
	  data[tmp] <- 1
	data <- log(data)
   }
  dmin <- min(data)
  if (missing(ig2)) {
    data <- cbind(data[,ig1],data[,-ig1])
  } else {
    data <- cbind(data[,ig1],data[,ig2])
  }
  n1 <- length(ig1)
  n2 <- ncol(data)-n1
  ig1 <- 1:n1
  d1 <- data[,ig1,drop=FALSE]
  if (n1>1) {
    m1 <- c(d1 %*% rep(1/n1,n1))
    s1 <- sqrt((d1-m1)^2 %*% rep(1/(n1-1),n1))
    d1[d1<=dmin] <- 0
    d1[d1>dmin] <- 1
    npos1 <- d1 %*% rep(1,n1)
  } else {
    m1 <- c(d1)
    s1 <- rep(0,length(d1))
    npos1 <- ifelse(d1>dmin,1,0)
  }
  d1 <- data[,-ig1,drop=FALSE]
  m2 <- c(d1 %*% rep(1/n2,n2))
  s2 <- if(n2>1) sqrt((d1-m2)^2 %*% rep(1/(n2-1),n2)) else rep(0,nrow(d1))
  d1[d1<=dmin] <- 0
  d1[d1>dmin] <- 1
  npos2 <- d1 %*% rep(1,n2)
  sub <- npos1+npos2 >= min.np
  data <- cbind(data[sub,ig1],data[sub,-ig1])
  if (ranks) data <- t(apply(data,1,rank))
  n <- nrow(data)
  if (nperm>0) {
    if (is.null(WHseed)) WHseed <- floor(30000*runif(3))+1
  } else {
    WHseed <- c(0,0,0)
    nn <- exp(sum(log(2:(n1+n2)))-sum(log(2:n1))-sum(log(2:n2)))
    cat('statistics will be computed for all',format(nn),'groupings\n')
  }
  Z <- .Fortran('ptn',d=as.single(data),n=as.integer(n),ng=as.integer(
    ncol(data)),ng1=as.integer(length(ig1)),stat=single(n),dm=single(n),
    pind=integer(n),p2=integer(n),p.lower=integer(n),p.upper=integer(n),
    nperm=as.integer(nperm),ix=as.integer(WHseed),nml=integer(n),
    nmr=integer(n), PACKAGE="permax")[5:14]
  if (nperm>0) endseed <- Z$ix
  Z <- data.frame(stat=Z$stat,pind=Z$pind/Z$nperm,p2=Z$p2/Z$nperm,
    p.lower=Z$p.lower/Z$nperm,p.upper=Z$p.upper/Z$nperm,nml=Z$nml,nmr=Z$nmr)
  m1 <- m1[sub]
  m2 <- m2[sub]
  if (logs){
    Z <- cbind(Z,m1=m1,m2=m2,s1=s1[sub],s2=s2[sub],np1=npos1[sub],
      np2=npos2[sub],mdiff=exp(m1)-exp(m2),mrat=exp(m1-m2))
  } else {
    Z <- cbind(Z,m1=m1,m2=m2,s1=s1[sub],s2=s2[sub],np1=npos1[sub],
      np2=npos2[sub],mdiff=m1-m2,mrat=m1/m2)
  }
  row.names(Z) <- dimnames(data)[[1]]
  class (Z) <- c('permax','data.frame')
  if (nperm>0) {
    attr(Z,'seed.start') <- WHseed
    attr(Z,'seed.end') <- endseed
  }
  Z
}

summary.permax <- function(object, data, nl=25, nr=25, ...) {
### prints the nl most significant in the lower tail and the nr most
###   significant in the upper tail, and returns the row.names of object
###   corresponding to these two groups in a list
###   if data (matrix) is specified and contains rows with dimnames matching
###   the row.names
### object is a dataframe (output from permax), data is a matrix
### nl and nr are the number selected in the lower and upper tails
    ans<-list(nl=nl, nr=nr)
    o <- order(object$stat)
    ans$objectl <- object[o[1:nl],]
    if (!missing(data)) {
        data <- as.matrix(data)
        ans$dataL <- data[match(row.names(ans$objectl), dimnames(data)[[1]],0),]
    }
    ans$objectr <- object[rev(o)[1:nr],]
    if (!missing(data))
        ans$dataR <- data[match(row.names(ans$objectr),dimnames(data)[[1]],0),]
    ans$lower <- row.names(ans$objectl)
    ans$upper <- row.names(ans$objectr)
    class(ans) <- "summary.permax"
    ans
}

print.summary.permax <- function(x, digits = max(3,
                                    getOption("digits") - 3), ...)
{
    print(x$objectl[,-c(2,3,5,6,7)], digits=digits)
    if( !is.null(x$dataL) )
        print(x$dataL, digits=digits)
    print(x$objectr[,-c(2,3,4,6,7)], digits=digits)
    if( !is.null(x$dataR) )
        print(x$dataR, digits=digits)
}

permsep <- function(data,ig1,nperm=0,ig2,WHseed=NULL) {
### perm dist of # of genes with complete separation
### data=data matrix; markers in rows, samples in columns
### ig1=columns of data in group 1
### nperm <= 0 => compute full permutation distribution
### nperm >0 => compute nperm random samples
### ig2 : column numbers for group 2.  If missing, all columns not in group
###   1 are assumed to be in group 2.
### WHseed = Initial random number seed (vector of 3 integers).  If missing,
###   generated from the runif() function.  Not needed if all permutations
###   are calculated.
### Printed Output: # genes with complete separation (all in one group
###   larger than all in the other, the proportion of permutations with
###   this many or more genes with complete separation (p-value)('permutation'
###   actually means a distinct rearrangement of columns into 2 groups), the
###   average number of genes per permutation with complete separation,
###   and the proportion of permutations with any genes with complete
###   separation.           (Note: for each gene there
###   will be 0, 1 or 2 rearrangements with complete separation, depending
###   on the number of unique values and the sizes of the two groups.  Adding
###   these numbers over genes and dividing by the number of rearrangements
###   gives the average number per permutation.  The value returned averages
###   only over the rearrangements actually used, though.)
### value returned: list with a vector indicating (with 1) which rows of data
###   have complete separation as the first component, and a vector containing
###   the printed output as the second component
  data <- as.matrix(data)
  if (missing(ig2)) {
    data <- cbind(data[,ig1],data[,-ig1])
  } else {
    data <- cbind(data[,ig1],data[,ig2])
  }
  n <- nrow(data)
  n1 <- length(ig1)
  if (nperm>0) {
    if (is.null(WHseed)) WHseed <- floor(30000*runif(3))+1
  } else {
    WHseed <- c(0,0,0)
    nn <- exp(sum(log(2:(ncol(data))))-sum(log(2:n1))-sum(log(2:(ncol(data)-n1))))
    cat('statistics will be computed for all',format(nn),'groupings\n')
  }
  Z <- .Fortran('ptc',d=as.single(data),n=as.integer(n),ng=as.integer(
    ncol(data)),ng1=as.integer(length(ig1)),ics=integer(nrow(data)),
    nperm=as.integer(nperm),dtcs=integer(4),ix=as.integer(WHseed),
                PACKAGE="permax")[5:8]
  if (nperm>0) endseed <- Z$ix
  Z$dtcs <- c(Z$dtcs[1],Z$dtcs[2:4]/Z$nperm)
  cat('# comp sep; prop perm with more; ave # per perm; prop perm with any\n')
  print(Z$dtcs)
  Z <- Z[-c(2,4)]
  names(Z$dtcs) <- c('obs.num','pval.no','avenum.per.perm','prop.with.any')
  names(Z$ics) <- dimnames(data)[[1]]
  if (nperm>0) {
    attr(Z,'seed.start') <- WHseed
    attr(Z,'seed.end') <- endseed
  }
  Z
}

rowperm <- function(x) {
# applies a separate permutation
# to the elements in each row of the input array x
  x <- as.matrix(x)
  dimx <- dim(x)
  dimnmx <- dimnames(x)
  x <- .Fortran('pa',as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),
     as.double(runif(nrow(x)*(ncol(x)-1))), PACKAGE="permax")[[1]]
  dim(x) <- dimx
  dimnames(x) <- dimnmx
  x
}

#plot.permax <- function(x, data, nl=25, nr=25, logs=TRUE, ig1=NULL,
#                        ig2=NULL, ...) {
# plots the expression levels for the most significant permax genes
# x = output from permax
# data = expression level array used as input to permax
# nl, nr = # of most extreme genes in lower (nl) and upper (nr) tails to plot
# logs = if true, function takes logs of values in data
# if ig1 is given a non-null value (which must be a vector of integers
#   corresponding to columns in data, then the rows of data will be
#   standardized so the mean of the columns in ig1 and the mean of the
#   columns in ig2 are equal in magnitude and opposite in sign (if ig2 is
#   not specified, it defaults to include all the columns not in ig1)
#   if ig1 is NULL, then the rows of data are standardized to have mean 0
#   In either case, the rows are also standardized to have variance 1
# A graphics device with support for appropriate image colors must be
#   specified prior to calling this function
#  o <- order(x$stat)
#  xl <- row.names(x)[o[1:nl]]
#  xr <- rev(row.names(x)[rev(o)[1:nr]])
#  plot.expr(data[c(xl,xr),],logs=logs,ig1,ig2,...)
#}

#plot.expr <- function(x, logs=TRUE, ig1=NULL, ig2=NULL, ...) {
## plots the expression levels for the rows of the arrary data
# x = expression level array
# logs = if true, function takes logs of values in data
# if ig1 is given a non-null value (which must be a vector of integers
#   corresponding to columns in data, then the rows of data will be
#   standardized so the mean of the columns in ig1 and the mean of the
#   columns in ig2 are equal in magnitude and opposite in sign (if ig2 is
#   not specified, it defaults to include all the columns not in ig1)
#   if ig1 is NULL, then the rows of data are standardized to have mean 0
#   In either case, the rows are also standardized to have variance 1
# A graphics device with support for appropriate image colors must be
#   specified prior to calling this function
# All rows of data are plotted in the order given, so appropriate sorting
#   and subsetting should be done prior to the call.
#  x <- as.matrix(x)
#  if (logs) x <- log(x)
#  m1 <- apply(x,1,mean)
#  x <- x-m1
#  s1 <- sqrt(apply(x^2,1,sum)/(ncol(x)-1))
#  x <- x/s1
## with unequal group sizes, standardization as above tends to guarantee
# that the smaller group takes more extreme values than the larger.
# the following shifts the overall mean so the mean of each group is
# equal in magnitude and opposite in sign.
#  if (!is.null(ig1)) {
#    if (is.null(ig2)) ig2 <- (1:ncol(x))[-ig1]
#    m1 <- apply(x[,ig1],1,mean)
#    m2 <- apply(x[,ig2],1,mean)
#    x <- x-(m1+m2)/2
#  }
#  image(1:ncol(x),1:nrow(x),t(x),xaxt='n',yaxt='n', xlab="",
#        ylab="",...)
#  text(-(ncol(x)/35), 1:nrow(x), dimnames(x)[[1]],xpd = TRUE, ...)
#  text(1:ncol(x),-(nrow(x)/9),dimnames(x)[[2]],srt=270,xpd=TRUE,...)
##  positioning is different in graphsheet() plots
#  text(1:ncol(x),-(nrow(x)/7),dimnames(x)[[2]],srt=270,xpd=TRUE,...)
#  invisible()
#}

permcor <- function(data, phen, nperm=1000, logs=TRUE, ranks=FALSE,
                    min.np=1, WHseed=NULL) {
### data=data matrix; markers in rows, samples in columns, gene codes used for
### matching should be in dimnames(data)[[1]]
### phen=vector of length ncol(data) giving the target attributes
###   (phenotype).  The rows most positively and negatively correlated
###   with phen are identified
### nperm <= 0 => compute full permutation distribution
### nperm >0 => compute nperm random permutation
### if logs=T then logs of values in data are used (if rank=F)
### if ranks=T then correlations are computed from within row ranks.
###   test)
### (In any case, the actualy values of phen are used)
### min.np : data will be subset to only rows with at least min.np
###   values > min(data)
### WHseed = Initial random number seed (vector of 3 integers).  If missing,
###   generated from the runif() function.  Not needed if all permutations
###   are calculated.
### Output is a data.frame with columns
###   stat: the correlation coefficients
###   pind: individual permutation p-values (2-sided)
###   p2: 2-sided p-value using dist of max overall rows
###   p.lower: 1-sided p-value for lower levels in group 1
###   p.upper: 1-sided p-value for higher levels in group 1
###   nml: # permutations where this row was the most significant for p.lower
###   nmr: # permutations where this row was the most sig for p.upper
###   np: # pos (actually number > minimum value in data) in each row
  data <- as.matrix(data)
  if (logs) data <- log(data)
  dmin <- min(data)
  ng <- ncol(data)
  d1 <- data-dmin
#  d1[d1<=dmin] <- 0
  d1[d1>0] <- 1
  npos <- d1 %*% rep(1,ng)
  data <- data[npos>=min.np,]
  if (ranks) data <- t(apply(data,1,rank))
  n <- nrow(data)
  if (nperm>0) {
    if (is.null(WHseed)) WHseed <- floor(30000*runif(3))+1
  } else {
    WHseed <- c(0,0,0)
    nn <- exp(sum(log(2:(ng))))
    cat('statistics will be computed for all',format(nn),'permutationss\n')
  }
  Z <- .Fortran('ptcor',d=as.single(data),n=as.integer(n),ng=as.integer(ng),
    phen=as.single(phen),stat=single(n),
    pind=integer(n),p2=integer(n),p.lower=integer(n),p.upper=integer(n),
    nperm=as.integer(nperm),ix=as.integer(WHseed),nml=integer(n),
    nmr=integer(n), PACKAGE="permax")[5:13]
  if (nperm>0) endseed <- Z$ix
  Z <- data.frame(stat=Z$stat,pind=Z$pind/Z$nperm,p2=Z$p2/Z$nperm,
    p.lower=Z$p.lower/Z$nperm,p.upper=Z$p.upper/Z$nperm,nml=Z$nml,
    nmr=Z$nmr,np=npos[npos>=min.np])
  row.names(Z) <- dimnames(data)[[1]]
  class (Z) <- c('permcor','permax','data.frame')
  if (nperm>0) {
    attr(Z,'seed.start') <- WHseed
    attr(Z,'seed.end') <- endseed
  }
  Z
}
