# By R Gray, DFCI
# Copyright (C) 2000, 2002 Robert Gray
# Distributed under the GNU Public License (see the file COPYING)

permax <- function(data,ig1,nperm=0,logs=TRUE, ranks=FALSE, min.np=1,
         ig2=NULL, WHseed=NULL, cluster=NULL, stratify=FALSE, weights=NULL,
         nl=50,nr=50,permute.cluster=FALSE,signed.rank=FALSE) {
#
  cl <- match.call()
  data <- as.matrix(data)
  if (logs) {
	tmp <- data<=0
	if(any(tmp))
	  data[tmp] <- 1
	data <- log(data)
   }
  if (!is.null(ig2)) {
# remove unused columns, & adjust ig1
    i2 <- ig1i <- rep(FALSE,ncol(data))
    ig1i[ig1] <- TRUE
    i2[c(ig1,ig2)] <- TRUE
    data <- data[,i2]
    if (!is.null(cluster)) cluster <- cluster[i2]
    ig1 <- (1:ncol(data))[ig1i[i2]]
  }
  dmin <- min(data)
  n1 <- length(ig1)
  n2 <- ncol(data)-n1
# compute summary statistics
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
  data <- data[sub,] 
  n <- nrow(data)
  if (!is.null(cluster)) {
    trt <- rep(2,n1+n2)
    trt[ig1] <- 1
    mclust <- table(cluster)
    nclust <- length(mclust)
    mct1 <- table(cluster,trt)[,1]
    if (permute.cluster) {
      if (stratify) stop('permute.cluster and stratify cannot both == TRUE')
      if (any(mclust != mct1 & mct1 != 0)) 
        stop('clusters cannot contain both groups when permute.cluster==TRUE')
      ipc <- 1
      o <- order(trt,cluster)
    } else {
      ipc <- 0
      o <- order(cluster,trt)
    }
    data <- data[,o]
    trt <- trt[o]
    ig1 <- (1:ncol(data))[trt==1]
    if (stratify) {
      istrt <- 1 
      if (is.null(weights)) weights <- rep(1,nclust)/nclust
      else {
        if (length(weights) != nclust) 
          stop(paste('weights must have length', format(nclust)))
      }
    } else {
      istrt <- 0
    }
    if (signed.rank) {
      if (max(mclust) !=2 | min(mclust) != 2) 
        stop('signed.rank requires paired data')
      if (stratify) stop('stratify and signed.rank cannot both = TRUE')
      d1 <- data[,ig1]-data[,-ig1]
      d2 <- t(apply(abs(d1),1,rank))
      d2[d1==0] <- 0 # ties don't contribute
      d2 <- ifelse(d1<0,-d2,d2)
      data[,ig1] <- d2 
      data[,-ig1] <- -d2
      irnk <- 2
    }
  } else {
    data <- cbind(data[,ig1],data[,-ig1])
    ig1 <- 1:n1
    nclust <- 1
    mclust <- n1+n2
    mct1 <- n1
    istrt <- 0
    ipc <- 0
  }
  if (!signed.rank) {
    if (ranks) {
      if (istrt == 1) {
        t1 <- c(0,cumsum(mclust))
        for (i in 1:nclust) {
          ii <- (t1[i]+1):t1[i+1]
          data[,ii] <- t(apply(data[,ii,drop=FALSE],1,rank))
        }
      } else {
        data <- t(apply(data,1,rank))
      }
      irnk <- 1
    } else {
      irnk <- 0
    }
  }
  if (nperm>0) {
    if (is.null(WHseed)) WHseed <- floor(30000*runif(3))+1
  } else {
    WHseed <- c(0,0,0)
    if (ipc==1) {
      nct1 <- sum(as.numeric(mct1>0))
      nn <- exp(sum(log(2:nclust))-sum(log(2:nct1))-sum(log(2:(nclust-nct1))))
    } else {
      nn <- 0
      for (i in 1:nclust) if (mclust[i]>mct1[i] & mct1[i]>0) nn <- nn +
        sum(log(1:mclust[i]))-sum(log(1:mct1[i]))-
        sum(log(1:(mclust[i]-mct1[i])))
      nn <- exp(nn)
    }
    cat('statistics will be computed for all',format(nn),'combinations\n')
  }
  Z <- .Fortran('ptnstd',d=as.single(data),n=as.integer(n),ng=as.integer(
    ncol(data)),ng1=as.integer(n1),stat=single(n),
    as.integer(nclust),as.integer(mclust),as.integer(mct1),
    as.integer(ig1),as.integer(irnk),as.integer(istrt),
    weights=as.single(weights),ipc=as.integer(ipc))[c('d','stat','weights')]
  Z2 <- Z$stat
  weights <- Z$weights
  data <- matrix(Z$d,nrow=n,dimnames=dimnames(data))
  o <- order(Z2)
  Z2 <- Z2[o]
  if (nl>n) nl <- round(n/2)
  if (nr>n) nr <- round(n/2)
  crit <- c(Z2[nl],Z2[n-nr+1])
  data <- data[o,]
  Z <- .Fortran('ptn',d=as.single(data),n=as.integer(n),ng=as.integer(
    ncol(data)),ng1=as.integer(n1),stat=as.single(Z2),pind.lower=integer(n),
    pind.upper=integer(n),p.lower=integer(n),p.upper=integer(n),
    nperm=as.integer(nperm),ix=as.integer(WHseed),as.integer(nclust),
    as.integer(mclust),as.integer(mct1),as.integer(c(ig1,rep(0,n2))),
    integer(ncol(data)+n2),as.integer(irnk),as.integer(istrt),
    as.single(weights),nlr=as.integer(c(nl,nr)),as.single(crit),dist=single(6),
    iflag=integer(nclust),ipc=as.integer(ipc),single(n))[c(6:11,22)]
  if (nperm>0) endseed <- Z$ix
  if (!ranks & !stratify & !signed.rank) {
    Z2 <- .Fortran('tst2',as.single(cbind(data[,ig1],data[,-ig1])),
             as.integer(n1),as.integer(n2),as.integer(n),single(n))[[5]]
  }
  dist <- c(nl,Z$dist[1:3],nr,Z$dist[4:6])
  names(dist) <- c('nl','prop.nl','prop.1l','ave.l','nr','prop.nr','prop.1r','ave.r')
  Z <- data.frame(stat=Z2,pind.lower=Z$pind.lower/Z$nperm,
      pind.upper=Z$pind.upper/Z$nperm,p.lower=Z$p.lower/Z$nperm,
      p.upper=Z$p.upper/Z$nperm)
  m1 <- m1[sub]
  m2 <- m2[sub]
  if (logs){
    d1 <- data.frame(m1=m1,m2=m2,s1=s1[sub],s2=s2[sub],np1=npos1[sub],
      np2=npos2[sub],mdiff=exp(m1)-exp(m2),mrat=exp(m1-m2))
  } else {
    d1 <- data.frame(m1=m1,m2=m2,s1=s1[sub],s2=s2[sub],np1=npos1[sub],
      np2=npos2[sub],mdiff=m1-m2,mrat=m1/m2)
  }
  Z <- cbind(Z,d1[o,])
  row.names(Z) <- dimnames(data)[[1]]
  class (Z) <- c('permax','data.frame')
  attr(Z,'dist') <- dist
  attr(Z,'call') <- cl
  if (nperm>0) {
    attr(Z,'seed.start') <- WHseed
    attr(Z,'seed.end') <- endseed
  }
  Z
}

summary.permax <- function(object, data, nl=25, nr=25, ...) {
    ans<-list(nl=nl, nr=nr, objectl=NULL, objectr=NULL, dataL=NULL, 
              dataR=NULL, lower=NULL, upper=NULL, dist=attr(object,'dist'), 
              call=attr(object,'call'))
    o <- order(object$stat)
    if (nl>0) {
      ans$objectl <- object[o[1:nl],]
      if (!missing(data)) {
        data <- as.matrix(data)
        ans$dataL <- data[match(row.names(ans$objectl), dimnames(data)[[1]],0),]
      }
    }
    if (nr>0) {
      ans$objectr <- object[rev(o)[1:nr],]
      if (!missing(data)) {
        data <- as.matrix(data)
        ans$dataR <- data[match(row.names(ans$objectr),dimnames(data)[[1]],0),]
      }
    }
    ans$lower <- row.names(ans$objectl)
    ans$upper <- row.names(ans$objectr)
    class(ans) <- "summary.permax"
    ans
}

print.summary.permax <- function(x, digits = max(3,
                                    getOption("digits") - 3), ...)
{
  cat('Call:\n')
  print(x$call)
  if (!is.null(x$dist)) {
    cat('Summary of the Null Permutation Distribution of # Positives:\n')
    print(x$dist)
  }
  print(x$objectl[,-c(2,3,5)], digits=digits)
  if( !is.null(x$dataL) ) print(x$dataL, digits=digits)
  print(x$objectr[,-c(2,3,4)], digits=digits)
  if( !is.null(x$dataR) ) print(x$dataR, digits=digits)
}

permsep <- function(data,ig1,nperm=0,ig2=NULL,WHseed=NULL, cluster=NULL, stratify=FALSE, permute.cluster=FALSE) {
#
  cl <- match.call()
  data <- as.matrix(data)
  if (!is.null(ig2)) {
# remove unused columns, & adjust ig1
    i2 <- ig1i <- rep(FALSE,ncol(data))
    ig1i[ig1] <- TRUE
    i2[c(ig1,ig2)] <- TRUE
    data <- data[,i2]
    if (!is.null(cluster)) cluster <- cluster[i2]
    ig1 <- (1:ncol(data))[ig1i[i2]]
  }
  n <- nrow(data)
  n1 <- length(ig1)
  n2 <- ncol(data)-n1
  if (!is.null(cluster)) {
    trt <- rep(2,n1+n2)
    trt[ig1] <- 1
    mclust <- table(cluster)
    nclust <- length(mclust)
    mct1 <- table(cluster,trt)[,1]
    if (permute.cluster) {
      if (stratify) stop('permute.cluster and stratify cannot both == TRUE')
      if (any(mclust != mct1 & mct1 != 0)) 
        stop('clusters cannot contain both groups when permute.cluster==TRUE')
      ipc <- 1
      o <- order(trt,cluster)
    } else {
      ipc <- 0
      o <- order(cluster,trt)
    }
    data <- data[,o]
    trt <- trt[o]
    cluster <- cluster[o]
    ig1 <- (1:ncol(data))[trt==1]
    ig2 <- (1:ncol(data))[trt != 1]
    if (stratify) {
      istrt <- 1 
    } else {
      istrt <- 0
    }
  } else {
    data <- cbind(data[,ig1],data[,-ig1])
    ig1 <- 1:n1
    ig2 <- n1+(1:n2)
    nclust <- 1
    mclust <- n1+n2
    mct1 <- n1
    istrt <- 0
    ipc <- 0
  }
  if (nperm>0) {
    if (is.null(WHseed)) WHseed <- floor(30000*runif(3))+1
  } else {
    WHseed <- c(0,0,0)
    if (ipc==1) {
      nct1 <- sum(as.numeric(mct1>0))
      nn <- exp(sum(log(1:nclust))-sum(log(1:nct1))-sum(log(1:(nclust-nct1))))
    } else {
      nn <- 0
      for (i in 1:nclust) if (mclust[i]>mct1[i] & mct1[i]>0) nn <- nn +
        sum(log(1:mclust[i]))-sum(log(1:mct1[i]))-
        sum(log(1:(mclust[i]-mct1[i])))
      nn <- exp(nn)
    }
    cat('statistics will be computed for all',format(nn),'combinations\n')
  }
  Z <- .Fortran('ptc',d=as.single(data),n=as.integer(n),ng=as.integer(
    ncol(data)),ng1=as.integer(length(ig1)),ics=integer(nrow(data)),
    nperm=as.integer(nperm),dtcs=integer(4),ix=as.integer(WHseed),
    as.integer(nclust),as.integer(mclust),as.integer(mct1),as.integer(ig1),
    as.integer(c(ig2,rep(0,n1+n2))),as.integer(istrt),iflag=integer(nclust),
    ipc=as.integer(ipc),igc=integer(nclust))[5:8]
  if (nperm>0) endseed <- Z$ix
  Z$dtcs <- c(Z$dtcs[1],Z$dtcs[2:4]/Z$nperm)
  cat('# attributes with complete separation:',Z$dtcs[1],'\n')
  cat('proportion of permutations with as many or more (p-value):',
      Z$dtcs[2],'\n')
  cat('average # per permutation:',Z$dtcs[3],'\n')
  cat('proportion of permutations with any:',Z$dtcs[4],'\n')
  Z <- Z[-c(2,4)]
  names(Z$dtcs) <- c('obs.num','pval.no','avenum.per.perm','prop.with.any')
  names(Z$ics) <- dimnames(data)[[1]]
  if (nperm>0) {
    attr(Z,'seed.start') <- WHseed
    attr(Z,'seed.end') <- endseed
  }
  attr(Z,'call') <- cl
  Z
}

rowperm <- function(x) {
# applies a separate permutation
# to the elements in each row of the input array x
  x <- as.matrix(x)
  dimx <- dim(x)
  dimnmx <- dimnames(x)
  x <- .Fortran('pa',as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),
     as.double(runif(nrow(x)*(ncol(x)-1))))[[1]]
  dim(x) <- dimx
  dimnames(x) <- dimnmx
  x
}

plot.permax <- function(x, data, nl=25, nr=25, logs=TRUE, ig1=NULL,
                        ig2=NULL, clmn.lab=dimnames(data)[[2]], 
                        row.lab=dimnames(data)[[1]], clmn.off=NULL, 
                        row.off=NULL, ...) {
  o <- order(x$stat)
  if (nl>0) xl <- row.names(x)[o[1:nl]] else xl <- NULL
  if (nr>0) xr <- rev(row.names(x)[rev(o)[1:nr]]) else xr <- NULL
  data <- data[c(xl,xr),]
  plot.expr(data,logs=logs,ig1,ig2,clmn.lab,row.lab,
            clmn.off,row.off,...)
}

plot.expr <- function(x, logs=TRUE, ig1=NULL, ig2=NULL, 
                      clmn.lab=dimnames(x)[[2]], row.lab=dimnames(x)[[1]],
                      clmn.off=NULL, row.off=NULL, ...) {
# All rows of data are plotted in the order given, so appropriate sorting
#   and subsetting should be done prior to the call.
  x <- as.matrix(x)
  if (logs) x <- log(x)
  m1 <- apply(x,1,mean)
  x <- x-m1
  s1 <- sqrt(apply(x^2,1,sum)/(ncol(x)-1))
  x <- x/s1
# with unequal group sizes, standardization as above tends to guarantee
# that the smaller group takes more extreme values than the larger.
# the following shifts the overall mean so the mean of each group is
# equal in magnitude and opposite in sign.
  if (!is.null(ig1)) {
    if (is.null(ig2)) ig2 <- (1:ncol(x))[-ig1]
    m1 <- apply(x[,ig1],1,mean)
    m2 <- apply(x[,ig2],1,mean)
    x <- x-(m1+m2)/2
  }
  image(1:ncol(x),1:nrow(x),t(x),xaxt='n',yaxt='n', xlab="",
        ylab="",...)
  if (is.null(row.off)) row.off <- -(ncol(x)/35)
  if (is.null(clmn.off)) clmn.off <- -(nrow(x)/12)
  text(row.off, 1:nrow(x), row.lab, xpd = TRUE)
  text(1:ncol(x), clmn.off, clmn.lab, srt=270, xpd=TRUE)
  invisible()
}

permcor <- function(data, phen, nperm=1000, logs=TRUE, ranks=FALSE,
        min.np=1, WHseed=NULL, cluster=NULL, stratify=FALSE, weights=NULL,
        permute.cluster=FALSE) {
#
  cl <- match.call()
  data <- as.matrix(data)
  if (logs) data <- log(data)
  dmin <- min(data)
  ng <- ncol(data)
  d1 <- data-dmin
  d1[d1>0] <- 1
  npos <- d1 %*% rep(1,ng)
  data <- data[npos>=min.np,]
  n <- nrow(data)
  if (!is.null(cluster)) {
    o <- order(cluster)
    data <- data[,o]
    phen <- phen[o]
    mclust <- table(cluster)
    nclust <- length(mclust)
    if (permute.cluster) {
      if (stratify) stop('permute.cluster and stratify cannot both == TRUE')
      uv <- tapply(phen,cluster,min)
      if (max(abs(tapply(phen,cluster,max)-uv))>0)
        stop('clusters cannot contain both groups when permute.cluster==TRUE')
      ipc <- 1
    } else {
      ipc <- 0
      uv <- 0
    }
    if (stratify) {
      istrt <- 1 
      if (is.null(weights)) weights <- rep(1/nclust,nclust)
      else {
        if (length(weights) != nclust) 
          stop(paste('weights must have length', format(nclust)))
        weights <- weights/sum(weights)
      }
    } else {
      istrt <- 0
    }
    if (stratify) istrt <- 1 else istrt <- 0
  } else {
    nclust <- 1
    mclust <- ng
    istrt <- 0
    ipc <- 0
    uv <- 0
  }
  if (ranks) {
    if (istrt == 1) {
      t1 <- c(0,cumsum(mclust))
      for (i in 1:nclust) {
        ii <- (t1[i]+1):t1[i+1]
        data[,ii] <- t(apply(data[,ii,drop=FALSE],1,rank))
      }
    } else {
      data <- t(apply(data,1,rank))
    }
  }
  if (nperm>0) {
    if (is.null(WHseed)) WHseed <- floor(30000*runif(3))+1
  } else {
    WHseed <- c(0,0,0)
    if (ipc==1) {
      nn <- exp(sum(log(1:nclust)))
    } else {
      nn <- 0
      for (i in 1:nclust) nn <- nn + sum(log(1:mclust[i]))
      nn <- exp(nn)
    }
    cat('statistics will be computed for all',format(nn),'permutations\n')
  }
  Z <- .Fortran('corstd',d=as.single(data),n=as.integer(n),ng=as.integer(
    ncol(data)),phen=as.single(phen),stat=single(n),
    as.integer(nclust),as.integer(mclust),as.integer(istrt),
    weights=as.single(weights))[c('d','phen','stat')]
  Z2 <- Z$stat
  data <- matrix(Z$d,nrow=n,dimnames=dimnames(data))
  phen <- Z$phen
  o <- order(Z2)
  Z2 <- Z2[o]
  data <- data[o,]
  Z <- .Fortran('ptcor',d=as.single(data),n=as.integer(n),ng=as.integer(ng),
    phen=as.single(phen),stat=as.single(Z2),pind.lower=integer(n),
    pind.upper=integer(n),p.lower=integer(n),p.upper=integer(n),
    nperm=as.integer(nperm),ix=as.integer(WHseed),
    ig=integer(2*ng),nclust=as.integer(nclust),
    mclust=as.integer(mclust),istrt=as.integer(istrt),as.single(weights),
    iflag=integer(nclust),ipc=as.integer(ipc),uv=as.single(uv),single(n))[6:11]
  if (nperm>0) endseed <- Z$ix
  Z <- data.frame(stat=Z2,pind.lower=Z$pind.lower/Z$nperm,
    pind.upper=Z$pind.upper/Z$nperm,p.lower=Z$p.lower/Z$nperm,
    p.upper=Z$p.upper/Z$nperm,np=npos[npos>=min.np][o])
  row.names(Z) <- dimnames(data)[[1]]
  class (Z) <- c('permcor','permax','data.frame')
  attr(Z,'call') <- cl
  if (nperm>0) {
    attr(Z,'seed.start') <- WHseed
    attr(Z,'seed.end') <- endseed
  }
  Z
}

impmv <- function(data,IFUN=mean,min.nonmiss=2,...) {
### imputes missing values in a row of d using the function FUN applied 
### to the non-missing values in the row
  imp <- function(x,IFUN,...) {
    i1 <- is.na(x)
    if (any(i1)) {
      nval <- sum(as.numeric(!i1))
      if (nval>=min.nonmiss) x[i1] <- IFUN(x[!i1],...)
    }
    x
  }
  dn <- dimnames(data)
  data <- t(apply(data,1,imp,IFUN=IFUN,...))
  dimnames(data) <- dn
  data
}
