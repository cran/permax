# attach the permax library before sourcing this file
#library(permax,lib.loc='/usr/stats/gray/chips')
#library(permax,lib.loc='~/R')
library(permax)
#dyn.load('../src/permax.so')
#source('../R/permax.R')
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
dimnames(exp1) <- list(paste('x',1:ngenes,sep=''),
                       paste('sample',1:ncol(exp1),sep=''))
exp1 <- round(exp1)
print(table(exp1<=20))
n1 <- apply(exp1>20,1,sum)
print(table(n1))
# permax
uu <- permax(exp1,1:5)
print(uus1 <- summary(uu))
# following should give same results
exp2 <- exp1[,c(1,2,7,10,3,8,11,6,5,9,12,15,13,14,4)]
exp3 <- cbind(ngenes:1,exp1,1:ngenes)
exp4 <- cbind(ngenes:1,exp2,1:ngenes)
uu2 <- permax(exp2,c(1,2,5,9,15))
uus2 <- summary(uu2)
uu3 <- permax(exp3,2:6,ig2=7:16)
uus3 <- summary(uu3)
uu4 <- permax(exp4,c(1,2,5,9,15)+1,ig2=c(4,5,7,8,9,11,12,13,14,15))
uus4 <- summary(uu4)
range(as.matrix(uus2$objectl)-as.matrix(uus1$objectl))
range(as.matrix(uus3$objectl)-as.matrix(uus1$objectl))
range(as.matrix(uus4$objectl)-as.matrix(uus1$objectl))
range(as.matrix(uus2$objectr)-as.matrix(uus1$objectr))
range(as.matrix(uus3$objectr)-as.matrix(uus1$objectr))
range(as.matrix(uus4$objectr)-as.matrix(uus1$objectr))

postscript(file='figs.ps')
par(cex=.7,mar=c(6,6,2,2))
plot(uu,exp1,ig1=1:5)
plot(uu,exp1,ig1=1:5,nr=40,nl=0)
plot.expr(exp1[1:20,])
dev.off()

u4 <- permax(exp1,1:5, nperm=1000)
print(summary(u4))
print(attr(u4,'seed.start'))
print(attr(u4,'seed.end'))

#clustered
clustind <- c(1,1,2,2,3,1,1,1,1,2,2,2,2,3,3)
uuc <- permax(exp1,1:5,nperm=5000,cluster=clustind,WHseed=attr(u4,'seed.end'))
summary(uuc,nl=5,nr=5) 
# stratified on cluster, using ranks; equal weights
uus <- permax(exp1,1:5,nperm=5000,cluster=clustind,stratify=TRUE,
          ranks=TRUE,WHseed=attr(uuc,'seed.end'))
summary(uus,nl=5,nr=5) 
print(attr(uus,'seed.end'))
# stratified on cluster, using ranks; average weighted by cluster sizes
uus <- permax(exp1,1:5,nperm=5000,cluster=clustind,stratify=TRUE,
           ranks=TRUE,weights=table(clustind),WHseed=attr(uus,'seed.end'))
summary(uus,nl=5,nr=5) 
print(attr(uus,'seed.start'))
print(attr(uus,'seed.end'))
summary(uus,nl=10,nr=0) 

clust2 <- c(1,1,2,2,2,3,3,4,4,5,5,5,6,6,6)
uu5 <- permax(exp1,1:5,nperm=5000,cluster=clust2,permute.cluster=TRUE,
           ranks=TRUE,WHseed=attr(uus,'seed.end'))
summary(uu5,nl=5,nr=5)
uu6 <- permax(exp1,1:5,nperm=0,cluster=clust2,permute.cluster=TRUE,ranks=TRUE)
print(uu6s <- summary(uu6,nl=5,nr=5))
apply(exp1[uu6s$upper,],1,rank)
uu6[uu6s$upper,]$pind.upper

ucor <- permcor(exp1[,3:8],1:6,nperm=0)
print(summary(ucor))
ucor <- permcor(exp1,1:15,WHseed=attr(uus,'seed.end'))
print(summary(ucor))
print(attr(ucor,'seed.end'))

ucorx <- permcor(exp1[,4:7],1:4,nperm=0)
summary(ucorx,nl=5,nr=5)
ucorx <- permcor(exp1[,4:7],1:4,nperm=0,cluster=c(1,2,1,2))
summary(ucorx,nl=5,nr=5)
ucorx <- permcor(exp1[,3:9],1:7,nperm=0,cluster=c(1,2,3,1,2,3,3))
summary(ucorx,nl=5,nr=5)
ucorx <- permcor(exp1[,3:9],1:7,nperm=10,cluster=c(1,2,3,1,2,3,3),WHseed=attr(uus,'seed.end'))
summary(ucorx,nl=5,nr=5)
ucorx <- permcor(exp1[,3:10],c(1,2,2,3,3,4,4,4),nperm=10,cluster=c(1,2,2,3,3,4,4,4),permute.cluster=TRUE,WHseed=attr(uus,'seed.end'))
summary(ucorx,nl=5,nr=5)
ucorx <- permcor(exp1[,3:10],c(1,2,2,3,3,4,4,4),nperm=0,cluster=c(1,2,2,3,3,4,4,4),permute.cluster=TRUE)
summary(ucorx,nl=5,nr=5)


postscript(file='figs2.ps')
par(cex=.7,mar=c(6,6,2,2))
plot(ucor,exp1) #columns of exp1 already sorted on Z
dev.off()

u8 <- permcor(exp1,1:15,cluster=clustind,WHseed=attr(ucor,'seed.end'))
summary(u8,nr=4,nl=4)
# correlations estimated within clusters; average weighted by cluster sizes
u8 <- permcor(exp1,1:15,cluster=clustind,stratify=TRUE,
              weights=table(clustind),WHseed=attr(u8,'seed.end'))
summary(u8,nr=4,nl=4)

# should give same p-values:
u1 <- permax(exp1[,c(1:3,5:8)],1:3)
u2 <- permcor(exp1[,c(1:3,5:8)],c(1,1,1,0,0,0,0),nperm=0)
print(range(u1$p.lower-u2$p.lower))
print(range(u1$p.upper-u2$p.upper))

uuu <- permsep(exp1,1:5)
print(uuu$dtcs)
uuu <- permsep(exp2,c(1,2,5,9,15))
uuu <- permsep(exp3,2:6,ig2=7:16)
uuu <- permsep(exp4,c(1,2,5,9,15)+1,ig2=c(4,5,7,8,9,11,12,13,14,15))
uuu <- permsep(exp1,1:5,nperm=1000)
print(attr(uuu,'seed.start'))
print(attr(uuu,'seed.end'))
u5 <- permsep(exp1,ig1=1:5,nperm=1000,WHseed=attr(uuu,'seed.end'))
print(attr(u5,'seed.start'))
print(attr(u5,'seed.end'))
uuu <- permsep(exp1,1:5,cluster=clustind,nperm=10000,WHseed=attr(u5,'seed.end'))
uuu <- permsep(exp1,1:5,nperm=2000,cluster=clustind,stratify=TRUE,
          WHseed=attr(uuu,'seed.end'))
uuu <- permsep(exp1,1:5,cluster=clust2,nperm=1000,WHseed=attr(uuu,'seed.end'),permute.cluster=TRUE)
uuu <- permsep(exp1,1:5,cluster=clust2,nperm=0,permute.cluster=TRUE)

x <- matrix(c(1:10,NA,NA),4,dimnames=list(format(1:4),c('a','b','c')))
impmv(x)
na.omit(impmv(x,min.nonmiss=3))
irn <- function(x,nc) rnorm(nc-length(x),mean(x),sqrt(var(x)))
impmv(x,IFUN=irn,nc=ncol(x)) 
impmv(t(x),IFUN=irn,nc=nrow(x)) 

x <- matrix(1:12,3)
rowperm(x)

# paired data -- checked the p-values and stats by hand
clust <- c(1,2,3,4,1,2,3,4)
dat <- rbind(c(1:8),c(8,1:7),c(5,1:7),c(8,7,1:6))
dimnames(dat) <- list(1:4,1:8)
permax(dat,1:4,cluster=clust,signed.rank=TRUE)
# note, order of the rows is flipped in the output of the following call
permax(dat,5:8,cluster=clust,signed.rank=TRUE)
