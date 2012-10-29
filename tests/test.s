# attach the permax library before sourcing this file
#library(permax,lib.loc='/usr/stats/gray/chips')

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
print(summary(uu))
#plots
ps.options(image.colors=cbind(seq(0,1,length=30),seq(0,1,length=30),seq(1,0,length=30)))
postscript(file='figs.ps')
#motif("-xrm 'sgraphMotif.colorSchemes : background : black; lines : yellow cyan magenta green MediumBlue red; text : white yellow cyan magenta green MediumBlue red; images : blue 30 yellow'")
#motif("-xrm 'sgraphMotif.defaultColorScheme : 3'") # cyan to red
#graphsheet(file='f:\\chips\\fig1.jpg',format='jpg',num.image.colors=2,
#  num.image.shades=30,image.color.table="0,0,255|255,255,0")

plot(uu,exp1,ig1=1:5,cex=.7)
plot.expr(exp1[1:20,])
dev.off()

u4 <- permax(exp1,1:5, nperm=1000)
print(summary(u4))
print(attr(u4,'seed.start'))
print(attr(u4,'seed.end'))

u5 <- permsep(exp1,ig1=1:5,nperm=1000,WHseed=attr(u4,'seed.end'))
print(attr(u5,'seed.start'))
print(attr(u5,'seed.end'))

ucor <- permcor(exp1[,3:8],1:6,nperm=0)
print(summary(ucor))
ucor <- permcor(exp1,1:15)
print(summary(ucor))
print(attr(ucor,'seed.end'))
#ucor <- permcor(exp1[1:50,1:5],1:5,nperm=500)
#for (i in 1:5) print(cor(log(exp1[i,1:5]),1:5))

# should give same p-values:
u1 <- permax(exp1[,c(1:3,5:8)],1:3)
u2 <- permcor(exp1[,c(1:3,5:8)],c(1,1,1,0,0,0,0),nperm=0)
print(range(u1$p.lower-u2$p.lower))
print(range(u1$p.upper-u2$p.upper))

uuu <- permsep(exp1,1:5)
print(uuu$dtcs)
uuu <- permsep(exp1,1:5,nperm=1000)
