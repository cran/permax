plot.permax <- function(x, data, nl=25, nr=25, logs=TRUE, ig1=NULL,
                        ig2=NULL, clmn.lab=dimnames(data)[[2]], 
                        row.lab=dimnames(data)[[1]], clmn.off=NULL, 
                        row.off=NULL, ...) {
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
# clmn.lab= labels for the columns in the array 
# row.lab= labels for the rows in the array 
# clmn.off= offset for printing the column labels (<0 to put labels
#   outside the plot) 
# row.off= offset for printing the row labels (<0 to put labels outside
#   the plot) 
# A graphics device with support for appropriate image colors must be
#   specified prior to calling this function
  o <- order(x$stat)
  xl <- row.names(x)[o[1:nl]]
  xr <- rev(row.names(x)[rev(o)[1:nr]])
  plot.expr(data[c(xl,xr),],logs=logs,ig1,ig2,clmn.lab,row.lab,
            clmn.off,row.off,...)
}

plot.expr <- function(x, logs=TRUE, ig1=NULL, ig2=NULL, 
                      clmn.lab=dimnames(x)[[2]], row.lab=dimnames(x)[[1]],
                      clmn.off=NULL, row.off=NULL, ...) {
# plots the expression levels for the rows of the arrary data
# x = expression level array
# logs = if true, function takes logs of values in data
# if ig1 is given a non-null value (which must be a vector of integers
#   corresponding to columns in data, then the rows of data will be
#   standardized so the mean of the columns in ig1 and the mean of the
#   columns in ig2 are equal in magnitude and opposite in sign (if ig2 is
#   not specified, it defaults to include all the columns not in ig1)
#   if ig1 is NULL, then the rows of data are standardized to have mean 0
#   In either case, the rows are also standardized to have variance 1
# clmn.lab= labels for the columns in the array 
# row.lab= labels for the rows in the array 
# clmn.off= offset for printing the column labels (<0 to put labels
#   outside the plot) 
# row.off= offset for printing the row labels (<0 to put labels outside
#   the plot) 
# A graphics device with support for appropriate image colors must be
#   specified prior to calling this function
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
  if (is.null(clmn.off)) clmn.off <- -(nrow(x)/9)
  text(row.off, 1:nrow(x), row.lab, xpd = TRUE)
  text(1:ncol(x), clmn.off, clmn.lab, srt=270, xpd=TRUE)
#  positioning is different in graphsheet() plots
#  text(1:ncol(x),-(nrow(x)/7),dimnames(x)[[2]],srt=270,xpd=TRUE,...)
  invisible()
}
