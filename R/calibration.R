# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July. 2016
# Version 1.1
# Licence GPL v3


.cutx <- function(x,n) {
  x <- round(x,3)
  brk <- (max(x) - min(x)) / n
  r <- rep(NA,length(x))
  l <- min(x)
  w <- which(x >= l & x <= (l+brk))
  r[w] <- 1
  for (i in 2:n) {
    b <- l + (i-1)*brk
    if(i < n) w <- which(x > b & x <= round((b+brk),3))
    else w <- which(x > b & x <= max(x))
    if ( length(w) > 0) r[w] <- i
  }
  r
}

#------
.bin <- function(x,th=0.99) {
  Loop <- TRUE
  i <-  2
  while (Loop) {
    y <- .cutx(x,i)
    if (cor(x,y,method='spearman') >= th) Loop <- FALSE
    i <- i + 1
  }
  return(i-1)
}
#--------------------
if (!isGeneric("calibration")) {
  setGeneric("calibration", function(x,p,nbin,weight,...)
    standardGeneric("calibration"))
}


setMethod("calibration", signature(x='vector',p='vector'),
          function(x,p,nbin,weight,...) {
            if (missing(weight)) weight <- TRUE
            w <- which(!is.na(p))
            p <- p[w]; x <- x[w]
            w <- which(!is.na(x))
            p <- p[w]; x <- x[w]
            
            if (missing(nbin)) nbin <- 10
            else if (nbin == 'seek') {
              nbin <- .bin(p)
            } else {
              if (!is.numeric(nbin)) {
                nbin <- 10
                warning('nbin is not numeric! the default value (i.e., 10) is used...')
              }
            }
            
            d <- data.frame(matrix(NA,nrow=nbin,ncol=3))
            colnames(d) <- c('pedicted_pobability_bin','observed_poportion','Weight')
            
            r <- max(p) - min(p)
            brk <- r / nbin
            for (i in 1:nbin) {
              b <- min(p) + (i - 1) * brk 
              if(i < nbin) w <- which(p >= b & p < (b+brk))
              else w <- which(p >= b & p <= (b+brk))
              
              if (length(w) > 0) {
                d[i,2] <-  length(which(x[w] == 1)) / length(x[w])
                d[i,1] <- b+(brk/2)
                d[i,3] <- length(w) / length(x)
              }
            }
            
            d <- d[which(!is.na(d[,1])),]
            
            if (!weight) d[,3] <- 1 / nrow(d)
            
            ds <- d[,1]
            ds <- ifelse(ds <= 0.5,ds + 0.5,ds - 0.5)
            
            st.f <- sqrt(sum(((d[,2] - ds)^2) * d[,3]))
            
            o <- new('.sdmCalibration')
            
            o@statistic <- 1-(sqrt(sum(((d[,1] - d[,2])^2) * d[,3]))/st.f)
            o@calibration <- d
            
            return(o)
          }
)

setMethod("calibration", signature(x='sdmEvaluate','missing'),
          function(x,p,nbin,weight,...) {
            
            if (missing(weight)) weight <- TRUE
            if (missing(nbin)) nbin <- 10
            
            calibration(x@observed,x@predicted,nbin=nbin,weight=weight,...)
          }
)

