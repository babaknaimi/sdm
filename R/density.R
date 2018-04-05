# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July. 2016
# Version 1.1
# Licence GPL v3


if (!isGeneric("density")) {
  setGeneric("density", function(x, ...)
    standardGeneric("density"))
}


setMethod('density', signature(x='sdmEvaluate'), 
          function(x,xlab,xlim,ylim,col,lwd,lty,main,...) {
            if (missing(xlab)) xlab <- 'Predicted Probabilities'
            if (missing(ylim)) ylim <- NULL
            if (missing(xlim)) xlim <- c(0,1)
            if (missing(lwd)) lwd <- 2
            if (missing(col)) col <- c('darkblue','red')
            else if (length(col) == 1) col <- c(col,col)
            if (missing(main)) main <- ''
            if (missing(lty)) lty <- 1
            
            
            p=x@predicted
            x=x@observed
            w <- which(!is.na(p))
            p <- p[w]; x <- x[w]
            w <- which(!is.na(x))
            p <- p[w]; x <- x[w]
            
            w1 <- which(x == 1)
            w2 <- which(x == 0)
            if (length(w1) > 0) {
              d1 <- density(p[w1])
              if (length(w2) > 0) {
                d2 <- density(p[w2],bw=d1$bw)
                if (is.null(ylim)) ylim <- c(0,max(c(d1$y,d2$y)))
                plot(d1,col=col[1],xlim=xlim,lwd=lwd,main=main,xlab=xlab,ylim=ylim,...)
                lines(d2,col=col[2],lwd=lwd)
                legend('topleft',legend=c('Presence','Absence'),col=col,lty=c(lty,lty),text.width = 0.2,lwd=lwd)
              } else {
                legend('topleft',legend=c('Presence'),col=col[1],lty=c(lty,lty),text.width = 0.2,lwd=lwd)
              }
            } else {
              if (length(w2) > 0) {
                plot(density(p[w2]),col=col[2],xlim=xlim,lwd=lwd,main=main,xlab=xlab,...)
                legend('topleft',legend=c('Absence'),col=col[2],lty=lty,text.width = 0.2,lwd=lwd)
              }
            }
          
          }
)
#---------
