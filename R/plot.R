# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Feb. 2018
# Version 2.3
# Licence GPL v3


if (!isGeneric("plot")) {
  setGeneric("plot", function(x,y,...)
    standardGeneric("plot"))
}  


setMethod("plot", signature(x='sdmEvaluate'),
          function(x,y,smooth=TRUE,...) {
            if (missing(y)) y <- 'roc'
            else {
              ar <- c('roc','sensitivity','specificity','TSS','Kappa','NMI','phi','ppv','npv','ccr','mcr','or','ommission','commission','predicted.prevalence')
              y <- .pmatch(y,ar)[1]
              if (is.na(y)) stop('when x is a sdmEvaluate object, y should be a statistic to plot, the name in the y is not recognised!')
            }
            
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            
            if (y == 'roc') {
              r <- .roc(x@observed,x@predicted)
              .plot.roc(r,auc=x@statistics$AUC,smooth=smooth,...)
            } else {
              r <- .getEvalThresholds(x@observed,x@predicted,y)
              mx <- r[which.max(r[,2]),1]
              if (smooth) {
                rr <- try(supsmu(r[,1],r[,2],bass=0),silent=TRUE)
                if (!inherits(rr,'try-error')) {
                  r <- rr
                  rm(rr)
                }
              }
              dot <- list(...)
              ndot <- names(dot)
              dot[['x']] <- r
              
              if (!'xlab' %in% ndot) dot[['xlab']] <- 'Thresholds'
              if (!'ylab' %in% ndot) dot[['ylab']] <- y
              if (!'main' %in% ndot) dot[['main']] <-  paste('max at:',round(mx,2))
              if (!'type' %in% ndot) dot[['type']] <- 'l'
              do.call(plot,dot)
            }
          }
)
#----------------

setMethod("plot", signature(x='sdmdata'),
          function(x,y,sp=NULL,test=FALSE,col=c('blue','red'),xlab,ylab,main,xlim,ylim,pch,...) {
            xy <- coordinates(x)
            
            if (missing(y)) {
              if (!is.null(xy)) y <- 'map'
              else y <- NULL
            } else y <- tolower(y)
            if (missing(test)) test <- FALSE
            
            
            if (missing(sp)) {
              sp <- x@species.names
            } else {
              if (numeric(sp)) sp <- x@species.names[sp]
            }
            np <- length(sp)
            if (y == 'map') {
              if (is.null(xy)) stop('the species train data does not contain coordinates, chnage y for the other plots')
              else xy <- data.frame(x@info@coords)
              
              if (missing(xlab)) xlab <- 'X'
              if (missing(ylab)) ylab <- 'Y'
              if (missing(pch)) pch <- c(16,16)
              if (length(col) == 1) col <- rep(col,2)
              if (length(pch) == 1) pch <- rep(pch,2)
              rx <- range(xy[,2])
              xadd <- (rx[2] - rx[1]) * 0.01
              ry <- range(xy[,3])
              yadd <- (ry[2] - ry[1])*0.01
              if (missing(xlim)) xlim <- rx + (xadd*c(-1,1))
              if (missing(ylim)) ylim <- ry + (yadd*c(-1,1))
              if (missing(main)) main <- 'presence/absence map'
              
              if (np > 16) {
                warning('Due to larger number of plots, only the first 16 species are used...')
                sp <- sp[1:16]
                par(mfrow=c(4,4))
              } else {
                w <- floor(sqrt(np))
                h <- ceiling(np/w)
                if (abs(w-h) > 1) {
                  w <- w+1
                  h <- h-1
                }
                par(mfrow=c(w,h))  
              }
              
              for (s in sp) {
                df <- as.data.frame(x,sp=s,grp=if(test) 'test' else NULL)
                w <- xy$rID %in% df$rID
                if (x@species[[sp]]@type == 'Presence-Absence') {
                  w1 <- df[,s] == 1
                  plot(xy[w1 & w,2:3],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=col[1],main=paste(main,'_',s),pch=pch[1],...)
                  w1 <- !w1
                  points(xy[w1 & w,2:3],col=col[2],pch=pch[2],...)
                } else {
                  plot(xy[w,2:3],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=col[1],main=paste(main,'_',s),pch=pch[1],...)
                }
              }
              
            }
            
          }
)


setMethod("plot", signature(x='.sdmCalibration'),
          function(x,y,...) {
            if (missing(y)) y <- NULL
            dot <- list(...)
            ndot <- names(dot)
            
            if (!'xlab' %in% ndot) dot[['xlab']] <- "Predicted Probability of Occurrence"
            if (!'ylab' %in% ndot) dot[['ylab']] <- "Proportion of Observed Occurrences"
            if (!'xlim' %in% ndot) dot[['xlim']] <-  c(0,1)
            if (!'ylim' %in% ndot) dot[['ylim']] <-  c(0,1)
            if (!'main' %in% ndot) dot[['main']] <- 'Calibration Plot'
            if (!'sub' %in% ndot) dot[['sub']] <- paste('statistic = ',round(x@statistic,3),sep='')
            if ('cex' %in% ndot) {
              cex <- dot[['cex']]
              w <- which(ndot == 'cex')
              dot <- dot[-w]
              ndot <- ndot[-w]
            } else cex <- 2
            
            if ('pch' %in% ndot) {
              pch <- dot[['pch']]
              w <- which(ndot == 'pch')
              dot <- dot[-w]
              ndot <- ndot[-w]
            } else pch <- 16
            
            dot[['x']] <- 2
            dot[['y']] <- 1
            do.call(plot,dot)
            abline(a=0,b=1,lty=2)
            points(x@calibration[,1:2],pch=16,cex=2)
            #legend(0.7,0.05,legend = paste('Calibration = ',round(o@statistic,3),sep=''),text.width=0.3)
          }
)
#-------


setMethod("plot", signature(x='.varImportance'),
          function(x,y,...) {
            if (missing(y)) y <- 'corTest'
            else {
              y <- y[1]
              if (is.character(y)) y <- .pmatch(y,c('corTest','AUCtest'))
              else if (is.numeric(y)) {
                if (!y %in% c(1,2)) {
                  y <- 'corTest'
                  warning('y should be 1 or 2, it is changed to 1 (i.e., corTest)')
                } else y <- c('corTest','AUCtest')[y]
              } else {
                y <- 'corTest'
                warning('y is not identified... default is used (i.e., corTest)')
              }
            }
            
            dot <- list(...)
            ndot <- names(dot)
            if (!'xlab' %in% ndot) dot[['xlab']] <- "Relative Variable Importance"
            if (!'horiz' %in% ndot) dot[['horiz']] <- TRUE
            if (!'names.arg' %in% ndot) dot[['names.arg']] <- x@variables
            if (!'col' %in% ndot) dot[['col']] <-'#DDE9EB'
            if (!'cex.names' %in% ndot) dot[['cex.names']] <- 0.8
            if (!'las' %in% ndot) dot[['las']] <- 1
            dot[['height']] <- x@varImportance[,y]
            
            do.call(barplot,dot)
          }
)

#-------
if (!isGeneric("boxplot")) {
  setGeneric("boxplot", function(x, ...)
    standardGeneric("boxplot"))
}	


setMethod('boxplot', signature(x='sdmEvaluate'), 
          function(x,notch = FALSE,col='#DDE9EB',names,...) {
            
            p <- x@predicted
            x <- x@observed
            
            w <- which(!is.na(p))
            p <- p[w]; x <- x[w]
            w <- which(!is.na(x))
            p <- p[w]; x <- x[w]
            
            w1 <- which(x == 1)
            w2 <- which(x == 0)
            
            if (length(w1) > 0) {
              if (length(w2) > 0) {
                if (missing(names)) names <- c('Absence','Presence')
                boxplot(p[w2],p[w1],notch=notch,names=names,col=col,...)
              } else {
                if (missing(names)) names <-'Absence'
                boxplot(p[w2],notch=notch,col=col,...)
              }
            } else {
              if (length(w2) > 0) {
                if (missing(names)) names <- 'Presence'
                boxplot(p[w1],notch=notch,names=names,col=col,...)
              }
            }
          }
)