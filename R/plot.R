# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  March 2019
# Version 2.7
# Licence GPL v3
#---------------------

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

setMethod("plot", signature(x='.varImportanceList'),
          function(x,y,gg=TRUE,confidence=TRUE,...) {
            if (missing(gg)) gg <- TRUE
            
            if (missing(confidence)) confidence <- TRUE
            
            if (gg && !.require('ggplot2')) gg <- FALSE
            
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
            if (gg) {
              p2 <- "ggplot(drcc,aes(x=Value,y=Response)) + geom_bar(stat = 'identity',fill=col)+ scale_y_continuous(name = ylab,limits=ylim) + facet_grid(.~variable,scale='free') + geom_errorbar(aes(ymin=lower, ymax=upper),width=.3,position=position_dodge(.9))"
              
              #p2 <- .eval(p2,env=environment())
            } else {
              dot <- list(...)
              ndot <- names(dot)
              if (!'xlab' %in% ndot) dot[['xlab']] <- "Relative Variable Importance"
              if (!'horiz' %in% ndot) dot[['horiz']] <- TRUE
              if (!'names.arg' %in% ndot) dot[['names.arg']] <- x@variables
              if (!'col' %in% ndot) dot[['col']] <-'#DDE9EB'
              if (!'cex.names' %in% ndot) dot[['cex.names']] <- 0.8
              if (!'las' %in% ndot) dot[['las']] <- 1
              dot[['height']] <- x@varImportanceMean[[y]][,2]
              
              .bar <- do.call(barplot,dot)
              
              if (confidence) {
                segments(.bar, x@varImportanceMean[[y]][,3], .bar,x@varImportanceMean[[y]][,4], lwd = 1.5)
                
                arrows(.bar, x@varImportanceMean[[y]][,3], .bar,x@varImportanceMean[[y]][,4], lwd = 1.5, angle = 90,code = 3, length = 0.05)
              }
              
            }
            
          }
)

#-------
setMethod("plot", signature(x='.responseCurve'),
          function(x,y,gg=TRUE,mean=TRUE,confidence=TRUE,xlab,ylab,ylim,lty,lwd,col,cex.axis,cex.lab,main,...) {
            if (gg && !.require('ggplot2')) gg <- FALSE
            
            if (missing(mean)) mean <- TRUE
            
            if (missing(confidence)) confidence <- TRUE
            
            if (missing(y) || is.null(y)) n <- x@variables
            else {
              n <- y
              n <- n[n %in% x@variables]
              if (length(n) == 0) stop('the specified variable(s) in n does not exist in the responseCurve object!')
            }
            
            if (!is.null(x@categorical)) {
              nF <- x@categorical
              nF <- nF[nF %in% n]
              if (length(nF) == 0) nF <- NULL
              n <- .excludeVector(n,nF)
            } else nF <- NULL
            #-------
            if (missing(xlab)) xlab <- 'Variables'
            if (missing(ylab)) ylab <- 'Response'
            if (missing(lty)) lty <- 1
            if (missing(col)) col <- '#00395F'
            if (missing(lwd)) lwd <- 1
            if (missing(cex.axis)) cex.axis <- 1
            if (missing(cex.lab)) cex.lab <- 1
            
            if (missing(main)) main <- 'Response Curve'
            if (missing(ylim)) ylim <- NULL
            
            #--------
            if (gg) {
              if (x@multi) {
                if (mean) {
                  if (confidence) {
                    drc <- data.frame(Value=0,Response=0,lower=0,upper=0,variable='a')[0,]
                    for (nn in n) {
                      .n <- length(x@response[[nn]][,1])
                      .m <- apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,mean,na.rm=TRUE)
                      .ci <- 1.96 * apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,sd,na.rm=TRUE) / sqrt(.n)
                      drc <- rbind(drc,data.frame(Value=x@response[[nn]][,1],Response=.m,lower=.m - .ci,upper=.m + .ci,variable=nn))
                    }
                    
                    p1 <- "ggplot(drc,aes(x=Value,y=Response)) + geom_line(colour=col,size=lwd,linetype=lty) + geom_ribbon(aes(ymin=lower, ymax=upper), linetype=1, alpha=0.2) + 
                      facet_grid(.~variable,scales='free_x') + labs(y = ylab,x = xlab) + ggtitle(main) +
                      theme(axis.text=element_text(size=rel(cex.axis)),axis.title=element_text(size=rel(cex.lab)),plot.title = element_text(hjust = 0.5))"
                    p1 <- .eval(p1,env=environment())
                    if (!is.null(nF)) {
                      drcc <- data.frame(Value=0,Response=0,variable='a')[0,]
                      for (nn in nF) {
                        .n <- length(x@response[[nn]][,1])
                        .m <- apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,mean,na.rm=TRUE)
                        .ci <- 1.96 * apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,sd,na.rm=TRUE) / sqrt(.n)
                        drcc <- rbind(drcc,data.frame(Value=x@response[[nn]][,1],Response=.m,lower=.m - .ci,upper=.m + .ci,variable=nn))
                      }
                      
                      p2 <- "ggplot(drcc,aes(x=Value,y=Response)) + geom_bar(stat = 'identity',fill=col) + facet_grid(.~variable,scale='free') + 
                      geom_errorbar(aes(ymin=lower, ymax=upper),width=.3,position=position_dodge(.9)) +
                      labs(y = ylab,x = xlab)"
                      
                      p2 <- .eval(p2,env=environment())
                      if (!.require('gridExtra')) {
                        warning('you need the package gridExtra to make the plots printed in a single page!')
                        return(list(p1,p2))
                      } else {
                        return(.eval("grid.arrange(p1,p2)",env=environment()))
                      }
                    } else return(p1)
                    
                  } else {
                    drc <- data.frame(Value=0,Response=0,variable='a')[0,]
                    for (nn in n) {
                      drc <- rbind(drc,data.frame(Value=x@response[[nn]][,1],Response=apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,mean,na.rm=TRUE),variable=nn))
                    }
                    
                    p1 <- "ggplot(drc,aes(x=Value,y=Response)) + geom_line(colour=col,size=lwd,linetype=lty) + 
                      facet_grid(.~variable,scales='free_x') + scale_y_continuous(limits=ylim) + labs(y = ylab,x = xlab)  + ggtitle(main) +
                      theme(axis.text=element_text(size=rel(cex.axis)),axis.title=element_text(size=rel(cex.lab)),plot.title = element_text(hjust = 0.5))"
                    
                    p1 <- .eval(p1,env=environment())
                    
                    if (!is.null(nF)) {
                      drcc <- data.frame(Value=0,Response=0,variable='a')[0,]
                      for (nn in nF) {
                        .n <- length(x@response[[nn]][,1])
                        .m <- apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,mean,na.rm=TRUE)
                        .ci <- 1.96 * apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,sd,na.rm=TRUE) / sqrt(.n)
                        drcc <- rbind(drcc,data.frame(Value=x@response[[nn]][,1],Response=.m,lower=.m - .ci,upper=.m + .ci,variable=nn))
                      }
                      
                      p2 <- "ggplot(drcc,aes(x=Value,y=Response)) + geom_bar(stat = 'identity',fill=col)+ scale_y_continuous(name = ylab,limits=ylim) + facet_grid(.~variable,scale='free') + geom_errorbar(aes(ymin=lower, ymax=upper),width=.3,position=position_dodge(.9))"
                      
                      p2 <- .eval(p2,env=environment())
                      
                      if (!.require('gridExtra')) {
                        warning('you need the package gridExtra to make the plots printed in a single page!')
                        return(list(p1,p2))
                      } else {
                        return(.eval("grid.arrange(p1,p2)",env=environment()))
                      }
                    } else return(p1)
                  }
                } else {
                  drc <- data.frame(variable='a',Value=0)[0,]
                  for (nn in n) {
                    drc <- rbind(drc,data.frame(variable=nn,Value=x@response[[nn]][,1],x@response[[nn]][,2:ncol(x@response[[nn]])]))
                  }
                }
                
                p1 <- ".p1 <- ggplot(drc,aes(x=Value)) + geom_line(aes_string(y=colnames(drc)[3]),colour=col,size=lwd,linetype=lty) +scale_y_continuous(name = ylab,limits = c(0,1)) + facet_grid(.~variable,scales='free_x')
                for (nn in colnames(drc)[4:ncol(drc)]) .p1 <- .p1 + geom_line(aes_string(y=nn),colour=col,size=lwd,linetype=lty)"
                .eval(p1,env=environment())
                
                if (!is.null(nF)) {
                  drcc <- data.frame(Value=0,Response=0,variable='a')[0,]
                  for (nn in nF) {
                    .n <- length(x@response[[nn]][,1])
                    .m <- apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,mean,na.rm=TRUE)
                    .ci <- 1.96 * apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,sd,na.rm=TRUE) / sqrt(.n)
                    drcc <- rbind(drcc,data.frame(Value=x@response[[nn]][,1],Response=.m,lower=.m - .ci,upper=.m + .ci,variable=nn))
                  }
                  
                  if (confidence) p2 <- "ggplot(drcc,aes(x=Value,y=Response)) + geom_bar(stat = 'identity',fill=col)+ scale_y_continuous(name = ylab,limits=ylim) + facet_grid(.~variable,scale='free') + geom_errorbar(aes(ymin=lower, ymax=upper),width=.3,position=position_dodge(.9))"
                  else p2 <- "ggplot(drcc,aes(x=Value,y=Response)) + geom_bar(stat = 'identity',fill=col)+ scale_y_continuous(name = ylab,limits=ylim) + facet_grid(.~variable,scale='free')"
                  
                  p2 <- .eval(p2,env=environment())
                  
                  if (!.require('gridExtra')) {
                    warning('you need the package gridExtra to make the plots printed in a single page!')
                    return(list(p1,p2))
                  } else {
                    return(.eval("grid.arrange(p1,p2)",env=environment()))
                  }
                } else return(p1)
              } else {
                drc <- data.frame(Value=0,Response=0,variable='a')[0,]
                for (nn in n) {
                  colnames(x@response[[nn]]) <- c('Value','Response')
                  drc <- rbind(drc,data.frame(x@response[[nn]],variable=nn))
                }
                p1 <- "ggplot(drc,aes(x=Value,y=Response)) + geom_line(colour=col,size=lwd,linetype=lty) + facet_grid(.~variable,scales='free_x') +
                scale_y_continuous(limits=ylim) + labs(y = ylab,x = xlab) + ggtitle(main) +
                theme(axis.text=element_text(size=rel(cex.axis)),axis.title=element_text(size=rel(cex.lab)),plot.title = element_text(hjust = 0.5))"
                p1 <- .eval(p1,env=environment())
                
                if (!is.null(nF)) {
                  drcc <- data.frame(Variable=0,Response=0)[0,]
                  for (nn in nF) {
                    colnames(x@response[[nn]]) <- c('Variable','Response')
                    drcc <- rbind(drcc,data.frame(x@response[[nn]],variable=nn))
                  }
                  p2 <- "ggplot(drcc,aes(x=Variable,y=Response)) + geom_bar(stat = 'identity',fill=col) + facet_grid(.~variable,scale='free')"
                  p2 <- .eval(p2,env=environment())
                  if (!.require('gridExtra')) {
                    warning('you need the package gridExtra to make the plots printed in a single page!')
                    return(list(p1,p2))
                  } else {
                    pp <- .eval("grid.arrange(p1,p2)",env=environment())
                    return(pp)
                  }
                } else return(p1)
              }
              
            } else {
              np <- length(n) + length(nF)
              
              if (np > 16) {
                warning('Due to larger number of variables, only the plots for the first 16 variables are generated!')
                if (length(n) >= 16) {
                  n <- n[1:16]
                  nF <- NULL
                } else {
                  nF <- nF[1:(16 - length(n))]
                }
                par(mfrow=c(4,4),mar=c(5,4,1,1))
              } else {
                w <- floor(sqrt(np))
                h <- ceiling(np/w)
                if (abs(w-h) > 1) {
                  w <- w+1
                  h <- h-1
                }
                par(mfrow=c(w,h),mar=c(5,4,1,1))
              }
              #===========
              if (x@multi) {
                for (nn in n) {
                  if (mean) {
                    .m <- apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,mean,na.rm=TRUE)
                    plot(x@response[[nn]][,1],.m,type='l',xlab=nn,col=col,main=main,cex.axis=cex.axis,cex.lab=cex.lab,ylab=ylab,ylim=ylim,...)
                    if (confidence) {
                      .n <- length(x@response[[nn]][,1])
                      .ci <- 1.96 * apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,sd,na.rm=TRUE) / sqrt(.n)
                      lines(x@response[[nn]][,1],.m - .ci,col='gray',lty=2,lwd=lwd)
                      lines(x@response[[nn]][,1],.m + .ci,col='gray',lty=2,lwd=lwd)
                    }
                  } else {
                    plot(x@response[[nn]][,1:2],type='l',xlab=nn,col=col,main=main,cex.axis=cex.axis,cex.lab=cex.lab,ylab=ylab,ylim=ylim,...)
                    for (i in 3:ncol(x@response[[nn]])) lines(x@response[[nn]][,1],x@response[[nn]][,i],col=col,lwd=lwd)
                  }
                }
                
                if (!is.null(nF)) {
                  for (nn in nF) {
                    .m <- apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,mean,na.rm=TRUE)
                    
                    .bar <- barplot(.m,xlab=nn,col=col,main=main,cex.axis=cex.axis,cex.lab=cex.lab,ylab=ylab,ylim=ylim)
                    
                    if (confidence) {
                      .n <- length(x@response[[nn]][,1])
                      .ci <- 1.96 * apply(x@response[[nn]][,2:ncol(x@response[[nn]])],1,sd,na.rm=TRUE) / sqrt(.n)
                      
                      segments(.bar, .m - .ci, .bar,.m + .ci, lwd = 1.5)
                      
                      arrows(.bar, .m - .ci, .bar,.m + .ci, lwd = 1.5, angle = 90,code = 3, length = 0.05)
                    }
                  }
                }
                
              } else {
                for (nn in n) {
                  plot(x@response[[nn]],type='l',xlab=nn,col=col,main=main,cex.axis=cex.axis,cex.lab=cex.lab,ylab=ylab,ylim=ylim,...)
                }
                if (!is.null(nF)) {
                  for (nn in nF) {
                    barplot(x@response[[nn]][,2],xlab=nn,col=col,main=main,cex.axis=cex.axis,cex.lab=cex.lab,ylab=ylab,ylim=ylim)
                  }
                }
              }
            }
          }
)

#--------

setMethod("plot", signature(x='.nicheRaster'),
          function(x,y=NULL,gg=TRUE,xlab,ylab,col,cex.axis,cex.lab,main,...) {
            if (missing(gg)) gg <- .require('ggplot2')
            else if (gg && !.require('ggplot2')) gg <- FALSE
            
            if (missing(xlab)) xlab <- x@names[1]
            if (missing(ylab)) ylab <- x@names[2]
            
            if (missing(col)) col <- c('darkred','red','yellow','green','darkgreen','darkblue')
            
            if (missing(cex.axis)) cex.axis <- 0.8
            if (missing(cex.lab)) cex.lab <- 1
            
            if (missing(main)) main <- paste0('Ecological Niche described by: ',x@names[1],' - ',x@names[2])
            
            .lab1 <- round(as.vector(seq(x@scaleParams[1,1], x@scaleParams[2,1], length.out = 6)),1)
            .lab2 <- round(as.vector(seq(x@scaleParams[1,2], x@scaleParams[2,2], length.out = 6)),1)
            #--------
            if (gg) {
              drc <- as.data.frame(x@nicheRaster,xy=TRUE)
              p1 <- "ggplot(drc,aes(x=x,y=y,fill=niche)) +geom_raster() + coord_quickmap() + 
                    scale_y_continuous(breaks=seq(0, 1, length.out = 6),name = ylab,labels=.lab1) + scale_x_continuous(breaks=seq(0, 1, length.out = 6),name = xlab,labels=.lab2) + ggtitle(main) +
                    scale_fill_gradientn(colours=col,na.value='white') + theme_bw() +
                    theme(axis.text=element_text(size=rel(cex.axis)),axis.title=element_text(size=rel(cex.lab)),plot.title = element_text(hjust = 0.5))"
              p1 <- .eval(p1,env=environment())
              return(p1)
            } else {
              if (length(col) < 10) col <- colorRampPalette(col)(100)
              
              plot(x@nicheRaster,col=col,xaxt='n',yaxt='n',xlab=xlab,ylab=ylab,main=main,cex.axis=cex.axis,cex.lab=cex.lab,...)
              
              axis(1, at=seq(0, 1, length.out = 6), labels = FALSE)
              text(seq(0, 1, length.out = 6),par("usr")[3] - 0.05, labels = .lab1, srt = 0, pos = 1, xpd = TRUE,cex=cex.axis)
              
              axis(2, at=seq(0, 1, length.out = 6), labels = FALSE)
              text(par("usr")[1]-0.1,seq(0, 1, length.out = 6)+0.05, labels = .lab2, srt = 0, pos = 1, xpd = TRUE,cex=cex.axis)
            }
          }
)

#--------
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