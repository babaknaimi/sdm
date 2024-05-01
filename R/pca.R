# Author: Babak Naimi, naimi.b@gmail.com
# Last Update :  April 2024
# Version 1.2
# Licence GPL v3
#---------------------

if (!isGeneric("pca")) {
  setGeneric("pca", function(x,scale,filename,...)
    standardGeneric("pca"))
}


setMethod('pca', signature(x='sdmdata'), 
          function(x,scale=FALSE,filename="",...) {
            if (missing(scale) || !is.logical(scale)) scale <- FALSE
            if (missing(filename)) filename <- ''
            n <- x@sdmFormula@vars@numeric$names # continuous variables
            
            if (length(n < 2)) stop('The input data, x, should have at least 2 numeric variables...!')
            
            x <- as.data.frame(x)[,n]
            
            if (scale) x <- scale(x)
            
            .p <- princomp(x,...)
            
            if (filename != "" && is.character(filename) ) {
              if (grepl('.csv$',filename)) write.csv(.p$scores,file=filename,row.names = FALSE)
              else warning('filename is ignored (it should be a CSV filename, e.g., filename="xxx.csv")')
            }
            
            new('.pcaObject',data=as.data.frame(.p$scores),pcaObject=.p,scaled=scale)
          }
)
#------
setMethod('pca', signature(x='data.frame'), 
          function(x,scale=FALSE,filename="",...) {
            
            if (missing(scale) || !is.logical(scale)) scale <- FALSE
            
            if (missing(filename)) filename <- ''
            
            w <- .where(is.factor,x)
            
            if (any(w)) {
              n <- names(w)[which(!w)]
              if (length(n > 1)) x <- x[,n]
              else stop('The input data, x, should have at least 2 numeric variables...!')
              
              warning('The factor variable(s) are excluded from the input data.frame (x)...!')
              
            }
            if (scale) x <- scale(x)
            .p <- princomp(x,...)
            
            if (filename != "" && is.character(filename) ) {
              if (grepl('.csv$',filename)) write.csv(x,file=filename,row.names = FALSE)
              else warning('filename is ignored (it should be a CSV filename, e.g., filename="xxx.csv")')
            }
              
            
            new('.pcaObject',data=as.data.frame(.p$scores),pcaObject=.p,scaled=scale)
          }
)
#----

setMethod('pca', signature(x='RasterStackBrick'), 
          function(x,scale=FALSE,filename="",...) {
            if (missing(scale) || !is.logical(scale)) scale <- FALSE
            if (missing(filename)) filename <- ''
            n <- nlayers(x)
            if (n < 2) stop('The input data, x, should have at least 2 numeric layers...!')
            if (scale) x <- scale(x)
            nc <- 1:ncell(x)
            nc <- nc[which(!is.na(x[[1]][]))]
            d <- x[nc]
            d <- princomp(d,...)
            for (i in 1:n) x[[i]][nc] <- d$scores[,i]
            names(x) <- colnames(d$scores)
            if (filename != "" && is.character(filename)) x <- writeRaster(x,filename=filename)
            new('.pcaObject',data=x,pcaObject=d,scaled=scale)
          }
)
#------

setMethod('pca', signature(x='SpatRaster'), 
          function(x,scale=FALSE,filename="",...) {
            if (missing(scale) || !is.logical(scale)) scale <- FALSE
            if (missing(filename)) filename <- ''
            n <- nlyr(x)
            if (n < 2) stop('The input data, x, should have at least 2 numeric layers...!')
            if (scale) x <- scale(x)
            nc <- cells(x)
            d <- x[nc]
            d <- princomp(d,...)
            for (i in 1:n) x[[i]][nc] <- d$scores[,i]
            names(x) <- colnames(d$scores)
            if (filename != "" && is.character(filename)) x <- writeRaster(x,filename=filename)
            new('.pcaObject',data=x,pcaObject=d,scaled=scale)
          }
)
