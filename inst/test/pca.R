# Author: Babak Naimi, naimi.b@gmail.com
# Date :  September 2015
# Version 1.0
# Licence GPL v3

if (!isGeneric("pca")) {
  setGeneric("pca", function(x,...)
    standardGeneric("pca"))
}


setMethod('pca', signature(x='singleSpecies'), 
          function(x,...) {
            n <- x@train@Features@featureNames
            x <- as.data.frame(x)[,n]
            d <- princomp(x,...)
            for (i in 1:length(n)) x[,i] <- d$scores[,i]
            colnames(x) <- colnames(d$scores)
            list(data=x,pca=d)
          }
)
          

setMethod('pca', signature(x='multipleSpecies'), 
          function(x,...) {
            n <- x@train@Features@featureNames
            x <- as.data.frame(x)[,n]
            d <- princomp(x,...)
            for (i in 1:length(n)) x[,i] <- d$scores[,i]
            colnames(x) <- colnames(d$scores)
            list(data=x,pca=d)
          }
)



setMethod('pca', signature(x='RasterStackBrick'), 
          function(x,...) {
            n <- nlayers(x)
            nc <- 1:ncell(x)
            nc <- nc[which(!is.na(x[[1]][]))]
            d <- x[nc]
            d <- princomp(d,...)
            for (i in 1:n) x[[i]][nc] <- d$scores[,i]
            names(x) <- colnames(d$scores)
            list(data=x,pca=d)
          }
)