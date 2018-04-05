# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Nov. 2016
# Version 1.1
# Licence GPL v3

#-------------
.gauss <- function(x,c,s) {
  exp((-(x-c)^2)/(2*s*s))
}
#------
.bioclimFit <- function(formula,data,c=2,weights=NULL,...) {
  varnames <- all.vars(formula)
  nsp <- deparse(formula[[2]])
  w <- data[,nsp] == 1
  bio <- new('.bioclimModel')
  x <- data[w,colnames(data) != nsp]
  nFact <- .where(is.factor,x)
  if (any(nFact)) x <- x[,-which(nFact),drop=FALSE]
  if (ncol(x) < 2) stop('At least two continous variables are needed to fit the model!')
  bio@features <- colnames(x)
  bio@min <- apply(x,2,min,na.rm=TRUE)
  bio@max <- apply(x,2,max,na.rm=TRUE)
  bio@median <- apply(x,2,median,na.rm=TRUE)
  bio@q25 <- apply(x,2,quantile,probs=0.25,na.rm=TRUE)
  bio@q75 <- apply(x,2,quantile,probs=0.75,na.rm=TRUE)
  bio@c <- c
  bio@weights <- weights
  bio
}

setMethod('predict', signature(object='.bioclimModel'), 
          function(object, newdata,...) {
            if (!all(object@features %in% colnames(newdata))) stop('One or more variables in the model do not exist in the data!')
            newdata <- newdata[,object@features]
            out <- matrix(nrow=nrow(newdata),ncol=ncol(newdata))
            s <- (object@q75 - object@q25) / object@c
            
            for (i in seq_along(object@features)) {
              out[,i] <- .gauss(newdata[,i],object@median[i],s[i])
            }
            w0 <- apply(out,1,function(x) any(x < 0.01))
            w0 <- ifelse(w0,0,1)
            if (is.null(object@weights) || length(object@weights) != length(object@features)) out <- apply(out,1,mean,na.rm=TRUE)
            else {
              object@weights <- object@weights / sum(object@weights)
              out <- apply(out,1,function(x) sum(x*object@weights))
            }
            out * w0
          }
)
###############
# Bioclimatic based on the package dismo:


.bioclimDismo <- function(formula,data,...) {
  nsp <- deparse(formula[[2]])
  w <- data[,nsp] == 1
  x <- data[w,colnames(data) != nsp]
  nFact <- .where(is.factor,x)
  if (any(nFact)) x <- x[,-which(nFact),drop=FALSE]
  if (ncol(x) < 2) stop('At least two continous variables are needed to fit the model!')
  dismo::bioclim(x=x,...)
}
