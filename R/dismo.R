# Author: Babak Naimi, naimi.b@gmail.com
# Date :  April 2018
# last update: April 2018
# Version 1.0
# Licence GPL v3

#-------------


# Profile-based methods used from the package dismo:


.bioclimDismo <- function(formula,data,...) {
  nsp <- deparse(formula[[2]])
  w <- data[,nsp] == 1
  x <- data[w,colnames(data) != nsp]
  nFact <- .where(is.factor,x)
  if (any(nFact)) x <- x[,-which(nFact),drop=FALSE]
  if (ncol(x) < 2) stop('At least two continous variables are needed to fit the model!')
  dismo::bioclim(x=x,...)
}

#-----------

.domainDismo <- function(formula,data,...) {
  nsp <- deparse(formula[[2]])
  w <- data[,nsp] == 1
  x <- data[w,colnames(data) != nsp]
  nFact <- .where(is.factor,x)
  if (any(nFact)) x <- x[,-which(nFact),drop=FALSE]
  if (ncol(x) < 2) stop('At least two continous variables are needed to fit the model!')
  dismo::domain(x=x,...)
}

#----------

.mahalDismo <- function(formula,data,...) {
  nsp <- deparse(formula[[2]])
  w <- data[,nsp] == 1
  x <- data[w,colnames(data) != nsp]
  nFact <- .where(is.factor,x)
  if (any(nFact)) x <- x[,-which(nFact),drop=FALSE]
  if (ncol(x) < 2) stop('At least two continous variables are needed to fit the model!')
  dismo::mahal(x=x,...)
}
