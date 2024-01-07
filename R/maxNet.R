# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Jan 2024
# last update: Jan 2024
# Version 1.1
# Licence GPL v3
#-------------


.maxNet <- function(formula,data,f=NULL,regfun = NULL ,regmult = 1,addsamplestobackground=TRUE,...) {
  nsp <- deparse(formula[[2]])
  p <- data[,nsp]
  x <- data[,colnames(data) != nsp]
  
  if (missing(f) || is.null(f)) f <- .eval('maxnet.formula(p, x)',environment())
  if (missing(regfun) || is.null(regfun)) regfun <- .eval('maxnet.default.regularization',environment())
  if (missing(regmult) || is.null(regmult)) regmult <- 1
  if (missing(addsamplestobackground) || is.null(addsamplestobackground)) addsamplestobackground <- TRUE
  
  .maxnet <- .eval('maxnet',environment())
  .maxnet(p=p,data=x,f=f,regfun=regfun,regmult=regmult,addsamplestobackground=addsamplestobackground,...)
}
