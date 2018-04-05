# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2016
# Version 2.7
# Licence GPL v3


# for internal use in .workLoad
# pv is pred.par
._varImp <- function(pv,pred,sp,nsim=5) {
  # if the datatype is different than sdmDataFrame, then it should be updated to support...
  ww <- names(pv[[2]])[which(names(pv[[2]]) != sp)]
  dd <- pv[[2]]
  obs <- pv[[2]][,sp]
  d1 <- pred(pv)
  vi <- vj <- rep(NA,nsim)
  varImp1 <- varImp2 <- rep(NA,length(ww))
  names(varImp1) <- names(varImp2) <- ww
  a1 <- .auc(obs,d1)
  for (v in ww) {
    for (i in 1:nsim) {
      pv[[2]][,v] <- dd[sample(nrow(dd)),v]
      d2 <- pred(pv)
      cr <- cor(d1,d2,use="complete.obs")
      if (cr < 0) cr <- 0
      vi[i] <- 1 - cr
      a2 <- .auc(obs,d2)
      a2 <- (a1-a2)*2
      if (a2 > 1) a2 <- 1
      else if (a2 < 0) a2 <- 0
      vj[i] <- a2
    }
    varImp1[v] <- round(mean(vi,na.rm=TRUE),4)
    varImp2[v] <- round(mean(vj,na.rm=TRUE),4)
    pv[[2]] <- dd
  }
  new('.varImportance',variables=ww,varImportance=data.frame(variables=ww,corTest=varImp1,AUCtest=varImp2))
}
#---------

.getVarImpObject <- function(x,id,wtest) {
  # stat can be 1 (threshold-independent) OR 2 (threshold-dependent)
  mi <- x@run.info
  w <- which(mi$modelID == id)
  if (length(w) == 1) {
    if (missing(wtest) || is.null(wtest)) wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
    else {
      wtest <- .pmatch(wtest,c('training','test.dep','test.indep'))[1]
      if (is.na(wtest)) wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
    }
    
    sp <- as.character(mi$species)[w]
    mo <- as.character(mi$method)[w]
    i <- as.character(mi$modelID)[w]
    x@models[[sp]][[mo]][[i]]@varImportance[[wtest]]
  } 
  
}



#--------
if (!isGeneric("getVarImp")) {
  setGeneric("getVarImp", function(x,id, wtest,...)
    standardGeneric("getVarImp"))
}  

setMethod('getVarImp', signature(x='sdmModels'),
          function(x, id, wtest,...) {
            if (missing(id)) id <- NULL
            if (missing(wtest)) wtest <- NULL
            
            if (!is.null(id) & length(id) == 1) {
              .getVarImpObject(x,id,wtest)
            } else {
              stop('This version only support extracting the variable importance for 1 model each time!')
            }
            
          }
)
