# Author: Babak Naimi, naimi.b@gmail.com
# Date:  June 2018
# Last update:  June 2018
# Version 1.0
# Licence GPL v3

.getModel.info <- function(x,w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL) {
  if (missing(w) || is.null(w)) {
    
    w1 <- w2 <- w3 <- w4 <- TRUE
    if (!is.null(species)) {
      species <- .pmatch(species,unique(as.character(x@run.info[,2])))
      w1 <- x@run.info[,2] %in% species
    }
    
    if (!is.null(method)) {
      method <- .methodFix(method)
      w2 <- x@run.info[,3] %in% method
    }
    
    
    if (!is.null(replication)) {
      if (length(x@replicates) != 0) {
        replication <- .replicate.methodFix(replication)
        w3 <- x@run.info[,4] %in% replication
      }
    }
    
    if (!is.null(run)) {
      if (!is.null(run)) {
        if (length(x@replicates) != 0) {
          r <- unlist(lapply(x@replicates[[1]],function(x) x$method))
          ru <- unique(r)
          names(ru) <- ru
          rID <- lapply(ru,function(x) which(r == x))
          w4 <- c()
          for (i in 1:length(rID)) {
            w4 <- c(w4,rID[[i]][c(1:length(rID[[i]])) %in% run])
          }
          
          w4 <- x@run.info[,5] %in% w4
        }
      }
    }
    x@run.info[w1 & w2 & w3 & w4,]
  } else x@run.info[x@run.info[,1] %in% w,]
}

#--------


.getModel.info2 <- function(x,w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,wtest=NULL) {
  # comparing to .getModel.info: In this, only one species is allowed!
  # x: sdmModels
  if (!is.null(w)) {
    mi <- .getModel.info(x,w)
  } else {
    mi <- .getModel.info(x)
    u <- as.character(unique(mi[,2]))
    m <- as.character(unique(mi[,3]))
    r <- unique(mi[,4])
    if (!is.null(species)) {
      if (length(species) > 1) {
        species <- species[1]
        warning('only the first species is considered!')
      }
      if (is.numeric(species)) {
        if (length(u) <= species) species <- u[species]
        else stop('The specified species is not recognised!')
      } else {
        species <- .pmatch(species,u)
        if (is.na(species)) stop('The specified species is not recognised!')
      }
    } else {
      if (length(u) > 1) stop('This object contains models for more than one species; in species argument spcify the name of species!')
      else species <- u
    }
    
    if (!is.null(method)) {
      method <- .sdmMethods$fixNames(method)
      wm <- method %in% m
      if (any(!wm)) {
        if (all(!wm)) stop('the specified methods do not exist in the object!')
        warning(paste('Methods',paste(method[!wm],collapse=', '),'do not exsit in the object, and are excluded!'))
        method <- method[wm]
      }
    } else method <- m
    
    if (!is.null(replication)) replication <- .replicate.methodFix(replication)
    else replication <- r
    
    mi <- .getModel.info(x,species=species,method=method,replication=replication,run=run)
  }
  
  mi
}
#---------

if (!isGeneric("getModelInfo")) {
  setGeneric("getModelInfo", function(x,...)
    standardGeneric("getModelInfo"))
}


setMethod('getModelInfo', signature(x='sdmModels'), 
          function(x,w,...) {
            if (missing(w)) w <- NULL
            .getModel.info(x,w,...)
            
          }
)
#-----------

if (!isGeneric("getModelId")) {
  setGeneric("getModelId", function(x,success,species,method,replication,run)
    standardGeneric("getModelId"))
}


setMethod('getModelId', signature(x='sdmModels'), 
          function(x,success=TRUE,species=NULL,method=NULL,replication=NULL,run=NULL) {
            if (missing(success)) success <- TRUE
            if (missing(species)) species <- NULL
            if (missing(method)) method <- NULL
            if (missing(replication)) replication <- NULL
            if (missing(run)) run <- NULL
            
            mi <- .getModel.info(x,w=NULL,species=species,method=method,replication=replication,run=run)
            
            if (nrow(mi) == 0) return(NULL)
            
            if (success) mi <- mi[mi$success,]
            
            if (nrow(mi) > 0) return(mi$modelID)
            else return(NULL)
          }
)
