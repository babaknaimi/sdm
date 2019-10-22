# Author: Babak Naimi, naimi.b@gmail.com
# Date:  June 2018
# Last update:  Oct. 2019
# Version 1.2
# Licence GPL v3


.factorFreq <- function(x) {
  x <- factor(x)
  u <- levels(x)
  names(u) <- u
  unlist(lapply(u,function(i) length(which(x == i))))
}



.responseCurve1 <- function(m,id,si=100,includeTest=FALSE) {
  mi <- m@run.info[m@run.info$modelID == id,]
  if (nrow(mi) != 1) stop('the model with the specified id does not exist!')
  
  if (!mi$success) stop('the model for the specified id did not fit successfully!')
  #--------
  rc <- new('.responseCurve')
  
  sp <- as.character(mi$species)
  
  
  if (includeTest) dt <- as.data.frame(m@data,sp=sp)
  else dt <- as.data.frame(m@data,grp=c('train'),sp=sp)
  
  nf <- m@setting@featuresFrame@vars
  
  rc@variables <- nf
  
  dt <- dt[,nf,drop=FALSE]
  
  nFact <- m@data@factors
  
  if (any(nFact %in% nf)) {
    nFact <- nFact[nFact %in% nf]
    nf <- .excludeVector(nf,nFact)
  } else nFact <- NULL
  
  rc@categorical <- nFact
  
  if (!is.null(nFact)) {
    dm <- data.frame(matrix(nrow=1,ncol=length(nf)+length(nFact)))
    colnames(dm) <- c(nf,nFact)
    
    for (n in nf) dm[,n] <- mean(dt[,n],na.rm=TRUE)
    
    for (n in nFact) {
      u <- .factorFreq(dt[,n])
      dm[,n] <- names(which.max(u)[1])
    }
    
    for (n in nFact) {
      dt[,n] <- factor(dt[,n])
      u <- unique(dt[,n])
      dc <- data.frame(matrix(nrow=length(u),ncol=2)) # the results of response curve for categorical variable
      dc[,1] <- u
      colnames(dc) <- c(n,'Response')
      for (i in seq_along(u)) {
        dv <- dm
        dv[,n] <- u[i]
        p <- predict(m,newdata=dv,w=id)
        dc[i,1] <- u[i]
        dc[i,2] <- p[1,1]
      }
      rc@response[[n]] <- dc
    }
    #------
    
    
    if (!is.null(nf)) {
      dc <- dm
      dm <- data.frame(matrix(nrow=si,ncol=length(nf)+length(nFact)))
      colnames(dm) <- c(nf,nFact)
      
      for (n in c(nf,nFact)) dm[,n] <- dc[,n]
      
      for (n in nf) {
        dv <- dm
        dv[,n] <- seq(min(dt[,n],na.rm=TRUE),max(dt[,n],na.rm=TRUE),length.out = si)
        p <- predict(m,newdata=dv,w=id)
        rc@response[[n]] <- data.frame(dv[,n],p[,1])
        colnames(rc@response[[n]]) <- c(n,'Response')
      }
    }
  } else {
    dm <- data.frame(matrix(nrow=si,ncol=length(nf)))
    colnames(dm) <- nf
    
    for (n in nf) dm[,n] <- mean(dt[,n],na.rm=TRUE)
    
    for (n in nf) {
      dv <- dm
      dv[,n] <- seq(min(dt[,n],na.rm=TRUE),max(dt[,n],na.rm=TRUE),length.out = si)
      p <- predict(m,newdata=dv,w=id)
      rc@response[[n]] <- data.frame(dv[,n],p[,1])
      colnames(rc@response[[n]]) <- c(n,'Response')
    }
  }
  
  rc
}
#-----------


.responseCurve <- function(m,id,si=100,includeTest=FALSE) {
  
  mi <- m@run.info[m@run.info$modelID %in% id,]
  
  mi <- mi[mi$success,]
  
  na <- paste0(mi$method,'_ID-',mi$modelID)
  #--------
  rc <- new('.responseCurve')
  
  rc@multi <- nrow(mi) > 1
  
  sp <- as.character(mi$species)
  
  
  if (includeTest) dt <- as.data.frame(m@data,sp=sp)
  else dt <- as.data.frame(m@data,grp=c('train'),sp=sp)
  
  nf <- m@setting@featuresFrame@vars
  
  rc@variables <- nf
  
  dt <- dt[,nf,drop=FALSE]
  
  nFact <- m@data@factors
  
  if (any(nFact %in% nf)) {
    nFact <- nFact[nFact %in% nf]
    nf <- .excludeVector(nf,nFact)
  } else nFact <- NULL
  
  rc@categorical <- nFact
  
  if (!is.null(nFact)) {
    dm <- data.frame(matrix(nrow=1,ncol=length(nf)+length(nFact)))
    colnames(dm) <- c(nf,nFact)
    
    for (n in nf) dm[,n] <- mean(dt[,n],na.rm=TRUE)
    
    for (n in nFact) {
      u <- .factorFreq(dt[,n])
      dm[,n] <- names(which.max(u)[1])
      dm[,n] <- factor(dm[,n],levels = levels(dt[,n]))
    }
    
    for (n in nFact) {
      #dt[,n] <- factor(dt[,n])
      u <- unique(dt[,n])
      dc <- data.frame(matrix(nrow=length(u),ncol=1+nrow(mi))) # the results of response curve for categorical variable
      dc[,1] <- u
      colnames(dc) <- c(n,na)
      for (i in seq_along(u)) {
        dv <- dm
        dv[,n] <- u[i]
        p <- predict(m,newdata=dv,w=mi$modelID)
        dc[i,1] <- u[i]
        for (j in 2:ncol(dc)) dc[i,j] <- p[1,(j-1)]
      }
      
      rc@response[[n]] <- dc
    }
    #------
    
    
    if (!is.null(nf)) {
      dc <- dm
      dm <- data.frame(matrix(nrow=si,ncol=length(nf)+length(nFact)))
      colnames(dm) <- c(nf,nFact)
      
      for (n in c(nf,nFact)) dm[,n] <- dc[,n]
      
      for (n in nf) {
        dv <- dm
        dv[,n] <- seq(min(dt[,n],na.rm=TRUE),max(dt[,n],na.rm=TRUE),length.out = si)
        dc <- data.frame(matrix(nrow=si,ncol=1+nrow(mi)))
        dc[,1] <- dv[,n]
        p <- predict(m,newdata=dv,w=mi$modelID)
        for (j in 2:ncol(dc)) dc[,j] <- p[,(j-1)]
        
        rc@response[[n]] <- dc
        colnames(rc@response[[n]]) <- c(n,na)
      }
    }
  } else {
    dm <- data.frame(matrix(nrow=si,ncol=length(nf)))
    colnames(dm) <- nf
    
    for (n in nf) dm[,n] <- mean(dt[,n],na.rm=TRUE)
    
    for (n in nf) {
      dv <- dm
      dv[,n] <- seq(min(dt[,n],na.rm=TRUE),max(dt[,n],na.rm=TRUE),length.out = si)
      dc <- data.frame(matrix(nrow=si,ncol=1+nrow(mi)))
      dc[,1] <- dv[,n]
      for (j in seq_along(mi$modelID)) {
        p <- predict(m,newdata=dv,w=mi$modelID[j])
        dc[,(j+1)] <- p[,1]
      }
      
      rc@response[[n]] <- dc
      colnames(rc@response[[n]]) <- c(n,na)
    }
  }
  
  rc
}
#-----------


if (!isGeneric("getResponseCurve")) {
  setGeneric("getResponseCurve", function(x,id,...)
    standardGeneric("getResponseCurve"))
}


setMethod("getResponseCurve", signature(x='sdmModels'),
          function(x,id,size=100,includeTest=FALSE,...) {
            
            if (missing(id)) {
              id <- x@run.info[which(x@run.info$success),1]
              if (length(id) == 1) cat(paste0('The id argument is not specified; id = ',id,' is considered. \n'))
              else if (length(id) > 1) cat(paste0('The id argument is not specified; The modelIDs of ', length(id),' successfully fitted models are assigned to id...! \n'))
              else stop('It seems that no successfully fitted models do exist in the sdmModels object!')
            }
            
            if (missing(size)) size <- 100
            
            if (missing(includeTest)) includeTest <- FALSE
            
            .responseCurve(x,id=id,si=size,includeTest=includeTest)
          }
)


if (!isGeneric("rcurve")) {
  setGeneric("rcurve", function(x,n,id,mean,confidence,gg,...)
    standardGeneric("rcurve"))
}  


setMethod("rcurve", signature(x='sdmModels'),
          function(x,n,id,mean,confidence,gg,size,includeTest,...) {
            
            # if (missing(id)) {
            #   id <- x@run.info[which(x@run.info$success)[1],1]
            #   if (length(which(x@run.info$success)) > 1) cat(paste0('The id is missing; the first successfully fitted model is considered, i.e., id = ',id,'\n'))
            # }
            if (missing(gg)) gg <- .require('ggplot2')
            else if (gg && !.require('ggplot2')) gg <- FALSE
            
            if (missing(id)) {
              id <- x@run.info[which(x@run.info$success),1]
              if (length(id) == 1) cat(paste0('The id argument is not specified; id = ',id,' is considered. \n'))
              else if (length(id) > 1) cat(paste0('The id argument is not specified; The modelIDs of ', length(id),' successfully fitted models are assigned to id...! \n'))
              else stop('It seems that no successfully fitted models do exist in the sdmModels object!')
            }
            
            if (missing(size)) size <- 100
            
            if (missing(includeTest)) includeTest <- FALSE
            
            if (missing(mean)) mean <- TRUE
            
            if (missing(confidence)) {
              if (length(id) == 1) confidence <- FALSE
              else confidence <- TRUE
            } else {
              if (is.logical(confidence)) {
                if (confidence && length(id) == 1) {
                  confidence <- FALSE
                  warning('confidence intervals can only be generated if modelIDs of more than one model are selected in the id argument!')
                }
              }
            }
            
            if (missing(n)) n <- NULL
            
            rc <- .responseCurve(x,id=id,si=size,includeTest=includeTest)
            
            plot(rc,y=n,gg=gg,mean=mean,confidence=confidence,...)
          }
)



setMethod("rcurve", signature(x='.responseCurve'),
          function(x,n,id,mean=TRUE,confidence=TRUE,gg=TRUE,...) {
            
            if (missing(gg)) gg <- .require('ggplot2')
            else if (gg && !.require('ggplot2')) gg <- FALSE
            
            if (missing(mean)) mean <- TRUE
            
            if (missing(confidence)) confidence <- TRUE
            
            if (missing(n)) n <- NULL
            
            plot(x,y=n,gg=gg,mean=mean,confidence=confidence,...)
          }
)
