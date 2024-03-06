# Author: Babak Naimi, naimi.b@gmail.com
# Date:  June 2018
# Last update:  March 2024
# Version 1.5
# Licence GPL v3


.factorFreq <- function(x) {
  x <- factor(x)
  u <- levels(x)
  names(u) <- u
  unlist(lapply(u,function(i) length(which(x == i))))
}



.responseCurve1 <- function(m,id,si=100,includeTest=FALSE,.fun=base::mean) {
  if (missing(.fun) || !is.function(.fun)) .fun <- base::mean
  
  mi <- m@run.info[m@run.info$modelID == id,]
  if (nrow(mi) != 1) stop('the model with the specified id does not exist!')
  
  if (!mi$success) stop('the model for the specified id did not fit successfully!')
  #--------
  rc <- new('.responseCurve')
  
  sp <- as.character(mi$species)
  
  
  if (includeTest) dt <- as.data.frame(m@data,sp=sp)
  else dt <- as.data.frame(m@data,grp=c('train'),sp=sp)
  
  .wPr <- which(dt[,sp] == 1)
  
  nf <- m@setting@featureFrame@predictors
  
  rc@variables <- nf
  
  dt <- dt[,nf,drop=FALSE]
  
  nFact <- names(m@setting@featureFrame@categorical)
  
  if (any(nFact %in% nf)) {
    nFact <- nFact[nFact %in% nf]
    nf <- .excludeVector(nf,nFact)
  } else nFact <- NULL
  
  rc@categorical <- nFact
  
  if (!is.null(nFact)) {
    dm <- data.frame(matrix(nrow=1,ncol=length(nf)+length(nFact)))
    colnames(dm) <- c(nf,nFact)
    
    for (n in nf) dm[,n] <- .fun(dt[.wPr,n],na.rm=TRUE)
    
    for (n in nFact) {
      u <- .factorFreq(dt[.wPr,n])
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
        p <- try(predict(m,newdata=dv,id=id),silent = TRUE)
        if (!inherits(p,'try-error')) {
          dc[i,2] <- p[1,1]
        }
        dc[i,1] <- u[i]
        
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
        
        p <- try(predict(m,newdata=dv,id=id),silent = TRUE)
        if (!inherits(p,'try-error')) {
          rc@response[[n]] <- data.frame(dv[,n],p[,1])
          colnames(rc@response[[n]]) <- c(n,'Response')
        }
      }
    }
  } else {
    dm <- data.frame(matrix(nrow=si,ncol=length(nf)))
    colnames(dm) <- nf
    
    for (n in nf) dm[,n] <- .fun(dt[.wPr,n],na.rm=TRUE)
    
    for (n in nf) {
      dv <- dm
      dv[,n] <- seq(min(dt[,n],na.rm=TRUE),max(dt[,n],na.rm=TRUE),length.out = si)
      p <- try(predict(m,newdata=dv,id=id),silent = TRUE)
      if (!inherits(p,'try-error')) {
        rc@response[[n]] <- data.frame(dv[,n],p[,1])
        colnames(rc@response[[n]]) <- c(n,'Response')
      }
      
    }
  }
  
  rc
}
#-----------

.responseCurveEns <- function(m,id,si=100,includeTest=FALSE,.fun=base::mean,setting=list(method='weighted',stat='auc')) {
  
  if (missing(setting) || is.null(setting)) setting <- list(method='weighted',stat='auc')
  
  if (missing(.fun) || !is.function(.fun)) .fun <- base::mean
  
  mi <- m@run.info
  if ('id' %in% names(setting)) {
    id <- setting[['id']]
    mi <- mi[mi$modelID %in% id,]
  } else id <-mi$modelID
  
  mi <- mi[mi$success,]
  
  
  sp <- unique(as.character(mi$species))
  
  if (length(sp) > 1) {
    warning('More than one species exists in the model object; only first species is considered!')
    sp <- sp[1]
    mi <- mi[mi$species == sp,]
    id <-mi$modelID
    setting[['id']] <- id
  }
  
  if (nrow(mi) < 2) stop('To get Response curve from Ensemble of models, at least two (successfully fitted) models are required in the sdmModels object!')
  
  na <- paste0(mi$method,'_ID-',mi$modelID)
  #--------
  rc <- new('.responseCurve',multi=FALSE)
  
  if (includeTest) dt <- as.data.frame(m@data,sp=sp)
  else dt <- as.data.frame(m@data,grp=c('train'),sp=sp)
  
  .wPr <- which(dt[,sp] == 1)
  
  nf <- m@setting@featureFrame@predictors
  
  rc@variables <- nf
  
  dt <- dt[,nf,drop=FALSE]
  
  nFact <- names(m@setting@featureFrame@categorical)
  
  if (any(nFact %in% nf)) {
    nFact <- nFact[nFact %in% nf]
    nf <- .excludeVector(nf,nFact)
  } else nFact <- NULL
  
  rc@categorical <- nFact
  
  if (!is.null(nFact)) {
    dm <- data.frame(matrix(nrow=1,ncol=length(nf)+length(nFact)))
    colnames(dm) <- c(nf,nFact)
    
    for (n in nf) dm[,n] <- .fun(dt[.wPr,n],na.rm=TRUE)
    
    for (n in nFact) {
      u <- .factorFreq(dt[.wPr,n])
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
        p <- try(ensemble(m,newdata=dv,setting=setting),silent = TRUE)
        if (!inherits(p,'try-error')) {
          dc[i,2] <- p[1,1]
        }
        dc[i,1] <- u[i]
        
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
        
        p <- try(ensemble(m,newdata=dv,setting=setting),silent = TRUE)
        if (!inherits(p,'try-error')) {
          rc@response[[n]] <- data.frame(dv[,n],p[,1])
          colnames(rc@response[[n]]) <- c(n,'Response')
        }
      }
    }
  } else {
    dm <- data.frame(matrix(nrow=si,ncol=length(nf)))
    colnames(dm) <- nf
    
    for (n in nf) dm[,n] <- .fun(dt[.wPr,n],na.rm=TRUE)
    
    for (n in nf) {
      dv <- dm
      dv[,n] <- seq(min(dt[,n],na.rm=TRUE),max(dt[,n],na.rm=TRUE),length.out = si)
      p <- try(ensemble(m,newdata=dv,setting=setting),silent = TRUE)
      if (!inherits(p,'try-error')) {
        rc@response[[n]] <- data.frame(dv[,n],p[,1])
        colnames(rc@response[[n]]) <- c(n,'Response')
      }
      
    }
  }
  
  rc
}
#-----------

.responseCurve <- function(m,id,si=100,includeTest=FALSE,.fun=base::mean) {
  if (missing(.fun) || !is.function(.fun)) .fun <- base::mean
  
  mi <- m@run.info[m@run.info$modelID %in% id,]
  
  mi <- mi[mi$success,]
  
  sp <- unique(as.character(mi$species))
  
  if (length(sp) > 1) {
    warning('More than one species exists in the model object; only first species is considered!')
    sp <- sp[1]
    mi <- mi[mi$species == sp,]
    id <-mi$modelID
  }
  
  na <- paste0(mi$method,'_ID-',mi$modelID)
  #--------
  rc <- new('.responseCurve')
  
  rc@multi <- nrow(mi) > 1
  
  if (includeTest) dt <- as.data.frame(m@data,sp=sp)
  else dt <- as.data.frame(m@data,grp=c('train'),sp=sp)
  
  .wPr <- which(dt[,sp] == 1)
  
  nf <- m@setting@featureFrame@predictors
  
  rc@variables <- nf
  
  dt <- dt[,nf,drop=FALSE]
  
  nFact <- names(m@setting@featureFrame@categorical)
  
  if (any(nFact %in% nf)) {
    nFact <- nFact[nFact %in% nf]
    nf <- .excludeVector(nf,nFact)
  } else nFact <- NULL
  
  rc@categorical <- nFact
  
  if (!is.null(nFact)) {
    dm <- data.frame(matrix(nrow=1,ncol=length(nf)+length(nFact)))
    colnames(dm) <- c(nf,nFact)
    
    for (n in nf) dm[,n] <- .fun(dt[.wPr,n],na.rm=TRUE)
    
    for (n in nFact) {
      u <- .factorFreq(dt[.wPr,n])
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
        p <- try(predict(m,newdata=dv,id=mi$modelID),silent = TRUE)
        if (!inherits(p,'try-error')) {
          for (j in 1:ncol(p)) dc[,j+1] <- p[,j]
        }
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
        p <- try(predict(m,newdata=dv,id=mi$modelID),silent = TRUE)
        if (!inherits(p,'try-error')) {
          for (j in 1:ncol(p)) dc[,j+1] <- p[,j]
        }
        
        rc@response[[n]] <- dc
        colnames(rc@response[[n]]) <- c(n,na)
      }
    }
  } else {
    dm <- data.frame(matrix(nrow=si,ncol=length(nf)))
    colnames(dm) <- nf
    
    for (n in nf) dm[,n] <- .fun(dt[.wPr,n],na.rm=TRUE)
    
    for (n in nf) {
      dv <- dm
      dv[,n] <- seq(min(dt[,n],na.rm=TRUE),max(dt[,n],na.rm=TRUE),length.out = si)
      dc <- data.frame(matrix(nrow=si,ncol=1+nrow(mi)))
      dc[,1] <- dv[,n]
      for (j in seq_along(mi$modelID)) {
        p <- try(predict(m,newdata=dv,id=mi$modelID[j]),silent = TRUE)
        if (!inherits(p,'try-error')) dc[,(j+1)] <- p[,1]
      }
      
      rc@response[[n]] <- dc
      colnames(rc@response[[n]]) <- c(n,na)
    }
  }
  
  rc
}

# .responseCurve <- function(m,id,si=100,includeTest=FALSE,.fun=base::mean) {
#   if (missing(.fun) || !is.function(.fun)) .fun <- base::mean
#   
#   mi <- m@run.info[m@run.info$modelID %in% id,]
#   
#   mi <- mi[mi$success,]
#   
#   sp <- unique(as.character(mi$species))
#   
#   if (length(sp) > 1) {
#     warning('More than one species exists in the model object; only first species is considered!')
#     sp <- sp[1]
#     mi <- mi[mi$species == sp,]
#     id <-mi$modelID
#   }
#   
#   na <- paste0(mi$method,'_ID-',mi$modelID)
#   #--------
#   rc <- new('.responseCurve')
#   
#   rc@multi <- nrow(mi) > 1
#   
#   
#   
#   if (includeTest) dt <- as.data.frame(m@data,sp=sp)
#   else dt <- as.data.frame(m@data,grp=c('train'),sp=sp)
#   
#   .wPr <- which(dt[,sp] == 1)
#   
#   nf <- m@setting@featureFrame@predictors
#   
#   rc@variables <- nf
#   
#   dt <- dt[,nf,drop=FALSE]
#   
#   nFact <- names(m@setting@featureFrame@categorical)
#   
#   if (any(nFact %in% nf)) {
#     nFact <- nFact[nFact %in% nf]
#     nf <- .excludeVector(nf,nFact)
#   } else nFact <- NULL
#   
#   rc@categorical <- nFact
#   
#   if (!is.null(nFact)) {
#     dm <- data.frame(matrix(nrow=1,ncol=length(nf)+length(nFact)))
#     colnames(dm) <- c(nf,nFact)
#     
#     for (n in nf) dm[,n] <- .fun(dt[.wPr,n],na.rm=TRUE)
#     
#     for (n in nFact) {
#       u <- .factorFreq(dt[,n])
#       dm[,n] <- names(which.max(u)[1])
#       dm[,n] <- factor(dm[,n],levels = levels(dt[,n]))
#     }
#     
#     for (n in nFact) {
#       #dt[,n] <- factor(dt[,n])
#       u <- unique(dt[,n])
#       dc <- data.frame(matrix(nrow=length(u),ncol=1+nrow(mi))) # the results of response curve for categorical variable
#       dc[,1] <- u
#       colnames(dc) <- c(n,na)
#       for (i in seq_along(u)) {
#         dv <- dm
#         dv[,n] <- u[i]
#         p <- predict(m,newdata=dv,id=mi$modelID)
#         dc[i,1] <- u[i]
#         for (j in 2:ncol(dc)) dc[i,j] <- p[1,(j-1)]
#       }
#       
#       rc@response[[n]] <- dc
#     }
#     #------
#     
#     
#     if (!is.null(nf)) {
#       dc <- dm
#       dm <- data.frame(matrix(nrow=si,ncol=length(nf)+length(nFact)))
#       colnames(dm) <- c(nf,nFact)
#       
#       for (n in c(nf,nFact)) dm[,n] <- dc[,n]
#       
#       for (n in nf) {
#         dv <- dm
#         dv[,n] <- seq(min(dt[,n],na.rm=TRUE),max(dt[,n],na.rm=TRUE),length.out = si)
#         dc <- data.frame(matrix(nrow=si,ncol=1+nrow(mi)))
#         dc[,1] <- dv[,n]
#         p <- predict(m,newdata=dv,id=mi$modelID)
#         for (j in 2:ncol(dc)) dc[,j] <- p[,(j-1)]
#         
#         rc@response[[n]] <- dc
#         colnames(rc@response[[n]]) <- c(n,na)
#       }
#     }
#   } else {
#     dm <- data.frame(matrix(nrow=si,ncol=length(nf)))
#     colnames(dm) <- nf
#     
#     for (n in nf) dm[,n] <- .fun(dt[,n],na.rm=TRUE)
#     
#     for (n in nf) {
#       dv <- dm
#       dv[,n] <- seq(min(dt[,n],na.rm=TRUE),max(dt[,n],na.rm=TRUE),length.out = si)
#       dc <- data.frame(matrix(nrow=si,ncol=1+nrow(mi)))
#       dc[,1] <- dv[,n]
#       for (j in seq_along(mi$modelID)) {
#         p <- predict(m,newdata=dv,id=mi$modelID[j])
#         dc[,(j+1)] <- p[,1]
#       }
#       
#       rc@response[[n]] <- dc
#       colnames(rc@response[[n]]) <- c(n,na)
#     }
#   }
#   
#   rc
# }
#-----------


if (!isGeneric("getResponseCurve")) {
  setGeneric("getResponseCurve", function(x,id,...)
    standardGeneric("getResponseCurve"))
}


setMethod("getResponseCurve", signature(x='sdmModels'),
          function(x,id,size=100,includeTest=FALSE,fun='mean',setting,...) {
            if (missing(fun)) fun <- base::mean
            else {
              if (is.character(fun)) {
                fun <- fun[1]
                if (fun == 'mean') fun <- base::mean
                else if (fun == 'median') fun <- median
                else if (fun == 'max') fun <- max
                else if (fun == 'min') fun <- min
                else {
                  warning('fun should be one of c("mean","median","max","min") or a function (default="mean" is considered)!')
                  fun <- base::mean
                }
              } else if (!is.function(fun)) {
                warning('fun should be one of c("mean","median","max","min") or a function (default="mean" is considered)!')
                fun <- base::mean
              }
            }
            
            if (missing(size)) size <- 100
            
            if (missing(includeTest)) includeTest <- FALSE
            
            
            if (missing(id)) {
              id <- x@run.info[which(x@run.info$success),1]
              if (length(id) == 1) cat(paste0('The id argument is not specified; id = ',id,' is considered. \n'))
              else if (length(id) > 1) cat(paste0('The id argument is not specified; The modelIDs of ', length(id),' successfully fitted models are assigned to id...! \n'))
              else stop('It seems that no successfully fitted models do exist in the sdmModels object!')
            } 
            
            if (is.numeric(id)) {
              .responseCurve(x,id=id,si=size,includeTest=includeTest,.fun=fun)
            } else if (is.character(id)) {
              if (!tolower(id[1]) %in% c('ens','ensemble','ensmble','ensmbl','en','ensembl','e')) stop('id should be either "ensemble" (character) or a numeric vector specifying model IDs!')
              id <- 'ensemble'
              .responseCurveEns(x,id=id,si=size,includeTest=includeTest,.fun=fun,setting = setting)
            }
          }
)


if (!isGeneric("rcurve")) {
  setGeneric("rcurve", function(x,n,id,mean,fun,confidence,gg,...)
    standardGeneric("rcurve"))
}  


setMethod("rcurve", signature(x='sdmModels'),
          function(x,n,id,mean,fun,confidence,gg,size,includeTest,setting,...) {
            
            if (missing(setting)) setting <- NULL
            
            if (missing(mean)) mean <- TRUE
            
            if (missing(fun)) fun <- base::mean
            else {
              if (is.character(fun)) {
                fun <- fun[1]
                if (fun == 'mean') fun <- base::mean
                else if (fun == 'median') fun <- median
                else if (fun == 'max') fun <- max
                else if (fun == 'min') fun <- min
                else {
                  warning('fun should be one of c("mean","median","max","min") or a function (default="mean" is considered)!')
                  fun <- base::mean
                }
              } else if (!is.function(fun)) {
                warning('fun should be one of c("mean","median","max","min") or a function (default="mean" is considered)!')
                fun <- base::mean
              }
            }
           
            if (missing(size)) size <- 100
            
            if (missing(includeTest)) includeTest <- FALSE
            
            if (missing(n)) n <- NULL
            
            if (missing(gg)) gg <- .require('ggplot2')
            else if (gg && !.require('ggplot2')) gg <- FALSE
            
            if (missing(id)) {
              id <- x@run.info[which(x@run.info$success),1]
              if (length(id) == 1) cat(paste0('The id argument is not specified; id = ',id,' is considered. \n'))
              else if (length(id) > 1) cat(paste0('The id argument is not specified; The modelIDs of ', length(id),' successfully fitted models are assigned to id...! \n'))
              else stop('It seems that no (successfully fitted) models do exist in the sdmModels object!')
            }
            #################
            if (is.numeric(id)) {
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
              #--------
              rc <- .responseCurve(x,id=id,si=size,includeTest=includeTest,.fun=fun)
            } else if (is.character(id)) {
              if (!tolower(id[1]) %in% c('ens','ensemble','ensmble','ensmbl','en','ensembl','e')) stop('id should be either "ensemble" (character) or a numeric vector specifying model IDs!')
              confidence <- FALSE
              id <- 'ensemble'
              rc <- .responseCurveEns(x,id=id,si=size,includeTest=includeTest,.fun=fun,setting = setting)
            }
            
            plot(rc,y=n,gg=gg,mean=mean,confidence=confidence,...)
          }
)



setMethod("rcurve", signature(x='.responseCurve'),
          function(x,n,id,mean=TRUE,fun,confidence=TRUE,gg=TRUE,...) {
            if (missing(fun)) fun <- NULL
            
            if (missing(gg)) gg <- .require('ggplot2')
            else if (gg && !.require('ggplot2')) gg <- FALSE
            
            if (missing(mean)) mean <- TRUE
            
            if (missing(confidence)) confidence <- TRUE
            
            if (missing(n)) n <- NULL
            
            plot(x,y=n,gg=gg,mean=mean,confidence=confidence,...)
          }
)
