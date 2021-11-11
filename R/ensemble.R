# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Oct. 2016
# Last Update :  Nov, 2021
# Version 3.0
# Licence GPL v3




.ent <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  n1 <- length(which(x == 1))
  n0 <- n - n1
  x <- c(n0,n1)/n
  if (any(x == 0)) {
    return (0)
  } else return(- sum(x*log(x)) / log(2))
}


#--------

.getWeights <- function(x, mi, wtest, id, stat, opt) {
  if (all(mi[,wtest[1]])) {
    weight <- getEvaluation(x,w = id,wtest = wtest[1],stat = stat,opt = opt)[,2]
  } else {
    id1 <- mi$modelID[mi[,wtest[1]]]
    if (length(wtest) > 1) {
      w <- which(!mi[,wtest[1]])
      if (any(mi[w,wtest[2]])) {
        id2 <- mi$modelID[w][mi[w,wtest[2]]]
      } else id2 <- NULL
      if (length(c(id1,id2)) < length(id)) {
        if (length(wtest) == 3) {
          w <- id[which(!id %in% c(id1,id2))]
          w <- which(mi$modelID %in% w)
          if (any(mi[w,wtest[3]])) {
            id3 <- mi$modelID[w][mi[w,wtest[3]]]
          } else id3 <- NULL
        } else id3 <- NULL
      } else id3 <- NULL
    } else id2 <- id3 <- NULL
    
    id <- id[which(id %in% c(id1,id2,id3))]
    mi <- mi[mi$modelID %in% id,]
    id <- mi$modelID
    w1 <- getEvaluation(x,w = id1,wtest = wtest[1],stat = stat,opt = opt)[,2]
    if (!is.null(id2)) w2 <- getEvaluation(x,w = id2,wtest = wtest[2],stat = stat,opt = opt)[,2]
    else w2 <- NULL
    if (!is.null(id3)) w3 <- getEvaluation(x,w = id3,wtest = wtest[2],stat = stat,opt = opt)[,2]
    else w3 <- NULL
    idw <- c(id1,id2,id3)
    w <- c(w1,w2,w3)
    weight <- w[sapply(id,function(x) which(idw == x))]
  }
  
  list(weight=weight,mi=mi)
  
}
#--------
if (!isGeneric("ensemble")) {
  setGeneric("ensemble", function(x,newdata,filename,setting,...)
    standardGeneric("ensemble"))
}	




setMethod('ensemble', signature(x='sdmModels',newdata='Raster'), 
          function(x, newdata, filename,setting,obj.size=1L,overwrite=FALSE,...) {
            if (missing(setting)) setting <- list(method='weighted',stat='AUC',wtest=NULL,power=1)
            
            if (missing(filename)) filename <- ''
            
            if (missing(overwrite)) overwrite <- FALSE
            
            if (filename != '' && !overwrite) {
              if (file.exists(filename)) stop('The specified filename does exist. You may use overwrite = TRUE or use a different filename...!')
            }
            
            if (missing(obj.size)) obj.size <- 1L
            
            if (filename == '') filename <- .generateName('sdm_ensemble')
            if (extension(filename) == '') filename <- paste(filename,'.grd',sep='')
            
            if (!is.list(setting)) stop('setting should be defined as a list object!')
            
            n <- names(setting)
            if (is.null(n) || '' %in% n) stop('one or more arguments in "setting" do not have names!')
            
            n <- .pmatch(n,c('method','stat','opt','id','wtest','weight','power'))
            if (any(is.na(n))) stop('At least one argument in setting is unkown!')
            #------------
            mi <- getModelInfo(x)
            mi <- mi[mi$success,]
            if (nrow(mi) == 0) stop('There is no model in the input object!')
            if (nrow(mi) < 2) stop('The input object has only 1 model, while at least 2 models are needed for ensemble!')
            #-----------id:
            if ('id' %in% n) {
              id <- setting[[which(n == 'id')]]
              if (!is.numeric(id)) {
                warning('id should be numeric. It is ignored, and all model IDs are considered!')
                id <- mi$modelID
              } 
            } else id <- mi$modelID
            
            
            
            w <- which(mi$modelID %in% id)
            if (length(w) < length(id)) {
              if (length(w) < 2) stop('At least 2 models should be selected for ensemble; change "id" numbers!')
              warning(paste(length(id) - length(w),'of the specified IDs do not exist in the input object, and are excluded from id'))
            }
            
            id <- mi$modelID[mi$modelID %in% id]
            mi <- mi[w,]
            #-------- method:
            if ('method' %in% n) {
              method <- tolower(setting[[which(n == 'method')]])
              if (method %in% c('weight','weighted','weited','w','wei','weig','weigh')) method <- 'weighted'
              else if (method %in% c('unweight','unweighted','unweited','unw','unwei','unweig','unweigh','u','uw')) method <- 'unweighted'
              else if (method %in% c('median','med','medi','mdian')) method <- 'median'
              else if (method %in% c('pamean','pa','pa_mean','pa-mean')) method <- 'pa'
              else if (method %in% c('mean-weighted','mean_weighted','m-w','mean-w','u-w','unweighted-weighted','unweighted_weighted')) method <- 'mean-weighted'
              else if (method %in% c('mean-unweighted','mean_unweighted','m-u','mean-u','u-u','unweighted-unweighted','unweighted_unweighted','mean-mean')) method <- 'mean-unweighted'
              else if (method %in% c('median-weighted','median_weighted','med-w','median-w','u-w')) method <- 'median-weighted'
              else if (method %in% c('median-unweighted','median_unweighted','med-u','median-u','median-mean')) method <- 'median-unweighted'
              else if (method %in% c('uncertainty','entropy','uncert','uncertain')) method <- 'entropy'
              else {
                warning('method is not identified, therefore the default ("weighted") is considered!')
                method <- 'weighted'
              }
            } else method <- 'weighted'
            #--------- stat:
            if ('stat' %in% n) {
              s <- c('AUC','COR','Deviance','sensitivity','specificity','TSS','Kappa','NMI','phi','ppv','npv','ccr')
              stat <- tolower(setting[[which(n == 'stat')]])
              if (length(stat) > 1) {
                warning('only one statistic should be used. The first statistic in stat is used')
                stat <- stat[1]
              }
              stat <- .pmatch(stat,s)
              if (is.na(stat)) stop('stat is unknown!')
            } else stat <- 'AUC'
            #--------- wtest:
            if ('wtest' %in% n) {
              wtest <- tolower(setting[[which(n == 'wtest')]])
              if (!is.null(wtest)) wtest <- .pmatch(wtest,c('training','test.dep','test.indep'))[1]
              if (is.null(wtest) && is.na(wtest)) {
                wtest <- colnames(mi)[9:7][which(c(length(which(mi[,9])),length(which(mi[,8])),length(which(mi[,7]))) > 0)]
              }
            } else wtest <- colnames(mi)[9:7][which(c(length(which(mi[,9])),length(which(mi[,8])),length(which(mi[,7]))) > 0)]
            #---------- opt
            if ('opt' %in% n) {
              opt <- setting[[which(n == 'opt')]][1]
              if (!is.null(opt)) {
                if (is.numeric(opt)) {
                  if (!opt %in% 1:10) {
                    opt <- 2
                    warning('opt (the criteria for optimum threshold) should be a number between 1:10; opt=2 is considered (i.e., max(se+sp))')
                  }
                } else {
                  th.criteria <- c("sp=se","max(se+sp)","min(cost)","minROCdist","max(kappa)","max(ppv+npv)","ppv=npv","max(NMI)","max(ccr)","prevalence")
                  opt <- .pmatch(opt,th.criteria)
                  opt <- opt[!is.na(opt)]
                  if (length(opt) == 0) {
                    warning('opt (the criteria for optimum threshold) is not understood! max(se+sp) is considered...')
                    opt <- 2
                  } else {
                    opt <- which(th.criteria == opt)
                  }
                } 
              } else opt <- 2
            } else opt <- 2
            #---------- weight
            if ('weight' %in% n) {
              weight <- setting[[which(n == 'weight')]]
              if (!is.numeric(weight)) stop('weight should be numeric!')
              if (length(weight) != length(id)) stop('the weight has a different length than the number of models in the input object (or the length of id)!')
            } else weight <- NULL
            
            
            #--------------
            if ('power' %in% n) {
              .power <- setting[[which(n == 'power')]]
              if (!is.numeric(.power) && method %in% c('weighted','mean-weighted','median-weighted')) {
                warning('power should be numeric. It is ignored, and the default (=1) is considered!')
                .power <- 1
              }
            } else .power <- 1
            #---------
            
            .n <- names(newdata)
            p <- NULL
            if (!all(x@setting@featuresFrame@vars %in% .n)) {
              if (length(.n) == length(id)) {
                cat('\n ......... the Raster object is used as the predicted probabilities...\n')
                p <- newdata
              } else stop('the data does not contain some or all of the variables that the model needs...')
            }
            
            nf <- nFact <- NULL
            
            nf <- .getFeatureNamesTypes(x@setting@featuresFrame)
            if (!is.null(nf) && 'factor' %in% nf[,2]) {
              nFact <- as.character(nf[nf[,2] == 'factor',1])
              nf <- nf$name
              nf <- .excludeVector(nf,nFact)
            } else nf <- nf$name
            
            
            .nd <- .raster2df(newdata,.getlevels(x@data)[nFact])
            w <- which(colnames(.nd) == "cellnr")
            .cellnr <- .nd$cellnr
            .nd <- .nd[,-w]
            
            if (!is.null(p)) p <- .nd
            
            .nr <- nrow(.nd) ## !
            #memreq <- (object.size(.nd)[[1]] / ncol(.nd))*.nr
            ####################
            memreq <- 1 # This needs to be updated to calculate the required memory, and if it is not fitted, .raster2df should be called for chunck of the raster data!
            ######################
            .ens <- raster(newdata)
            
            
            
            ####################
            if (method == 'weighted') {
              if (is.null(weight)) {
                weight <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = stat, opt = opt)
                mi <- weight$mi
                id <- mi$modelID
                weight <- weight$weight
              }
              f <- function(x) {
                sum(weight * x)
              }
              
              f2 <- function(x) {
                if (!all(is.na(x))) {
                  .w <- which(!is.na(x))
                  .we <- weight[.w]
                  .we <- .we / sum(.we)
                  sum(.we * x[.w])
                } else NA
              }
              
              weight <- weight ^ .power
              weight <- weight / sum(weight)
              
              if (is.null(p)) p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              else {
                if (length(id) != ncol(p)) p <- p[,id]
              }
              
              .p <- apply(p,1,f)
              .wn <- which(is.na(.p))
              if (length(.wn) > 1) {
                .p[.wn] <- apply(p[.wn,],1,f2)
              }
              
              rm(p); gc()
              
              .ens[.cellnr] <- .p
              
              # if (memreq <= (obj.size * 1073741824)) {
              #   
              #   
              # } else {
              #   memdiv <- ceiling(obj.size / (memreq /  .nr))
              #   ii <- ceiling(.nr/memdiv)
              #   
              #   for (i in 1:ii) {
              #     .ii <- (i-1) * memdiv + c(1:memdiv)
              #     .ii <- .ii[.ii %in% 1:.nr]
              #     p <- predict(x,newdata=.nd[.ii,],w=id,obj.size=obj.size,...)
              #     .p <- apply(p,1,f)
              #     .wn <- which(is.na(.p))
              #     if (length(.wn) > 1) {
              #       .p[.wn] <- apply(p[.wn,],1,f2)
              #     }
              #     
              #     rm(p); gc()
              #     
              #     .ens[.cellnr[.ii]] <- .p
              #   }
              # }
              
            } else if (method == 'unweighted') {
              
              if (is.null(p)) p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              .p <- apply(p,1,mean,na.rm=TRUE)
              
              rm(p); gc()
              
              .ens[.cellnr] <- .p
              
              # if (memreq <= (obj.size * 1073741824)) {
              #   p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              #   .p <- apply(p,1,mean,na.rm=TRUE)
              #   
              #   rm(p); gc()
              #   
              #   .ens[.cellnr] <- .p
              #   
              # } else {
              #   memdiv <- ceiling(obj.size / (memreq /  .nr))
              #   ii <- ceiling(.nr/memdiv)
              #   
              #   for (i in 1:ii) {
              #     .ii <- (i-1) * memdiv + c(1:memdiv)
              #     .ii <- .ii[.ii %in% 1:.nr]
              #     p <- predict(x,newdata=.nd[.ii,],w=id,obj.size=obj.size,...)
              #     .p <- apply(p,1,mean,na.rm=TRUE)
              #     rm(p); gc()
              #     .ens[.cellnr[.ii]] <- .p
              #   }
              # }
            }  else if (method == 'median') {
              if (is.null(p)) p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              .p <- apply(p,1,median,na.rm=TRUE)
              
              rm(p); gc()
              
              .ens[.cellnr] <- .p
            } else if (method == 'pa') {
              .th <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = 'threshold', opt = opt)
              mi <- weight$mi
              id <- mi$modelID
              .th <- .th$weight
              
              #-----
              if (is.null(p)) p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              else {
                if (length(id) < ncol(p)) p <- p[,id]
              }
              
              f <- function(x) {
                .x <- x >= .th
                x[] <- 0
                x[.x] <- 1
                sum(x) / length(x)
              }
              
              f2 <- function(x) {
                if (!all(is.na(x))) {
                  .w <- which(!is.na(x))
                  x <- x[.w]
                  .x <- x >= .th[.w]
                  x[] <- 0
                  x[.x] <- 1
                  sum(x) / length(x)
                } else NA
              }
              
              .p <- apply(p,1,f)
              .wn <- which(is.na(.p))
              if (length(.wn) > 1) {
                .p[.wn] <- apply(p[.wn,],1,f2)
              }
              
              rm(p); gc()
              
              .ens[.cellnr] <- .p
            } else if (method == 'mean-weighted') {
              
              if (is.null(weight)) {
                weight <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = stat, opt = opt)
                mi <- weight$mi
                id <- mi$modelID
                weight <- weight$weight
              }
              
              #weight <- weight / sum(weight)
              
              if (length(weight) != nrow(mi)) stop('something is wrong...! [length of weights is not equal to the number of models...]')
              
              if (is.null(p)) {
                
                p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,mean=TRUE,...)
                
                if (ncol(p) == 1) {
                  warning('It seems that only one method was used to fit the models, therefore, the predictions are combined using unweighted method...!')
                  .ens <- p[,1]
                } else {
                  weight2 <- c()
                  mi$method <- as.character(mi$method)
                  for (.u in unique(mi$method)) {
                    w <- which(mi$method == .u)
                    weight2 <- c(weight2,mean(weight[w],na.rm=TRUE))
                  }
                  weight2 <- weight2 ^ .power
                  weight <- weight2 / sum(weight2)
                  
                  if (ncol(p) != length(weight)) {
                    warning('"mean-unweighted" is used...!!! \n something is wrong...: length of weight vector is not equal to the number of modelling methods...!')
                    .p <- apply(p,1,mean,na.rm=TRUE)
                    
                    rm(p); gc()
                    
                    .ens[.cellnr] <- .p
                  } else {
                    
                    f <- function(x) {
                      sum(weight * x)
                    }
                    
                    f2 <- function(x) {
                      if (!all(is.na(x))) {
                        .w <- which(!is.na(x))
                        .we <- weight[.w]
                        .we <- .we / sum(.we)
                        sum(.we * x[.w])
                      } else NA
                    }
                    
                    .p <- apply(p,1,f)
                    .wn <- which(is.na(.p))
                    if (length(.wn) > 1) {
                      .p[.wn] <- apply(p[.wn,],1,f2)
                    }
                    
                    rm(p); gc()
                    
                    .ens[.cellnr] <- .p
                  }
                }
              } else {
                
                .me <- as.character(mi$method)
                
                if (length(unique(.me)) == 1) {
                  warning('It seems that only one method was used to fit the models, therefore, "unweighted" method is used...!')
                  .p <- apply(p,1,mean,na.rm=TRUE)
                } else {
                  .p <- data.frame(matrix(NA,nrow=nrow(p),ncol=length(unique(.me))))
                  colnames(.p) <- unique(.me)
                  
                  weight2 <- c()
                  
                  for (.u in unique(.me)) {
                    w <- which(mi$method == .u)
                    if (length(w) > 0) {
                      if (length(w) == 1) .p[,.u] <- p[,w]
                      else .p[,.u] <- apply(p[,w],1,mean,na.rm=TRUE)
                      
                      weight2 <- c(weight2,mean(weight[w],na.rm=TRUE))
                    }
                  }
                  weight2 <- weight2 ^ .power
                  weight <- weight2 / sum(weight2)
                  ######
                  if (ncol(.p) == 1) {
                    warning('It seems that only one method was used to fit the models, therefore, the predictions are combined using unweighted method...!')
                    .ens <- .p[,1]
                  } else {
                    if (ncol(.p) != length(weight)) {
                      warning('"mean-unweighted" is used...!!! \n something is wrong...: length of weight vector is not equal to the number of modelling methods...!')
                      .p <- apply(.p,1,mean,na.rm=TRUE)
                      
                      rm(p); gc()
                      
                      .ens[.cellnr] <- .p
                    } else {
                      
                      f <- function(x) {
                        sum(weight * x)
                      }
                      
                      f2 <- function(x) {
                        if (!all(is.na(x))) {
                          .w <- which(!is.na(x))
                          .we <- weight[.w]
                          .we <- .we / sum(.we)
                          sum(.we * x[.w])
                        } else NA
                      }
                      
                      .pp <- apply(.p,1,f)
                      .wn <- which(is.na(.pp))
                      if (length(.wn) > 1) {
                        .pp[.wn] <- apply(.p[.wn,],1,f2)
                      }
                      
                      rm(p); gc()
                      
                      .ens[.cellnr] <- .pp
                    }
                  }
                }
              }
            }  else if (method == 'mean-unweighted') {
              if (is.null(p)) {
                p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,mean=TRUE,...)
                .p <- apply(p,1,mean,na.rm=TRUE)
              } else {
                .me <- as.character(mi$method)
                
                if (length(unique(.me)) == 1) {
                  warning('It seems that only one method was used to fit the models, therefore, "unweighted" method is used...!')
                  .p <- apply(p,1,mean,na.rm=TRUE)
                } else {
                  .p <- data.frame(matrix(NA,nrow=nrow(p),ncol=length(unique(.me))))
                  colnames(.p) <- unique(.me)
                  
                  for (.u in unique(.me)) {
                    w <- which(mi$method == .u)
                    if (length(w) > 0) {
                      if (length(w) == 1) .p[,.u] <- p[,w]
                      else .p[,.u] <- apply(p[,w],1,mean,na.rm=TRUE)
                    }
                  }
                  .p <- apply(.p,1,mean,na.rm=TRUE)
                }
                
              }
              
              rm(p); gc()
              
              .ens[.cellnr] <- .p
            } else if (method == 'median-weighted') {
              if (is.null(weight)) {
                weight <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = stat, opt = opt)
                mi <- weight$mi
                id <- mi$modelID
                weight <- weight$weight
              }
              
              if (length(weight) != nrow(mi)) stop('something is wrong...! [length of weights is not equal to the number of models...]')
              
              if (is.null(p)) {
                p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              } else {
                if (length(id) < ncol(p)) p <- p[,id]
              }
              
              .me <- as.character(mi$method)
              .p <- data.frame(matrix(NA,nrow=nrow(p),ncol=length(unique(.me))))
              colnames(.p) <- unique(.me)
              
              if (ncol(.p) == 1) {
                warning('It seems that only one method was used to fit the models, therefore, the predictions are combined using median method...!')
                .ens <- apply(p,1,median,na.rm=TRUE)
              } else {
                weight2 <- c()
                
                for (.u in unique(.me)) {
                  w <- grep(.u,colnames(p))
                  if (length(w) > 0) {
                    weight2 <- c(weight2,mean(weight[which(mi$method == .u)],na.rm=TRUE))
                    
                    if (length(w) == 1) .p[,.u] <- p[,w]
                    else .p[,.u] <- apply(p[,w],1,median,na.rm=TRUE)
                  }
                }
                #------
                rm(p); gc()
                
                weight2 <- weight2 ^ .power
                weight <- weight2 / sum(weight2)
                
                if (ncol(.p) != length(weight)) {
                  warning('"median-unweighted" is used...!!! \n something is wrong...: length of weight vector is not equal to the number of modelling methods...!')
                  .p <- apply(.p,1,mean,na.rm=TRUE)
                  .ens <- .p
                } else {
                  
                  f <- function(x) {
                    sum(weight * x)
                  }
                  
                  f2 <- function(x) {
                    if (!all(is.na(x))) {
                      .w <- which(!is.na(x))
                      .we <- weight[.w]
                      .we <- .we / sum(.we)
                      sum(.we * x[.w])
                    } else NA
                  }
                  #---------
                  .p2 <- apply(.p,1,f)
                  .wn <- which(is.na(.p2))
                  if (length(.wn) > 1) {
                    .p2[.wn] <- apply(.p[.wn,],1,f2)
                  }
                  
                  rm(.p); gc()
                  
                  .ens[.cellnr] <- .p2
                }
              }
            } else if (method == 'median-unweighted') {
              if (is.null(p)) {
                p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
                
                .me <- as.character(mi$method)
                
                if (length(unique(.me)) == 1) {
                  warning('It seems that only one method was used to fit the models, therefore, the predictions are combined using median method...!')
                  .ens <- apply(p,1,median,na.rm=TRUE)
                } else {
                  .p <- data.frame(matrix(NA,nrow=nrow(p),ncol=length(unique(.me))))
                  colnames(.p) <- unique(.me)
                  
                  for (.u in unique(.me)) {
                    w <- grep(.u,colnames(p))
                    if (length(w) > 0) {
                      if (length(w) == 1) .p[,.u] <- p[,w]
                      else .p[,.u] <- apply(p[,w],1,median,na.rm=TRUE)
                    }
                  }
                  #------
                }
              
              } else {
                .me <- as.character(mi$method)
                
                if (length(unique(.me)) == 1) {
                  warning('It seems that only one method was used to fit the models, therefore, the predictions are combined using median method...!')
                  .ens <- apply(p,1,median,na.rm=TRUE)
                } else {
                  .p <- data.frame(matrix(NA,nrow=nrow(p),ncol=length(unique(.me))))
                  colnames(.p) <- unique(.me)
                  
                  for (.u in unique(.me)) {
                    w <- which(mi$method == .u)
                    if (length(w) > 0) {
                      if (length(w) == 1) .p[,.u] <- p[,w]
                      else .p[,.u] <- apply(p[,w],1,median,na.rm=TRUE)
                    }
                  }
                }
              }
              rm(p); gc()
              .ens[.cellnr] <- apply(.p,1,mean,na.rm=TRUE)
            } else if (method == 'entropy') {
              .th <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = 'threshold', opt = opt)
              mi <- weight$mi
              id <- mi$modelID
              .th <- .th$weight
              #-----
              if (is.null(p)) p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              else {
                if (length(id) != ncol(p)) p <- p[,id]
              }
              
              f <- function(x) {
                .x <- x >= .th
                x[] <- 0
                x[.x] <- 1
                .ent(x)
              }
              
              .p <- apply(p,1,f)
              
              rm(p); gc()
              
              .ens[.cellnr] <- .p
              
            }
            
            
            if (filename != '') .ens <- writeRaster(.ens,filename,overwrite = overwrite)
            
            .ens
          }
)

#-------


setMethod('ensemble', signature(x='sdmModels',newdata='data.frame'), 
          function(x, newdata, filename,setting,obj.size=1L,overwrite=FALSE,...) {
            if (missing(setting)) setting <- list(method='weighted',stat='AUC',wtest=NULL,power=1)
            if (missing(filename)) filename=''
            
            if (missing(overwrite)) overwrite <- FALSE
            
            if (filename != '' && !overwrite) {
              if (file.exists(filename)) stop('The specified filename does exist. You may use overwrite = TRUE or use a different filename...!')
            } 
            
            if (missing(obj.size)) obj.size <- 1L
            
            if (filename == '') filename <- .generateName('sdm_ensemble')
            if (extension(filename) == '') filename <- paste(filename,'.grd',sep='')
            
            if (!is.list(setting)) stop('setting should be defined as a list object!')
            
            n <- names(setting)
            if (is.null(n) || '' %in% n) stop('one or more arguments in setting do not have a name!')
            
            n <- .pmatch(n,c('method','stat','opt','id','wtest','weight','power'))
            if (any(is.na(n))) stop('At least one argument in setting is unkown!')
            #------------
            mi <- getModelInfo(x)
            mi <- mi[mi$success,]
            if (nrow(mi) == 0) stop('There is no model in the input object!')
            if (nrow(mi) < 2) stop('The input object has only 1 model, while at least 2 models are needed for ensemble!')
            #-----------id:
            if ('id' %in% n) {
              id <- setting[[which(n == 'id')]]
              if (!is.numeric(id)) {
                warning('id should be numeric. It is ignored, and all model IDs are considered!')
                id <- mi$modelID
              } 
            } else id <- mi$modelID
            
            w <- which(mi$modelID %in% id)
            if (length(w) < length(id)) {
              if (length(w) < 2) stop('At least 2 models should be selected for ensemble; change "id" numbers!')
              warning(paste(length(id) - length(w),'of the specified IDs do not exist in the input object, and are excluded from id'))
            }
            
            id <- mi$modelID[mi$modelID %in% id]
            mi <- mi[w,]
            #-------- method:
            if ('method' %in% n) {
              method <- tolower(setting[[which(n == 'method')]])
              if (method %in% c('weight','weighted','weited','w','wei','weig','weigh')) method <- 'weighted'
              else if (method %in% c('unweight','unweighted','unweited','unw','unwei','unweig','unweigh','u','uw','mean')) method <- 'unweighted'
              else if (method %in% c('median','med','medi','mdian')) method <- 'median'
              else if (method %in% c('pamean','pa','pa_mean','pa-mean')) method <- 'pa'
              else if (method %in% c('mean-weighted','mean_weighted','m-w','mean-w','u-w','unweighted-weighted','unweighted_weighted')) method <- 'mean-weighted'
              else if (method %in% c('mean-unweighted','mean_unweighted','m-u','mean-u','u-u','unweighted-unweighted','unweighted_unweighted','mean-mean')) method <- 'mean-unweighted'
              else if (method %in% c('median-weighted','median_weighted','med-w','median-w','u-w')) method <- 'mean-weighted'
              else if (method %in% c('median-unweighted','median_unweighted','med-u','median-u','median-mean')) method <- 'median-unweighted'
              else if (method %in% c('uncertainty','entropy','uncert','uncertain')) method <- 'entropy'
              else {
                warning('method is not identified, therefore the default ("weighted") is considered!')
                method <- 'weighted'
              }
            } else method <- 'weighted'
            #--------- stat:
            if ('stat' %in% n) {
              s <- c('AUC','COR','Deviance','sensitivity','specificity','TSS','Kappa','NMI','phi','ppv','npv','ccr')
              stat <- tolower(setting[[which(n == 'stat')]])
              if (length(stat) > 1) {
                warning('only one statistic should be used. The first statistic in stat is used')
                stat <- stat[1]
              }
              stat <- .pmatch(stat,s)
              if (is.na(stat)) stop('stat is unknown!')
            } else stat <- 'AUC'
            #--------- wtest:
            if ('wtest' %in% n) {
              wtest <- tolower(setting[[which(n == 'wtest')]])
              if (!is.null(wtest)) wtest <- .pmatch(wtest,c('training','test.dep','test.indep'))[1]
              if (is.null(wtest) && is.na(wtest)) {
                wtest <- colnames(mi)[9:7][which(c(length(which(mi[,9])),length(which(mi[,8])),length(which(mi[,7]))) > 0)]
              }
            } else wtest <- colnames(mi)[9:7][which(c(length(which(mi[,9])),length(which(mi[,8])),length(which(mi[,7]))) > 0)]
            #---------- opt
            if ('opt' %in% n) {
              opt <- setting[[which(n == 'opt')]][1]
              if (!is.null(opt)) {
                if (is.numeric(opt)) {
                  if (!opt %in% 1:10) {
                    opt <- 2
                    warning('opt (the criteria for optimum threshold) should be a number between 1:10; opt=2 is considered (i.e., max(se+sp))')
                  }
                } else {
                  th.criteria <- c("sp=se","max(se+sp)","min(cost)","minROCdist","max(kappa)","max(ppv+npv)","ppv=npv","max(NMI)","max(ccr)","prevalence")
                  opt <- .pmatch(opt,th.criteria)
                  opt <- opt[!is.na(opt)]
                  if (length(opt) == 0) {
                    warning('opt (the criteria for optimum threshold) is not understood! max(se+sp) is considered...')
                    opt <- 2
                  } else {
                    opt <- which(th.criteria == opt)
                  }
                } 
              } else opt <- 2
            } else opt <- 2
            #---------- weight
            if ('weight' %in% n) {
              weight <- setting[[which(n == 'weight')]]
              if (!is.numeric(weight)) stop('weight should be numeric!')
              if (length(weight) != length(id)) stop('the weight has a different length than the number of models in the input object (or the length of id)!')
            } else weight <- NULL
            #--------------
            if ('power' %in% n) {
              .power <- setting[[which(n == 'power')]]
              if (!is.numeric(.power) && method %in% c('weighted','mean-weighted','median-weighted')) {
                warning('power should be numeric. It is ignored, and the default (=1) is considered!')
                .power <- 1
              }
            } else .power <- 1
            #---------
            
            .n <- colnames(newdata)
            
            if (!all(x@setting@featuresFrame@vars %in% .n)) stop('the data does not contain some or all of the variables that the model needs...')
            
            nf <- nFact <- NULL
            
            nf <- .getFeatureNamesTypes(x@setting@featuresFrame)
            if (!is.null(nf) && 'factor' %in% nf[,2]) {
              nFact <- as.character(nf[nf[,2] == 'factor',1])
              nf <- nf$name
              nf <- .excludeVector(nf,nFact)
            } else nf <- nf$name
            
            
            .nd <- newdata[,c(nf,nFact)]
            .nr <- nrow(.nd) ## !
            #memreq <- (object.size(.nd)[[1]] / ncol(.nd))*.nr
            
            ####################
            memreq <- 1 # This needs to be updated to calculate the required memory
            ######################
            
            .ens <- rep(NA,nrow(.nd))
            
            
            
            ####################
            if (method == 'weighted') {
              if (is.null(weight)) {
                weight <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = stat, opt = opt)
                mi <- weight$mi
                id <- mi$modelID
                weight <- weight$weight
              }
              f <- function(x) {
                sum(weight * x)
              }
              
              f2 <- function(x) {
                if (!all(is.na(x))) {
                  .w <- which(!is.na(x))
                  .we <- weight[.w]
                  .we <- .we / sum(.we)
                  sum(.we * x[.w])
                } else NA
              }
              weight <- weight ^ .power
              weight <- weight / sum(weight)
              
              if (memreq <= (obj.size * 1073741824)) {
                p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
                .p <- apply(p,1,f)
                .wn <- which(is.na(.p))
                if (length(.wn) > 1) {
                  .p[.wn] <- apply(p[.wn,],1,f2)
                }
                
                rm(p); gc()
                
                .ens <- .p
                
              } else {
                memdiv <- ceiling(obj.size / (memreq /  .nr))
                ii <- ceiling(.nr/memdiv)
                
                for (i in 1:ii) {
                  .ii <- (i-1) * memdiv + c(1:memdiv)
                  .ii <- .ii[.ii %in% 1:.nr]
                  p <- predict(x,newdata=.nd[.ii,],w=id,obj.size=obj.size,...)
                  .p <- apply(p,1,f)
                  .wn <- which(is.na(.p))
                  if (length(.wn) > 1) {
                    .p[.wn] <- apply(p[.wn,],1,f2)
                  }
                  
                  rm(p); gc()
                  
                  .ens[.ii] <- .p
                }
              }
              
            } else if (method == 'unweighted') {
              
              if (memreq <= (obj.size * 1073741824)) {
                p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
                .p <- apply(p,1,mean,na.rm=TRUE)
                
                rm(p); gc()
                
                .ens <- .p
                
              } else {
                memdiv <- ceiling(obj.size / (memreq /  .nr))
                ii <- ceiling(.nr/memdiv)
                
                for (i in 1:ii) {
                  .ii <- (i-1) * memdiv + c(1:memdiv)
                  .ii <- .ii[.ii %in% 1:.nr]
                  p <- predict(x,newdata=.nd[.ii,],w=id,obj.size=obj.size,...)
                  .p <- apply(p,1,mean,na.rm=TRUE)
                  rm(p); gc()
                  .ens[.ii] <- .p
                }
              }
            } else if (method == 'median') {
              p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              .p <- apply(p,1,median,na.rm=TRUE)
              
              rm(p); gc()
              
              .ens <- .p
            } else if (method == 'pa') {
              .th <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = 'threshold', opt = opt)
              mi <- weight$mi
              id <- mi$modelID
              .th <- .th$weight
              #-----
              p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              f <- function(x) {
                .x <- x >= .th
                x[] <- 0
                x[.x] <- 1
                sum(x) / length(x)
              }
              
              f2 <- function(x) {
                if (!all(is.na(x))) {
                  .w <- which(!is.na(x))
                  x <- x[.w]
                  .x <- x >= .th[.w]
                  x[] <- 0
                  x[.x] <- 1
                  sum(x) / length(x)
                } else NA
              }
              
              .p <- apply(p,1,f)
              .wn <- which(is.na(.p))
              if (length(.wn) > 1) {
                .p[.wn] <- apply(p[.wn,],1,f2)
              }
              
              rm(p); gc()
              
              .ens <- .p
            } else if (method == 'mean-weighted') {
              if (is.null(weight)) {
                weight <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = stat, opt = opt)
                mi <- weight$mi
                id <- mi$modelID
                weight <- weight$weight
              }
              
              if (length(weight) != nrow(mi)) stop('something is wrong...! [length of weights is not equal to the number of models...]')
              
              p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,mean=TRUE,...)
              if (ncol(p) == 1) {
                warning('It seems that only one method was used to fit the models, therefore, the predictions are combined using unweighted method...!')
                .ens <- p[,1]
              } else {
                weight2 <- c()
                mi$method <- as.character(mi$method)
                for (.u in unique(mi$method)) {
                  w <- which(mi$method == .u)
                  weight2 <- c(weight2,mean(weight[w],na.rm=TRUE))
                }
                weight2 <- weight2 ^ .power
                weight <- weight2 / sum(weight2)
                
                if (ncol(p) != length(weight)) {
                  warning('"mean-unweighted" is used...!!! \n something is wrong...: length of weight vector is not equal to the number of modelling methods...!')
                  .p <- apply(p,1,mean,na.rm=TRUE)
                  
                  rm(p); gc()
                  
                  .ens <- .p
                } else {
                  
                  f <- function(x) {
                    sum(weight * x)
                  }
                  
                  f2 <- function(x) {
                    if (!all(is.na(x))) {
                      .w <- which(!is.na(x))
                      .we <- weight[.w]
                      .we <- .we / sum(.we)
                      sum(.we * x[.w])
                    } else NA
                  }
                  
                  .p <- apply(p,1,f)
                  .wn <- which(is.na(.p))
                  if (length(.wn) > 1) {
                    .p[.wn] <- apply(p[.wn,],1,f2)
                  }
                  
                  rm(p); gc()
                  
                  .ens <- .p
                }
              }
              
            } else if (method == 'mean-unweighted') {
              p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,mean=TRUE,...)
              .p <- apply(p,1,mean,na.rm=TRUE)
              
              rm(p); gc()
              
              .ens <- .p
            } else if (method == 'median-weighted') {
              if (is.null(weight)) {
                weight <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = stat, opt = opt)
                mi <- weight$mi
                id <- mi$modelID
                weight <- weight$weight
              }
              
              if (length(weight) != nrow(mi)) stop('something is wrong...! [length of weights is not equal to the number of models...]')
              
              p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              .me <- as.character(mi$method)
              .p <- data.frame(matrix(NA,nrow=nrow(p),ncol=length(unique(.me))))
              colnames(.p) <- unique(.me)
              
              if (ncol(.p) == 1) {
                warning('It seems that only one method was used to fit the models, therefore, the predictions are combined using median method...!')
                .ens <- apply(p,1,median,na.rm=TRUE)
              } else {
                weight2 <- c()
                
                for (.u in unique(.me)) {
                  w <- grep(.u,colnames(p))
                  if (length(w) > 0) {
                    weight2 <- c(weight2,mean(weight[which(mi$method == .u)],na.rm=TRUE))
                    
                    if (length(w) == 1) .p[,.u] <- p[,w]
                    else .p[,.u] <- apply(p[,w],1,median,na.rm=TRUE)
                  }
                }
                #------
                rm(p); gc()
                
                weight2 <- weight2 ^ .power
                weight <- weight2 / sum(weight2)
                
                if (ncol(.p) != length(weight)) {
                  warning('"median-unweighted" is used...!!! \n something is wrong...: length of weight vector is not equal to the number of modelling methods...!')
                  .p <- apply(.p,1,mean,na.rm=TRUE)
                  .ens <- .p
                } else {
                  
                  f <- function(x) {
                    sum(weight * x)
                  }
                  
                  f2 <- function(x) {
                    if (!all(is.na(x))) {
                      .w <- which(!is.na(x))
                      .we <- weight[.w]
                      .we <- .we / sum(.we)
                      sum(.we * x[.w])
                    } else NA
                  }
                  #---------
                  .p2 <- apply(.p,1,f)
                  .wn <- which(is.na(.p2))
                  if (length(.wn) > 1) {
                    .p2[.wn] <- apply(.p[.wn,],1,f2)
                  }
                  
                  rm(.p); gc()
                  
                  .ens <- .p2
                }
              }
            } else if (method == 'median-unweighted') {
              p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              .me <- as.character(mi$method)
              
              if (length(unique(.me)) == 1) {
                warning('It seems that only one method was used to fit the models, therefore, the predictions are combined using median method...!')
                .ens <- apply(p,1,median,na.rm=TRUE)
              } else {
                .p <- data.frame(matrix(NA,nrow=nrow(p),ncol=length(unique(.me))))
                colnames(.p) <- unique(.me)
                
                for (.u in unique(.me)) {
                  w <- grep(.u,colnames(p))
                  if (length(w) > 0) {
                    if (length(w) == 1) .p[,.u] <- p[,w]
                    else .p[,.u] <- apply(p[,w],1,median,na.rm=TRUE)
                  }
                }
                #------
              
              rm(p); gc()
              
              .ens <- apply(.p,1,mean,na.rm=TRUE)
              }
            } else if (method == 'entropy') {
              .th <- .getWeights(x = x, mi = mi, wtest = wtest, id = id, stat = 'threshold', opt = opt)
              mi <- weight$mi
              id <- mi$modelID
              .th <- .th$weight
              #-----
              p <- predict(x,newdata=.nd,w=id,obj.size=obj.size,...)
              f <- function(x) {
                .x <- x >= .th
                x[] <- 0
                x[.x] <- 1
                .ent(x)
              }
              
              .p <- apply(p,1,f)
              
              rm(p); gc()
              
              .ens <- .p
            }
            
            if (filename != '') write.csv(data.frame(probability=.ens),filename,row.names = FALSE)
            
            .ens
          }
)
