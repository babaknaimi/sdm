# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  September 2024
# Version 3.8
# Licence GPL v3


#--------------------
# .raster2data.table <- function(r) {
#   if (inherits(r,'RasterBrick'))  {
#     o <- data.table(r@data@values)
#     o$cellnr <- 1:ncell(r)
#   } else if (inherits(r,'RasterLayer')) {
#     o <- data.table(r@data@values)
#     colnames(o) <- names(r)
#     o$cellnr <- 1:ncell(r)
#   } else {
#     o <- data.table(as.data.frame(r))
#     o$cellnr <- 1:ncell(r)
#   }
#   o
# }
#---------

.getlevels <- function(x) {
  o <- NULL
  if (inherits(x,'sdmdata')) {
    if (!is.null(x@factors)) {
      o <- x@factors
      names(o) <- o
      o <- lapply(o,function(i) levels(x@features[[i]]))
    }
  } else if (inherits(x,'sdmModels')) {
    
    if (!is.null(x@data@factors)) {
      o <- x@data@factors
      names(o) <- o
      o <- lapply(o,function(i) levels(x@data@features[[i]]))
    }
  } else if (inherits(x,'data.frame')) {
    f <- .where(is.factor,x)
    if (any(f)) {
      f <- names(f)[f]
      o <- vector('list',length(f))
      names(o) <- f
      for (i in seq_along(o)) o[[i]] <- levels(x[[f[i]]])
    }
  }
  o
}
#-------------

.getTotal.object.size <- function() {
  # only the size of objects in R_GlobalEnv
  paste(as.character(round(sum(unlist(lapply(ls(all.names = TRUE,envir=parent.frame()),function(x) object.size(get(x)))))/1024/1024/1024,4)),'Gb')
}

#-----------
# .raster2df <- function(x,nFact) {
#   n <- names(x)
#   d <- data.frame(getValues(x))
#   colnames(d) <- n
#   d <- data.frame(cellnr=1:ncell(x),d)
#   bb <- rep(TRUE,nrow(d))
#   for (i in 2:ncol(d)) bb <- bb & !is.na(d[,i])
#   if (length(which(bb)) == 0) stop('raster object has no data...!')
#   d <- d[bb,]
#   rm(bb)
#   
#   if (!missing(nFact) && !is.null(nFact)) {
#     for (i in seq_along(nFact)) {
#       d[,nFact[i]] <- factor(d[,nFact[i]])
#     }
#   }
#   d
# }
#----------------
.raster2df <- function(x,level) {
  d <- data.frame(cellnr=1:ncell(x),as.data.frame(x))
  bb <- rep(TRUE,nrow(d))
  for (i in 2:ncol(d)) bb <- bb & !is.na(d[,i])
  if (length(which(bb)) == 0) stop('raster object has no data...!')
  d <- d[bb,]
  rm(bb)

  if (!missing(level) && !is.null(level) && all(names(level) %in% names(d))) {
    n <- names(level)
    for (i in seq_along(n)) {
      l <- level[[i]]
      u <- sort(unique(d[,n[i]]))
      d[,n[i]] <- factor(d[,n[i]])
      levels(d[,n[i]]) <- l
      # if (any(u %in% c(1:length(l)))) {
      #   u <- u[u %in% c(1:length(l))]
      #   l <- l[u]
      #   d[,n[i]] <- factor(d[,n[i]])
      #   levels(d[,n[i]]) <- l
      # } else stop('the grid values in categorical rasters does not match with the factor levels in the model')
    }
  }
  d
}
#----------------

.generateName <- function(x) {
  paste(c(x,'_',sample(c(letters,1:9),6,replace=T)),collapse='')
}
#-----
.domClass <- function(v) {
  # decreasing sort of dominant classes in a character vector
  v <- as.character(v)
  u <- unique(v)
  names(u) <- u
  u <- unlist(lapply(u,function(x) length(which(v == x))))
  names(u)[order(u,decreasing = TRUE)]
}
#---------
if (!isGeneric("predict")) {
  setGeneric("predict", function(object, ...)
    standardGeneric("predict"))
}	
# 
# .predict=function(obj,pred,pred.par,dt=dt) {
#   pred.par[[1]] <- obj
#   pred.par[[2]] <- dt
#   m <- try(pred(pred.par),silent=TRUE)
#   options(warn=0)
#   m
# }


.generateWLP <- function(x,newdata,w=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,parallelSetting=NULL) {
  
  
  mi <- .getModel.info(x,w=w,species=species,method=method,replication=replication,run=run)
  
  s <- mi$success
  
  if (!all(s)) {
    if (!any(s)) stop('There is no model objects that were successfully fitted!')
    if (length(which(!s)) == 1) {
      warning(paste('1 model from the total of',length(s),'models was NOT successfully fitted, so it is excluded in prediction!'))
    } else {
      warning(paste(length(which(!s)),'models from the total of',length(s),'models were NOT successfully fitted, so they are excluded in prediction!'))
    }
    
    mi <- mi[s,]
    
  }
  mi <- mi[,1:5]
  #-----------
  nr <- nrow(mi)
  if (nr == 0) stop('the specified models do not exist!')
  species <- as.character(unique(mi[,2]))
  
  m <- unique(as.character(mi[,3]))
  #----
  w <- new('.workloadPredict',runTasks=mi,parallelSetting=parallelSetting)
  
  
  #w$newdata$raster <- NULL
  
  if (inherits(newdata,'data.frame')) {
    n <- colnames(newdata)
    if (!all(x@setting@featureFrame@predictors %in% n)) stop('the newdata does not contain some or all of the predictor variables required by the model...!')
    
  } else if (inherits(newdata,'Raster') || inherits(newdata,'SpatRaster')) {
    n <- names(newdata)
    if (!all(x@setting@featureFrame@predictors %in% n)) stop('the newdata does not contain some or all of the predictor variables required by the model...!')
    
  } else stop('newdata should be a Raster* (or SpatRaster) object or a data.frame...!')
  
  
  w$funs <- .sdmMethods$getPredictFunctions(m)
  w$arguments <- .sdmMethods$getPredictArguments(m)
  w$dataObject.names <- unique(unlist(lapply(x@setting@methods, .sdmMethods$getDataArgumentNames)))
  for (mo in m) {
    wc <- unlist(lapply(w$arguments[[mo]]$params,function(x) is.character(x)))
    
    if (any(!wc)) {
      if (!all(unlist(lapply(w$arguments[[mo]]$params[!wc],function(x) is.function(x))))) stop(paste('parameter definition for the model',m,'in the model container is not correctly defined!'))
      for (n in names(w$arguments[[mo]]$params[!wc])) {
        #if (!all(names(formals(w$arguments$predict[[mo]]$params[[n]])) %in% reserved.names)) stop(paste('the input argument for the function generates the parameter for model',m,'is unknown (not in the reseved objects)'))
        w$params[[n]] <- do.call(w$arguments[[mo]]$params[[n]],w$generateParams(names(formals(w$arguments[[mo]]$params[[n]])),sp))
        w$arguments$predict[[mo]]$params[[n]] <- n
      }
    }
  }
  
  #######
  w$obj <- x@models
  w$runTasks$species <- as.character(w$runTasks$species)
  w$runTasks$method <- as.character(w$runTasks$method)
  sp <- as.character(mi[,2])
  nw <- unique(sp)
  #w$runTasks$speciesID <- unlist(lapply(sp,function(x) {which(nw == x)}))
  
  m <- as.character(mi[,3])
  nw <- names(w$funs)
  #w$runTasks$methodID <- unlist(lapply(m,function(x) {which(nw == x)}))
  w$runTasks$mIDChar <- as.character(mi[,1])
  w
}

setMethod('predict', signature(object='sdmModels'), 
          function(object, newdata, filename="",id=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,mean=FALSE,overwrite=FALSE,parallelSetting=NULL,err=FALSE,...) {
            if (missing(newdata)) stop('mewdata is missing...')
            if (missing(overwrite)) overwrite <- FALSE
            if (missing(filename)) filename <- ""
            
            if (missing(parallelSetting)) parallelSetting <- NULL
            #---
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            #---
            if (missing(method) || is.null(method)) method <- object@setting@methods
            pkgs <- .sdmMethods$getPackageNames(method)
            
            
            .sdm...temp <- NULL; rm(.sdm...temp)
            pos <- 1
            
            tmp <- sapply(pkgs, function(x) '.temp' %in% x)
            
            if (any(tmp)) {
              if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
                ww <- ls(.sdmMethods$userFunctions)
                rm(list=ww,pos=1)
                rm(.sdm...temp,pos=1)
              }
              #----
              for (i in which(tmp)) pkgs[[i]] <- pkgs[[i]][which(pkgs[[i]] != '.temp')]
              #----
              tmp<- ls(.sdmMethods$userFunctions)
              
              if (length(tmp) > 0) {
                assign('.sdm...temp',c(),envir = as.environment(pos))
                for (ww in tmp) {
                  if (!exists(ww,where=1)) {
                    assign(ww,.sdmMethods$userFunctions[[ww]],envir = as.environment(pos))
                    .sdm...temp <<- c(.sdm...temp,ww)
                  }
                }
              }
            }
            
            ww <- .loadLib(pkgs)
            if (!all(ww)) {
              if (!any(ww)) {
                stop(paste('There is no installed packages rquired by the selected methods. Package names:',paste(unlist(pkgs),collapse=', ')))
              } else {
                warning(paste('There is no installed packages rquired by the methods:',paste(object@setting@methods[!ww],collapse=', '),'; These methods are excluded! The packages need to be installed for these methods:',paste(unlist(pkgs[!ww]),collapse=', ')))
                method <- method[ww] 
              }
            }
            
            if (filename != "") {
              
              if (file.exists(filename) && !overwrite) stop('filename exists... (use overwrite=TRUE or use a different filename!)')
            }
            #----
            if (!missing(parallelSetting) && is.list(parallelSetting)) {
              nparallel <- names(parallelSetting)
              a <- c('ncore','doParallel','method','cluster','hosts','fork','type','strategy')
              nparallel <- .pmatch(nparallel,a)
              w <- which(!is.na(nparallel))
              if (length(w) > 0) {
                parallelSetting <- parallelSetting[w]
                nparallel <- nparallel[w]
                names(parallelSetting) <- nparallel
              }
              #--
              if ('cluster' %in% nparallel && inherits(parallelSetting$cluster,'cluster')) .parSetting@cl <- parallelSetting$cluster
              #--
              .parSetting <- new('.parallelSetting')
              if ('ncore' %in% nparallel) .parSetting@ncore <- min(c(parallelSetting$ncore,parallel::detectCores()))
              else {
                if (!is.null(.parSetting@cl)) .parSetting@ncore <- length(.parSetting@cl)
                else .parSetting@ncore <- max(c(floor(parallel::detectCores() * 0.5),1))
              }
              #--
              if ('method' %in% nparallel) {
                if (!parallelSetting$method %in% c('parallel','foreach','future')) {
                  warning('The parallel method is not recognised; the default value ("parallel") is used!')
                  .parSetting@method <- 'parallel'
                } else .parSetting@method <- 'parallel'
              } else .parSetting@method <- 'parallel'
              #--
              if ('type' %in% nparallel) .parSetting@type <- parallelSetting$type
              #--
              if ('strategy' %in% nparallel) {
                parallelSetting$strategy <- tolower(parallelSetting$strategy)[1]
                if (!parallelSetting$strategy %in% c('data','model','auto')) {
                  warning('The parallel strategy is not recognised (should be one of c("auto","data","model")); the default, "auto", is used!')
                  .parSetting@strategy <- 'auto'
                } else .parSetting@strategy <- parallelSetting$strategy
              } else .parSetting@strategy <- 'auto'
              #--
              if ('doParallel' %in% nparallel && is.expression(parallelSetting$doParallel)) .parSetting@doParallel <- parallelSetting$doParallel
              #--
              
              if ('hosts' %in% nparallel && is.character(parallelSetting$hosts)) .parSetting@hosts <- parallelSetting$hosts
              
              #----
              if ('fork' %in% nparallel) {
                if (is.logical(parallelSetting$fork)) {
                  .parSetting@fork <- parallelSetting$fork
                  if (parallelSetting$fork && .is.windows()) {
                    warning('"fork" in the parallel setting cannot be TRUE on Windows Operating Systems; It is changed to FALSE!')
                    .parSetting@fork <- FALSE
                  }
                } else {
                  if (is.null(parallelSetting$fork)) .parSetting@fork <- !.is.windows()
                  else {
                    warning('"fork" in parallel setting should be logical; the default value is used!')
                    .parSetting@fork <- !.is.windows()
                  }
                }
              } else .parSetting@fork <- !.is.windows()
              #--
              parallelSetting <- .parSetting
            } else {
              if (!is.null(parallelSetting)) {
                if (!inherits(parallelSetting,'.parallelSetting')) {
                  warning('parallelSetting should be provided as a list (it is ignored)!')
                  parallelSetting <- NULL
                }
              }
            }
            #----
            if (inherits(newdata,'Raster')) newdata <- rast(newdata)
            #-----
            w <- .generateWLP(x = object,newdata=newdata,w=id,species=species,method=method,replication=replication,run=run,parallelSetting=parallelSetting)
            
            mid <- w$runTasks$modelID
            #----
            nr <- nrow(w$runTasks)
            ###############################
            
            if (!is.null(parallelSetting) && parallelSetting@ncore > 1) {
              .require('parallel')
              
              if (parallelSetting@method == 'future') {
                if (!.require('future.apply')) {
                  warning('To use "future" as the parallel method, packages of "future" and "future.apply" are needed but they are not available, so the method is changed to "parallel"')
                  parallelSetting@method <- "parallel"
                }
              } else if (parallelSetting@method == 'foreach') {
                if (!.require('foreach')) {
                  warning('To use "foreach" as the parallel method, the "foreach" package is needed but it is not available, so the method is changed to "parallel"')
                  parallelSetting@method <- "parallel"
                }
              }
              
              
              if (parallelSetting@method == 'future') {
                if (parallelSetting@fork && !.is.windows()) {
                  if (!is.null(parallelSetting@cl)) {
                    .eval('plan(cluster,workers=parallelSetting@cl,gc=TRUE)',environment())
                  } else {
                    .eval('plan(multicore,workers=parallelSetting@ncore,gc=TRUE)',environment())
                  }
                  
                } else {
                  if (!is.null(parallelSetting@cl)) {
                    .eval('plan(cluster,workers=parallelSetting@cl,gc=TRUE)',environment())
                  } else {
                    .eval('plan(multisession,workers=parallelSetting@ncore,gc=TRUE)',environment())
                  }
                  
                }
                .lapply <- .fuLapply
                .cl <- NULL
                #on.exit({.eval('plan(sequential,gc=TRUE)',environment());gc()}, add=TRUE)
              } else {
                if (!is.null(parallelSetting@cl)) {
                  .cl <- parallelSetting@cl
                  clusterExport(parallelSetting@cl,c('w','mid'),envir = environment())
                  
                  .cle <- try(clusterEvalQ(parallelSetting@cl,{
                    library(sdm)
                    sdm:::.addMethods()
                    sdm:::.pkgLoad(names(w$funs))
                    
                  }),silent = TRUE)
                  
                  if (inherits(.cle,'try-error')) {
                    parallelSetting@ncore
                    warning('parallelising is not working for the predict function (it is ignored)!')
                    cat('\n ncore is changed to 1!')
                    .lapply <- .Lapply
                  } else {
                    .lapply <- .parLapply
                  }
                  
                  
                } else {
                  if (parallelSetting@fork && !.is.windows()) {
                    .cl <- makeForkCluster(parallelSetting@ncore)
                  } else {
                    if (!is.null(parallelSetting@hosts)) {
                      .cl <- try(makePSOCKcluster(parallelSetting@hosts),silent = TRUE)
                      # it should be updated given settings for remote connections
                      # https://parallelly.futureverse.org/reference/makeClusterPSOCK.html
                      if (inherits(.cl,'try-error')) {
                        cat('\n Error in connecting to remote servers:' ,print(.cl))
                        cat('\n cores on the local machine is considered!')
                        .cl <- makePSOCKcluster(parallelSetting@ncore)
                      }
                    } else .cl <- makePSOCKcluster(parallelSetting@ncore)
                    #-----------
                    clusterExport(.cl,c('w','mid'),envir = environment())
                    
                    .cle <- try(clusterEvalQ(.cl,{
                      library(sdm)
                      sdm:::.addMethods()
                      sdm:::.pkgLoad(names(w$funs))
                      
                    }),silent = TRUE)
                    
                    if (inherits(.cle,'try-error')) {
                      parallelSetting@ncore
                      stopCluster(.cl)
                      warning('parallelising is not working for the predict function (it is ignored)!')
                      cat('\n ncore is changed to 1!')
                      .lapply <- .Lapply
                    } else .lapply <- .parLapply
                    
                  }
                  parallelSetting@cl <- .cl
                }
                
                #on.exit({stopCluster(.cl);gc()}, add=TRUE)
              }
              #----------
              if (parallelSetting@method == 'foreach' && .require('foreach') && (.require('doParallel') | .eval('getDoParRegistered()',environment()))) {
                if (!.eval('getDoParRegistered()',environment())) {
                  if (!is.null(parallelSetting@doParallel) && is.expression(parallelSetting@doParallel)) {
                    .tmp <- try(eval(parallelSetting@doParallel),silent = TRUE)
                    if (inherits(.tmp,'try-error') || !.eval('getDoParRegistered()',environment())) {
                      eval(expression({registerDoParallel(parallelSetting@cl)}))
                    }
                  } else eval(expression({registerDoParallel(parallelSetting@cl)}))
                }
                .lapply <- .feLapply
              }
            }
            ###############################
            success <- rep(FALSE,nr)
            rnames <- fullnames <- c()
            errLog <- list()
            
            if (nr > 1 && mean && !is.na(w$runTasks$replication[1]) && length(unique(w$runTasks$replicationID)) > 1) {
              .mi <- w$runTasks[0,c(2:4)]
              for (.sp in unique(w$runTasks$species)) {
                .tmp <- w$runTasks[w$runTasks$species == .sp,]
                for (.me in unique(.tmp$method)) {
                  .tmp2 <- .tmp[.tmp$method == .me,]
                  .mi <- rbind(.mi,data.frame(species=.sp,method=.me,replication=unique(.tmp2$replication)))
                }
              }
              nr <- nrow(.mi)
              
              for (j in 1:nr) {
                rnames <- c(rnames,paste0('sp_',.mi$species[j],'__m_',.mi$method[j],paste0('__re_',paste(strsplit(.mi$replication[j],'')[[1]][1:4],collapse=''))))
                fullnames <- c(fullnames,paste0('species_',.mi$species[j],'-method_',.mi$method[j],paste0('-replication (Mean)_',.mi$replication[j])))
              }
              success <- rep(FALSE,nr)
            } else {
              mean <- FALSE
              for (j in 1:nr) {
                rnames <- c(rnames,paste('id_',w$runTasks$modelID[j],'__sp_',w$runTasks$species[j],'__m_',w$runTasks$method[j],if (!is.na(w$runTasks$replication[j])) paste('__re_',paste(strsplit(w$runTasks$replication[j],'')[[1]][1:4],collapse=''),sep=''),sep=''))
                fullnames <- c(fullnames,paste('id_',w$runTasks$modelID[j],'-species_',w$runTasks$species[j],'-method_',w$runTasks$method[j],if (!is.na(w$runTasks$replication[j])) paste('-replication_',w$runTasks$replication[j],sep=''),sep=''))
              }
            }
            #---------
            
            #options(warn=-1)
            
            if (inherits(newdata,'SpatRaster')) {
              .out <- rast(newdata,nlyrs=nr)
              .nc <- ncol(newdata)
              names(.out) <- rnames
              metags(.out) <- cbind(rnames,fullnames)
              
              if (filename == "" && .canProcessInMemory(.out[[1]],max(nlyr(.out), nlyr(newdata))*2)) {
                .d <- as.data.frame(newdata,cells=TRUE,na.rm=TRUE)
                .c <- .d$cell
                .d <- object@setting@featureFrame@featureGenerator(.d) # features
                #---
                if (!is.null(parallelSetting) && parallelSetting@ncore > 1) {
                  if (parallelSetting@strategy == 'auto') {
                    if (length(unique(w$runTasks$method)) >= parallelSetting@ncore &&  (nrow(.d) / parallelSetting@ncore) < 20000) parallelSetting@strategy <- 'model'
                    else parallelSetting@strategy <- 'data'
                  }
                  #---------
                  if (parallelSetting@strategy == 'data') {
                    .ds <- split(.d, rep(1:parallelSetting@ncore, each=ceiling(nrow(.d)/parallelSetting@ncore), length.out=nrow(.d)))
                    .ds <- try(.lapply(.ds,function(x) w$predictMID(mid,.frame = x),cl = .cl),silent = TRUE)
                    
                    if (!inherits(.ds,'try-error')) {
                      .d <- vector('list',length(mid))
                      for (j in 1:length(.ds[[1]])) {
                        for (k in 1:length(.ds)) {
                          .d[[j]] <- c(.d[[j]],.ds[[k]][[j]])
                        }
                      }
                    } else {
                      warning('error in parallelized procedure (non-parallel predict is re-executed...)!')
                      .d <- w$predictMID(IDs = mid,.frame = .d)
                    }
                    rm(.ds); gc()
                  } else {
                    .ds <- try(.lapply(mid,function(i,.d) w$predictID(i,.frame = .d),cl=.cl,.d=.d),silent = TRUE)
                    if (!inherits(.ds,'try-error')) .d <- .ds
                    else {
                      warning('error in parallelized procedure (non-parallel predict is re-executed...)!')
                      .d <- w$predictMID(IDs = mid,.frame = .d)
                    }
                  }
                } else {
                  .d <- w$predictMID(IDs = mid,.frame = .d)
                }
                #----
                .success <- sapply(.d, function(x) !inherits(x,'try-error'))
                names(.d) <- names(.success) <- mid
                #---
                if (mean) {
                  for (j in 1:nrow(.mi)) {
                    .w <- as.character(mid[w$runTasks$species == .mi$species[j] & w$runTasks$method == .mi$method[j] & w$runTasks$replication == .mi$replication[j]])
                    
                    if (length(which(.success[.w])) > 1) {
                      set.values(.out,.c,rowMeans(data.frame(.d[.w[.success[.w]]]),na.rm=TRUE),layer=j)
                      success[j] <- TRUE
                    } else if (length(which(.success[.w])) == 1) {
                      set.values(.out,.c,.d[[.w[.success[.w]]]],layer=j)
                      success[j] <- TRUE
                    } else {
                      errLog <- c(errLog,paste0('Error in predictions for method ',.mi$method[j],' and replication ',.mi$replication[j]))
                    }
                  }
                } else {
                  for (j in 1:length(.d)) {
                    if (.success[j]) {
                      set.values(.out,.c,.d[[j]],layer=j)
                      success[j] <- TRUE
                    } else {
                      errLog <- c(errLog,.d[[j]])
                    }
                  }
                }
                
                
              } else {
                readStart(newdata)
                on.exit(readStop(newdata))
                
                b <- writeStart(.out, filename=filename, overwrite=overwrite, n=max(nlyr(.out), nlyr(newdata))*2, sources=sources(newdata),...)
                for (i in 1:b$n) {
                  .d <- readValues(newdata, b$row[i], b$nrows[i], 1, .nc, TRUE, TRUE)
                  .c <- .getCells(.nc,b$row[i], b$nrows[i])
                  .wna <- which(!apply(.d,1,function(x) any(is.na(x))))
                  if (length(.wna) > 0) {
                    .d <- object@setting@featureFrame@featureGenerator(.d[.wna,]) # features
                    if (!is.null(parallelSetting) && parallelSetting@ncore > 1) {
                      if (parallelSetting@strategy == 'auto') {
                        if (length(unique(w$runTasks$method)) >= parallelSetting@ncore &&  (nrow(.d) / parallelSetting@ncore) < 20000) parallelSetting@strategy <- 'model'
                        else parallelSetting@strategy <- 'data'
                      }
                      #---------
                      if (parallelSetting@strategy == 'data') {
                        .ds <- split(.d, rep(1:parallelSetting@ncore, each=ceiling(nrow(.d)/parallelSetting@ncore), length.out=nrow(.d)))
                        .ds <- try(.lapply(.ds,function(x) w$predictMID(mid,.frame = x),cl = .cl),silent = TRUE)
                        
                        if (!inherits(.ds,'try-error')) {
                          .d <- vector('list',length(mid))
                          for (j in 1:length(.ds[[1]])) {
                            for (k in 1:length(.ds)) {
                              .d[[j]] <- c(.d[[j]],.ds[[k]][[j]])
                            }
                          }
                        } else {
                          warning('error in parallelized procedure (non-parallel predict is re-executed...)!')
                          .d <- w$predictMID(IDs = mid,.frame = .d)
                        }
                        rm(.ds); gc()
                      } else {
                        .ds <- try(.lapply(mid,function(i,.d) w$predictID(i,.frame = .d),cl=.cl,.d=.d),silent = TRUE)
                        if (!inherits(.ds,'try-error')) .d <- .ds
                        else {
                          warning('error in parallelized procedure (non-parallel predict is re-executed...)!')
                          .d <- w$predictMID(IDs = mid,.frame = .d)
                        }
                      }
                    } else {
                      .d <- w$predictMID(IDs = mid,.frame = .d)
                    }
                    
                    .tmp <- rep(NA,length(.c))
                    .v <- c()
                    
                    .success <- sapply(.d, function(x) !inherits(x,'try-error'))
                    names(.d) <- names(.success) <- mid
                    
                    if (mean) {
                      for (j in 1:nrow(.mi)) {
                        .vtmp <- NULL
                        .w <- as.character(mid[w$runTasks$species == .mi$species[j] & w$runTasks$method == .mi$method[j] & w$runTasks$replication == .mi$replication[j]])
                        
                        if (length(which(.success[.w])) > 1) {
                          .vtmp <- .tmp
                          .vtmp[.wna] <- rowMeans(data.frame(.d[.w[.success[.w]]]),na.rm=TRUE)
                          .v <- c(.v,.vtmp)
                          success[j] <- TRUE
                        } else if (length(which(.success[.w])) == 1) {
                          .vtmp <- .tmp
                          .vtmp[.wna] <- .d[[.w[.success[.w]]]]
                          .v <- c(.v,.vtmp)
                          success[j] <- TRUE
                        } else {
                          .v <- c(.v,.tmp)
                          errLog <- c(errLog,paste0('Error in predictions for method ',.mi$method[j],' and replication ',.mi$replication[j]))
                        }
                      }
                      
                    } else {
                      .vtmp <- NULL
                      for (j in 1:length(.d)) {
                        if (.success[j]) {
                          .vtmp <- .tmp
                          .vtmp[.wna] <- .d[[j]]
                          .v <- c(.v,.vtmp)
                          success[j] <- TRUE
                        } else {
                          .v <- c(.v,.tmp)
                          errLog <- c(errLog,.d[[j]])
                        }
                      }
                    }
                    #---
                    if (length(.v) != prod(b$nrows[i], .nc, nr)) {
                      msg <- "the number of values returned by the predict does not match the number of cells in newdata!"
                      writeStop(.out)
                      rm(.out)
                      if (filename != "") unlink(filename)
                      gc()
                      stop(paste("predict", msg))
                    }
                    writeValues(.out, .v, b$row[i], b$nrows[i])
                    rm(.tmp,.vtmp,.d,.wna,.c,.v); gc()
                  }
                  
                }
                
                writeStop(.out)
                
              }
              #----
              if (!all(success)) {
                if (!any(success)) {
                  rm(.out)
                  if (filename != "") unlink(filename)
                  gc()
                  stop('Error in prediction....!')
                }
                warning(paste0('Failed prediction for ',length(which(!success)),' models (out of ',length(success),')...!'))
                .out <- subset(.out,which(success),filename=filename,overwrite=TRUE)
                mid <- mid[success]
                w$runTasks <-  w$runTasks[success,]
              }
              #----------
            } else if (inherits(newdata,'data.frame')) {
              .out <- matrix(NA,nrow=nrow(newdata),ncol=0)
              #----
              .d <- object@setting@featureFrame@featureGenerator(newdata) # features
              #############
              # parallel:
              if (!is.null(parallelSetting) && parallelSetting@ncore > 1) {
                if (parallelSetting@strategy == 'auto') {
                  if (length(unique(w$runTasks$method)) >= parallelSetting@ncore &&  (nrow(.d) / parallelSetting@ncore) < 20000) parallelSetting@strategy <- 'model'
                  else parallelSetting@strategy <- 'data'
                }
                #---------
                if (parallelSetting@strategy == 'data') {
                  .ds <- split(.d, rep(1:parallelSetting@ncore, each=ceiling(nrow(.d)/parallelSetting@ncore), length.out=nrow(.d)))
                  .ds <- try(.lapply(.ds,function(x) w$predictMID(mid,.frame = x),cl = .cl),silent = TRUE)
                  
                  if (!inherits(.ds,'try-error')) {
                    .d <- vector('list',length(mid))
                    for (j in 1:length(.ds[[1]])) {
                      for (k in 1:length(.ds)) {
                        .d[[j]] <- c(.d[[j]],.ds[[k]][[j]])
                      }
                    }
                  } else {
                    warning('error in parallelized procedure (non-parallel predict is re-executed...)!')
                    .d <- w$predictMID(IDs = mid,.frame = .d)
                  }
                  rm(.ds); gc()
                } else {
                  .ds <- try(.lapply(mid,function(i,.d) w$predictID(i,.frame = .d),cl=.cl,.d=.d),silent = TRUE)
                  if (!inherits(.ds,'try-error')) .d <- .ds
                  else {
                    warning('error in parallelized procedure (non-parallel predict is re-executed...)!')
                    .d <- w$predictMID(IDs = mid,.frame = .d)
                  }
                }
              } else {
                .d <- w$predictMID(IDs = mid,.frame = .d)
              }
              ##############
              
              .success <- sapply(.d, function(x) !inherits(x,'try-error'))
              
              if (mean) {
                for (j in 1:nrow(.mi)) {
                  .vtmp <- NULL
                  .w <- mid[w$runTasks$species == .mi$species[j] & w$runTasks$method == .mi$method[j] & w$runTasks$replication == .mi$replication[j]]
                  
                  if (length(which(.success[.w])) > 1) {
                    .vtmp <- rowMeans(data.frame(.d[.w[.success[.w]]]),na.rm=TRUE)
                    .out <- cbind(.out,.vtmp)
                    success[j] <- TRUE
                  } else if (length(which(.success[.w])) == 1) {
                    .vtmp <- .d[[.w[.success[.w]]]]
                    .out <- cbind(.out,.vtmp)
                    success[j] <- TRUE
                  } else {
                    errLog <- c(errLog,paste0('Error in predictions for method ',.mi$method[j],' and replication ',.mi$replication[j]))
                  }
                  
                  
                  
                  for (k in 1:length(.w)) {
                    if (.success[.w[k]]) {
                      if (is.null(.vtmp)) .vtmp <- .d[[.w[k]]]
                      else .vtmp <- .vtmp + .d[[.w[k]]]
                    }
                  }
                  #---
                  if (!is.null(.vtmp)) {
                    .vtmp <- .vtmp / length(.w[.success[.w]])
                    .out <- cbind(.out,.vtmp)
                    success[j] <- TRUE
                  } else {
                    errLog <- c(errLog,paste0('Error in predictions for method ',.mi$method[j],' and replication ',.mi$replication[j]))
                  }
                }
                #---
                
              } else {
                for (j in 1:length(.d)) {
                  if (.success[j]) {
                    .out <- cbind(.out,.d[[j]])
                    success[j] <- TRUE
                  } else {
                    errLog <- c(errLog,.d[[j]])
                  }
                }
              }
              #---
              if (!all(success)) {
                if (!any(success)) {
                  rm(.out)
                  stop('Error in prediction....!')
                }
                warning(paste0('Failed prediction for ',length(which(!success)),' models (out of ',length(success),')...!'))
                rnames <- rnames[success]
                
                mid <- mid[success]
                w$runTasks <-  w$runTasks[success,]
              }
              #---
              colnames(.out) <- rnames
              .out <- data.frame(.out)
              if (filename != "") write.csv(.out,filename,row.names = FALSE)
            }
            
            
            #------------
            
            if (err && length(errLog) > 0) {
              for (i in seq_along(errLog)) cat(errLog[[i]],'\n')
            }
            
            if (".sdm...temp" %in% ls(pattern='^.sdm..',pos=1,all.names = TRUE)) {
              .w <- ls(.sdmMethods$userFunctions)
              rm(list=.w,pos=1)
              rm(.sdm...temp,pos=1)
            }
            #---
            
            
            if (!is.null(parallelSetting)) {
              if (parallelSetting@method == 'future')  {
                .w <- try(.eval('plan(sequential,gc=TRUE)',environment()),silent = TRUE)
              } else  {
                .w <- try(stopCluster(.cl),silent = TRUE)
              }
            }
            gc()
            
            return(.out)
          }
)