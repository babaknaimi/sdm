# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Jan 2024
# Version 1.2
# Licence GPL v3
#--------


if (!isGeneric("sdmSetting")) {
  setGeneric("sdmSetting", function(formula,data,methods,interaction.depth=1,n=1,replication=NULL,
                                    cv.folds=NULL,test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
                                    var.selection=FALSE,modelSettings=NULL,seed=NULL,parallelSetting=NULL,...)
    standardGeneric("sdmSetting"))
}

setMethod('sdmSetting', signature(formula='ANY','sdmdata','character'), 
          function(formula,data,methods,interaction.depth=1,n=1,replication=NULL,
                   cv.folds=NULL,test.percent=NULL,bg=NULL,bg.n=NULL,var.importance=NULL,response.curve=TRUE,
                   var.selection=FALSE,modelSettings=NULL,seed=NULL,parallelSetting=NULL,...) {
            
            if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
            
            dot <- list(...)
            sobj <- NULL
            if (length(dot) > 0) {
              ndot <- names(dot)
              if ('' %in% ndot) {
                for (i in seq_along(which(ndot == ''))) {
                  if (inherits(dot[[i]],'.sdmCorSetting')) {
                    sobj <- dot[[i]]
                    break
                  }
                }
                dot <- dot[-which(ndot == '')]
                ndot <- names(dot)
              }
              
              #a <- c('interaction.depth','replication','cv.folds','test.percent','bg','bg.n','var.importance','response.curve','var.selection','ncore','modelSettings','seed','setting','parallelSetting')
              #ndot <- .pmatch(ndot,a)
              #w <- !is.na(ndot)
              # if (length(w) > 0) {
              #   dot <- dot[w]
              #   ndot <- ndot[w]
              #   names(dot) <- ndot
              # }
              
              if ('setting' %in% names(dot) && inherits(dot[['setting']],'.sdmCorSetting')) {
                sobj <- dot[['setting']]
                dot <- dot[-which(ndot == 'setting')]
                ndot <- names(dot)
              }
            }
            #--------
            
            m <- .methodFix(methods)
            if (any(is.na(m))) stop(paste('methods',paste(methods[is.na(m)],collapse=', '),'do not exist!'))
            m <- unique(m)
            #---------
            s <- new('.sdmCorSetting',methods=m)
            #---------
            if (missing(formula)) {
              if (!is.null(sobj)) {
                if (all(sobj@sdmFormula@vars@names %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
                else s@sdmFormula <- data@sdmFormula
              } else s@sdmFormula <- data@sdmFormula
              
            } else if (inherits(formula,'sdmFormula')) s@sdmFormula <- formula
            else if (inherits(formula,'formula')) {
              s@sdmFormula <- .exFormula(formula,as.data.frame(data)[,-1])
              if (is.null(s@sdmFormula@data.terms) ) {
                if (!is.null(data@sdmFormula@data.terms)) {
                  .tmp <- sapply(data@sdmFormula@data.terms,class)
                  if (any(c(".scaleSetting",".pcaSetting") %in% .tmp)) {
                    s@sdmFormula@data.terms <- data@sdmFormula@data.terms[.tmp %in% c(".scaleSetting",".pcaSetting")]
                  }
                }
              } else {
                if (!is.null(data@sdmFormula@data.terms)) {
                  .tmp <- sapply(data@sdmFormula@data.terms,class)
                  if (any(c(".scaleSetting",".pcaSetting") %in% .tmp)) {
                    .tmp2 <- sapply(s@sdmFormula@data.terms,class)
                    w <- c(".scaleSetting",".pcaSetting")[c(".scaleSetting",".pcaSetting") %in% .tmp]
                    if (any(!w %in% .tmp2)) {
                      w <- w[!w %in% .tmp2]
                      s@sdmFormula@data.terms <- c(s@sdmFormula@data.terms,data@sdmFormula@data.terms[.tmp %in% w])
                    }
                  }
                }
              }
            } else if (inherits(formula,'.sdmCorSetting')) {
              sobj <- formula
              if (all(sobj@sdmFormula@vars@names %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
              else s@sdmFormula <- data@sdmFormula
            } else {
              if (!is.null(sobj)) {
                if (all(sobj@sdmFormula@vars@names %in% data@features.name)) s@sdmFormula <- sobj@sdmFormula
                else s@sdmFormula <- data@sdmFormula
              } else s@sdmFormula <- data@sdmFormula
            }
            
            s@featureFrame <- .getFeatureFrame(s@sdmFormula,data = as.data.frame(data)[,-1]) #.getFeaturetype(data,s@sdmFormula)  
            #---------
            s@distribution <- .getSpeciesDistribution(data,sp=s@sdmFormula@vars@species)
            #---------
            if (!is.null(test.percent)) s@test.percentage <- test.percent
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@test.percentage)) s@test.percentage <- sobj@test.percentage
              }
            }
            #---------
            if (!missing(parallelSetting) && !is.null(parallelSetting) && is.list(parallelSetting)) {
              nparallel <- names(parallelSetting)
              a <- c('ncore','doParallel','method','cluster','hosts','fork','type')
              nparallel <- .pmatch(nparallel,a)
              w <- which(!is.na(nparallel))
              if (length(w) > 0) {
                parallelSetting <- parallelSetting[w]
                nparallel <- nparallel[w]
                names(parallelSetting) <- nparallel
              }
              #--
              if ('ncore' %in% nparallel) s@parallelSetting@ncore <- parallelSetting$ncore
              else {
                if (!is.null(sobj) && length(sobj@parallelSetting@ncore) == 1) s@parallelSetting@ncore <- sobj@parallelSetting@ncore
                else s@parallelSetting@ncore <- max(c(floor(parallel::detectCores() * 0.5),1))
              }
              #--
              if ('method' %in% nparallel) {
                if (parallelSetting$method %in% c('parallel','foreach','future')) s@parallelSetting@method <- parallelSetting$method
                else {
                  warning('The parallel method is not recognised; the default value ("parallel") is used!')
                  s@parallelSetting@method <- 'parallel'
                }
              } else s@parallelSetting@method <- 'parallel'
              #--
              if ('fork' %in% nparallel) {
                if (is.logical(parallelSetting$fork)) {
                  if (parallelSetting$fork && .is.windows()) {
                    warning('"fork" in the parallel setting cannot be TRUE on Windows Operating Systems; It is changed to FALSE!')
                    s@parallelSetting@fork <- FALSE
                  } else s@parallelSetting@fork <- parallelSetting$fork
                } else {
                  warning('"fork" in parallel setting should be logical; the default value is used!')
                  s@parallelSetting@fork <- !.is.windows()
                }
              } else s@parallelSetting@fork <- !.is.windows()
              #--
              if ('strategy' %in% nparallel) {
                parallelSetting$strategy <- tolower(parallelSetting$strategy)[1]
                if (!parallelSetting$strategy %in% c('species','method','replicate','simple','auto')) {
                  warning('The parallel strategy is not recognised (should be one of c("auto","species","method","replicate","simple")); the default, "auto", is used!')
                  s@parallelSetting@strategy <- 'auto'
                } else s@parallelSetting@strategy <- parallelSetting$strategy
              } else s@parallelSetting@strategy <- 'auto'
              #---
              if ('type' %in% nparallel) s@parallelSetting@type <- parallelSetting$type
              #--
              if ('doParallel' %in% nparallel && is.expression(parallelSetting$doParallel)) s@parallelSetting@doParallel <- parallelSetting$doParallel
              #--
              if ('cluster' %in% nparallel && inherits(parallelSetting$cluster,'cluster')) s@parallelSetting@cl <- parallelSetting$cluster
              #--
              if ('hosts' %in% nparallel && is.character(parallelSetting$hosts)) s@parallelSetting@hosts <- parallelSetting$hosts
              
            } else {
              if (!is.null(sobj)) s@parallelSetting <- sobj@parallelSetting
              else {
                if (length(dot) > 0 && 'ncore' %in% ndot) {
                  if (is.numeric(dot[['ncore']])) {
                    s@parallelSetting@ncore <- dot[['ncore']][1]
                    s@parallelSetting@method <- 'parallel'
                    s@parallelSetting@fork <- !.is.windows()
                  } else s@parallelSetting@ncore <- 1
                } else s@parallelSetting@ncore <- 1
              } 
            }
            
            #---------
            if (!is.null(replication)) {
              nx <- .replicate.methodFix(replication)
              if (any(is.na(nx))) warning(paste(paste(replication[is.na(nx)],collapse=', '),'methods in replication are not found [They are ignored!]'))
              replication <- nx[!is.na(nx)]
              s@replicate <- replication
            } else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@replicate)) s@replicate <- sobj@replicate
              }
              if (is.null(s@replicate) && !is.null(s@test.percentage)) {
                s@replicate <- "subsampling"
              }
            }
            
            s@n.replicates <- n
            if (!is.null(sobj) && !is.null(sobj@n.replicates)) s@n.replicates <- sobj@n.replicates
            
            if ("subsampling" %in% s@replicate) {
              if (is.null(s@test.percentage)) s@test.percentage <- 30
            }
            
            if (!is.null(cv.folds)) s@cv.folds <- cv.folds
            else {
              if (!is.null(sobj) && !is.null(sobj@cv.folds)) s@cv.folds <- sobj@cv.folds
              if (is.null(s@cv.folds) && "cross_validation" %in% s@replicate) s@cv.folds <- 5
            }
            
            if (!is.null(s@cv.folds) && !"cross_validation" %in% s@replicate) {
              s@replicate <- c("cross_validation",s@replicate)
            }
            #---------
            if (!is.null(bg)) s@pseudo.absence.methods <- bg
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@pseudo.absence.methods)) s@pseudo.absence.methods <- sobj@pseudo.absence.methods
              }
            }
            if (!is.null(bg.n)) s@n.pseudo.absence <- bg.n
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@n.pseudo.absence)) s@n.pseudo.absence <- sobj@n.pseudo.absence
              }
              if (is.null(s@n.pseudo.absence) && !is.null(s@pseudo.absence.methods)) {
                s@n.pseudo.absence <- 1000
              }
            }
            #---------
            if (!is.null(var.importance)) s@varImportance.methods <- var.importance
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@varImportance.methods)) s@varImportance.methods <- sobj@varImportance.methods
              }
            }
            #---------
            if (response.curve) s@response.curve <- TRUE
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@response.curve) && sobj@response.curve) s@response.curve <- sobj@response.curve
              } else s@response.curve <- FALSE
            }
            #---------
            if (var.selection) s@var.selection <- TRUE
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@var.selection) && sobj@var.selection) s@var.selection <- sobj@var.selection
              } else s@var.selection <- FALSE
            }
            #---------
            #s@interaction.depth <- interaction.depth
            #if (interaction.depth ==1 && !is.null(sobj) && !is.null(sobj@interaction.depth)) s@interaction.depth <- sobj@interaction.depth
            #---------
            
            if (!is.null(interaction.depth)) s@interaction.depth <- interaction.depth
            else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@interaction.depth)) s@interaction.depth <- sobj@interaction.depth
              }
            }
            #---------
            if (!is.null(modelSettings) && inherits(modelSettings,'list')) {
              .ms <- names(modelSettings)
              if (!is.null(.ms)) {
                .ms <- .methodFix(.ms)
                if (!all(.ms %in% s@methods)) warning(paste('the models in the modelSettings:',paste0(names(modelSettings)[!.ms %in% s@methods],collapse = ', '),'are not selected in the methods, or do not exitst!'))
                w <- which(.ms %in% s@methods)
                if (length(w) > 0) {
                  .ms <- .ms[w]
                  modelSettings <- modelSettings[w]
                  names(modelSettings) <- .ms
                  ww <- c()
                  for (i in seq_along(.ms)) {
                    if(!inherits(modelSettings[[.ms[i]]],'list')) ww <- c(ww,i)
                  }
                  
                  if (length(ww) > 0) {
                    if (length(ww) < length(modelSettings)) {
                      warning(paste('the modelSettings for the items:',paste(.ms[ww],collapse = ','),'are not a list, and so they are ignored!'))
                      modelSettings <- modelSettings[-ww]
                    } else {
                      warning('the arguments for each method in the modelSettings should be introduced using a list; modelSettings is ignored!')
                      modelSetting <- NULL
                    }
                  }
                } else modelSetting <- NULL
                
              } else warning('modelSettings is not in the right structure, so it is ignored!')
              
              if (!is.null(modelSettings)) {
                s@modelSettings <- modelSettings
              }
            } else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@modelSettings)) s@modelSettings <- sobj@modelSettings
              }
            }
            #---------
            if (!is.null(seed)) {
              if (is.logical(seed)) seed <- sample(100000,1)
              else if (!is.numeric(seed)) seed <- NULL
              s@seed <- seed
            } else {
              if (!is.null(sobj)) {
                if (!is.null(sobj@seed)) s@seed <- sobj@seed
              }
            }
            #-------------
            s
          }
)
