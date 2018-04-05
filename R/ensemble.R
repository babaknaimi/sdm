# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Oct. 2016
# Version 2.0
# Licence GPL v3

if (!isGeneric("ensemble")) {
  setGeneric("ensemble", function(x,newdata,filename,method,...)
    standardGeneric("ensemble"))
}	


setMethod('ensemble', signature(x='sdmModels',newdata='Raster'), 
          function(x, newdata, filename,setting,...) {
            if (missing(setting)) setting <- list(method='weighted',stat='AUC',wtest=NULL)
            if (missing(filename)) filename=''
            
            if (filename == '') filename <- .generateName('sdm_ensemble')
            if (extension(filename) == '') filename <- paste(filename,'.grd',sep='')
            
            if (!is.list(setting)) stop('setting should be defined as a list object!')
            
            n <- names(setting)
            if (is.null(n) || '' %in% n) stop('one or more arguments in setting do not have name name!')
            
            n <- .pmatch(n,c('method','stat','opt','id','wtest','weight'))
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
            
            if (method == 'weighted') {
              if (is.null(weight)) {
                if (all(mi[,wtest[1]])) {
                  weight <- getEvaluation(x,w = id,wtest = wtest[1],stat = stat,opt = opt)[,2]
                } else {
                  id1 <- mi$modelID[mi[,wtest[1]]]
                  if (length(wtest) > 1) {
                    w <- which(!mi[,wtest[1]])
                    if (any(mi[w,wtest[2]])) {
                      id2 <- mi$modelID[w][mi[w,wtest[2]]]
                    }
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
              }
              f <- function(x,...) {
                sum(weight * x)
              }
              weight <- weight / sum(weight)
              p <- predict(x,newdata,w=id,filename=filename,...)
              calc(p,f,filename=filename,overwrite=TRUE)
            } else if (method == 'unweighted') {
              p <- predict(x,newdata,w=id,filename=filename,mean=TRUE,...)
              calc(p,mean,filename=filename,overwrite=TRUE)
            }
          }
)
