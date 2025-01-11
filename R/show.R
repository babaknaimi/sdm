# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Jan. 2025
# Version 2.7
# Licence GPL v3

setMethod ('show' , 'sdmdata',
           function ( object ) {
             cat('class                                 :' , class(object), '\n')
             cat('===========================================================','\n')
             cat('number of species                     : ' , length(object@species.names) , '\n')
             cat('species names                         : ' , if (length(object@species.names) > 3) paste(c(object@species.names[1:3],'...'),collapse=', ') else paste(object@species.names,collapse=', ') , '\n')
             cat('number of features                    : ' , length(object@features.name), '\n')
             cat('feature names                         : ' , if (length(object@features.name) > 0) {
               if (length(object@features.name) > 3) paste(c(object@features.name[1:3],'...'),collapse=', ') else paste(object@features.name,collapse=', ')
             } else NA, '\n')
             if (length(object@factors) > 0) {
               cat('which feature is categorical (factor) : ' , 
                 if (length(object@factors) > 3) paste(c(object@factors[1:3],'...'),collapse=', ') else paste(object@factors,collapse=', ')
               , '\n') 
             }
             
             typ <- c()
             for (n in object@species.names) typ <- c(typ,object@species[[n]]@type)
             if (length(unique(typ)) == 1) typ <- unique(typ)
             hasTest <- 'test' %in% .getGroupNames(object,levels=TRUE)
             hasTrain <- 'train' %in% .getGroupNames(object,levels=TRUE)
             cat('type                                  : ', if (length(typ) > 3) paste(c(typ[1:3],'...'),collapse=', ') else paste(typ,collapse=', '), '\n')
             cat('has independent test data?             : ' , hasTest, '\n')
             if (hasTrain) {
               cat('number of records                     : ', if (hasTest) paste('train-> ',length(.getGroupIndex(object,'train')),"; ",'test-> ',length(.getGroupIndex(object,'test')),sep='') else length(.getSpeciesIndex(object)),'\n')
             } else
             cat('number of records                     : ', if (is.null(object@features)) '0' else nrow(object@features),'\n')
             cat('has Coordinates?                      : ' , !is.null(object@info) && !is.null(object@info@coords)  , '\n')
           }
)
#-----------
setMethod ('show' , '.sdmCorSetting',
           function (object) {
             cv.n <- 1
             r.n <- 1
             r <- 1
             
             cat('class                                 :' , class(object), '\n')
             cat('========================================================','\n')
             cat('modelling methods                     : ' , paste(object@methods,collapse=', '), '\n')
             sn <- object@sdmFormula@vars@species
             vn <- object@sdmFormula@vars@names
             vn <- .excludeVector(vn,sn)
             cat('species names                         : ' , if (length(sn) > 4) paste0(paste(c(sn[1:4],'...'),collapse=', '),' (',length(sn),' species)') else paste(sn,collapse=', '), '\n')
             cat('predictor names                       : ' , if (length(vn) > 4) paste0(paste(c(vn[1:4],'...'),collapse=', '),' (',length(vn),' variables)') else paste(vn,collapse=', '), '\n')
             #cat('feature types                         : ' , paste(unique(unlist(lapply(object@featureFrame@feature.types,function(x) x@type))),collapse=', '), '\n')
             if (!is.null(object@replicate)) {
               r.n <- length(object@replicate)
               r <- object@n.replicates
               cat('replicate.methods (data partitioning) : ' , paste(object@replicate,collapse=','), '\n')
               if ('subsampleing' %in% object@replicate)
                 cat('test percentage                       : ' , object@test.percentage, '\n')
               cat('number of replicates (each method)    : ' , r, '\n')
               if ("cross-validation" %in% object@replicate) {
                 cat('n.folds in cross-validation           : ' ,object@cv.folds , '\n')
               }
             }
             n <- r * (r.n-1+object@cv.folds)
             cat('------------------------------------------\n')
             cat('number of runs                        : ' , paste(n,' for each model, ',n*length(object@methods),' in total... (per species)' ,sep=''), '\n')
             
           }
)
#---------------------

setMethod ('show' , 'sdmModels',
           function (object) {
             if (!.sdmOptions$getOption('sdmLoaded')) .addMethods()
             mi <- object@run.info
             for (i in c(2,3)) mi[,i] <- as.character(mi[,i])
             sp <- unique(mi$species)
             mo <- unique(mi$method)
             n <- 1
             
             cat('class                                 :' , class(object), '\n')
             cat('========================================================','\n')
             cat('number of species                     : ' , length(sp), '\n')
             cat('number of modelling methods           : ' , length(mo), '\n')
             cat('names of modelling methods            : ' , paste(mo,collapse=', '), '\n')
             
             if (!is.na(mi[1,4])) {
               r.n <- length(unique(mi$replication))
               n <- length(object@replicates[[1]])
               
               
               cat('replicate.methods (data partitioning) : ' , paste(unique(mi$replication),collapse=','), '\n')
               cat('number of replicates (each method)    : ' , object@setting@n.replicates, '\n')
               cat('toral number of replicates per model  : ' , paste(n,'(per species)'), '\n')
               if ('subsampling' %in% unique(mi$replication))
                 cat('test percentage (in subsampling)      : ' , object@setting@test.percentage, '\n')
               
             }
             cat('------------------------------------------\n')
             cat('model run success percentage (per species)  :\n')
             cat('------------------------------------------\n')
             
             p2 <- function(x) {
               o <- rep(NA,length(sp))
               for (i in seq_along(sp)) {
                 w1 <- mi$species == sp[i]
                 w2 <- mi$method == x
                 w <- mi$success[w1 & w2]
                 o[i] <- (length(which(w))/length(w))*100
               }
               o
             }
             
             a <- ''
             a2 <- c()
             
             lb <- max(unlist(lapply(mo,function(x) length(strsplit(x,'')[[1]]))))
             
             if (lb < 10) lb <- 10
             else lb <- lb+2
             
             mAdd <- 0
             if (lb > 10) mAdd <- lb-10
             
             for (i in sp) {
               aa <- length(unlist(strsplit(i,'')))
               la <- (15 + mAdd)
               if (aa >= (15 + mAdd)) la <- aa + 3 + mAdd
               a <- paste(a,i,paste(rep(' ',la - aa),collapse=''),collapse='')
               a2 <- c(a2,la)
             }
             a2 <- a2 - 3
             cat(paste('method        ',a,collapse='|'),'\n')
             cat(paste(rep('-',length(unlist(strsplit(a,'')))+5),collapse=''),'\n')
             for (i in seq_along(mo)) {
               a3 <- c()
               a4 <- ''
               p2.mo <- p2(mo[i])
               if (length(a2) > 1) {
                 for (j in 1:(length(a2)-1)) a3 <- c(a3,paste(paste(rep(' ',a2[j]/2),collapse=''),'|',paste(rep(' ',a2[j+1]/2),collapse=''),collapse='')) 
               } else {
                 a3 <- paste(rep(' ',a2),collapse='') 
               }
               
               a3 <- c(paste(rep(' ',a2[1]/2-2),collapse=''),a3)
               for (j in 1:length(a2)) a4 <- paste(c(a4,a3[j],p2.mo[j]),collapse='')
               
               cat(paste(mo[i],paste(rep(' ',lb - length(unlist(strsplit(mo[i],'')))),collapse=''),' : ',paste(rep(' ',3),collapse=''),a4,sep='')  , '  %\n')
             }
             #wtest <- colnames(mi)[9:7][which(as.matrix(mi[1,c(9,8,7)]))[1]]
             wtest <- colnames(mi)[9:7][which(apply(mi[,c(9,8,7)],2,any))[1]]
             
             cat('\n###################################################################\n')
             if (n > 1) {
               if (wtest == 'test.dep')
                 cat('model Mean performance (per species), using test dataset (generated using partitioning):\n')
               else if (wtest == 'test.indep')
                 cat('model Mean performance (per species), using independent test dataset:\n')
               else
                 cat('model Mean performance (per species), using training test dataset:\n')
             } else {
               if (wtest == 'test.dep')
                 cat('model performance (per species), using test dataset (generated using partitioning):\n')
               else if (wtest == 'test.indep')
                 cat('model performance (per species), using independent test dataset:\n')
               else
                 cat('model performance (per species), using training test dataset:\n')
             } 
             
             cat('-------------------------------------------------------------------------------\n')
             p <- function(x,sp,s1=c('AUC','COR','Deviance'),s2='TSS') {
               a <- c()
               w1 <- mi$species == sp
               w2 <- mi$method == x
               id <- mi$modelID[w1 & w2 & mi$success]
               if (length(id) > 0) {
                 o <- ._getPerformance(object,id,wtest = NULL,s1=s1,s2=s2,opt = 2)
                 if (is.null(s2)) {
                   a <- mean(sapply(o,function(x) x[[s1[1]]]),na.rm=TRUE)
                   a <- c(a,mean(sapply(o,function(x) x[[s1[2]]]),na.rm=TRUE))
                   a <- c(a,mean(sapply(o,function(x) x[[s1[3]]]),na.rm=TRUE))
                   a <- c(a,mean(sapply(o,function(x) x[[s1[4]]]),na.rm=TRUE))
                 } else {
                   a <- mean(sapply(o,function(x) x[[s1[1]]]),na.rm=TRUE)
                   a <- c(a,mean(sapply(o,function(x) x[[s1[2]]]),na.rm=TRUE))
                   a <- c(a,mean(sapply(o,function(x) x[[s2[1]]]),na.rm=TRUE))
                   a <- c(a,mean(sapply(o,function(x) x[[s1[3]]]),na.rm=TRUE))
                 }
                 
                 a <- round(a,2)
                 a <- as.character(a)
               } else a <- as.character(c(NA,NA,NA,NA))
               b <- c()
               for (i in a) b <- c(b,paste(i,paste(rep(' ',7 - length(unlist(strsplit(i,'')))),collapse=''),collapse=''))
               paste(b,collapse='|     ')
             }
             
             if (object@setting@distribution[1] == 'binomial') {
               for (spp in sp) {
                 
                 cat('\n ## species   : ',spp,'\n')
                 cat('=========================\n\n')
                 if (lb > 10) cat(paste('methods',paste(rep(' ',(lb/3+1)),collapse=''),' : ',paste(rep(' ',lb/3),collapse=''),sep='')  , paste(c('AUC','COR','TSS','Deviance'),collapse='     |     '), '\n')
                 else cat(paste('methods',paste(rep(' ',3),collapse=''),' : ',paste(rep(' ',3),collapse=''),sep='')  , paste(c('AUC','COR','TSS','Deviance'),collapse='     |     '), '\n')
                 cat('-------------------------------------------------------------------------\n')
                 for (i in seq_along(mo)) {
                   cat(paste(mo[i],paste(rep(' ',lb - length(unlist(strsplit(mo[i],'')))),collapse=''),' : ',paste(rep(' ',3),collapse=''),sep='')  , p(mo[i],spp), '\n')
                 }
               }
             } else {
               for (spp in sp) {
                 
                 cat('\n ## species   : ',spp,'\n')
                 cat('=========================\n\n')
                 if (lb > 10) cat(paste('methods',paste(rep(' ',(lb/3+1)),collapse=''),' : ',paste(rep(' ',lb/3),collapse=''),sep='')  , paste(c('RMSE','COR','MAE','Deviance'),collapse='     |     '), '\n')
                 else cat(paste('methods',paste(rep(' ',3),collapse=''),' : ',paste(rep(' ',3),collapse=''),sep='')  , paste(c('RMSE','COR','MAE','Deviance'),collapse='     |     '), '\n')
                 cat('-------------------------------------------------------------------------\n')
                 for (i in seq_along(mo)) {
                   cat(paste(mo[i],paste(rep(' ',lb - length(unlist(strsplit(mo[i],'')))),collapse=''),' : ',paste(rep(' ',3),collapse=''),sep='')  , p(mo[i],spp,s1=c('RMSE','COR','MAE','Deviance'),s2=NULL), '\n')
                 }
               }
             }
             
           }
)

#----------

setMethod ('show' , '.varImportanceList',
           function (object) {
             cat('Relative Variable Importance List \n')
             cat('=============================================================','\n')
             cat('method              : Permutation based on two metrics (Pearson Correlation and AUC)\n')
             cat('number of variables : ' , length(object@variables), '\n')
             cat('variable names      : ' , if (length(object@variables) > 3) paste(c(object@variables,'...'),collapse=', ') else paste(object@variables,collapse=', '), '\n')
             cat('number of models    : ' , length(object@varImportanceList), '\n')
             cat('=============================================================','\n')
             cat('Summary of relative variable importance \n')
             cat('----------------------------------------------','\n')
             cat('Based on Correlation metric: \n')
             cat('----------------------------------------------','\n')
             for (n in object@variables) {
               .ns <- unlist(strsplit(n,''))
               w <- which(object@varImportanceMean$corTest$variables == n)
               v <- round(object@varImportanceMean$corTest$corTest[w]*50)
               ci1 <- round(object@varImportanceMean$corTest$lower[w]*50)
               ci2 <- round(object@varImportanceMean$corTest$upper[w]*50)
               vp <- round(object@varImportanceMean$corTest$corTest[w]*100,1)
               
               if (length(.ns) >= 20) {
                 xx <- c(.ns[1:20],': ',rep('*',v))
                 xx[21+ci1] <- '['
                 xx[21+ci2] <- ']'
                 if (any(is.na(xx))) xx[is.na(xx)] <- '-'
                 xx <- c(xx,' (',vp,' %)')
                 cat(paste(xx,collapse=''),'\n')
               } else {
                 xx <- c(.ns,rep(' ',20 - length(.ns)), ': ',rep('*',v))
                 xx[21+ci1] <- '['
                 xx[21+ci2] <- ']'
                 if (any(is.na(xx))) xx[is.na(xx)] <- '-'
                 xx <- c(xx,' (',vp,' %)')
                 cat(paste(xx,collapse=''),'\n')
               }
             }
             
             cat('=============================================================','\n')
             cat('Based on AUC metric: \n')
             cat('----------------------------------------------','\n')
             for (n in object@variables) {
               .ns <- unlist(strsplit(n,''))
               w <- which(object@varImportanceMean$AUCtest$variables == n)
               v <- round(object@varImportanceMean$AUCtest$AUCtest[w]*50)
               ci1 <- round(object@varImportanceMean$AUCtest$lower[w]*50)
               ci2 <- round(object@varImportanceMean$AUCtest$upper[w]*50)
               vp <- round(object@varImportanceMean$AUCtest$AUCtest[w]*100,1)
               
               if (length(.ns) >= 20) {
                 xx <- c(.ns[1:20],': ',rep('*',v))
                 xx[21+ci1] <- '['
                 xx[21+ci2] <- ']'
                 if (any(is.na(xx))) xx[is.na(xx)] <- '-'
                 xx <- c(xx,' (',vp,' %)')
                 cat(paste(xx,collapse=''),'\n')
               } else {
                 xx <- c(.ns,rep(' ',20 - length(.ns)), ': ',rep('*',v))
                 xx[21+ci1] <- '['
                 xx[21+ci2] <- ']'
                 if (any(is.na(xx))) xx[is.na(xx)] <- '-'
                 xx <- c(xx,' (',vp,' %)')
                 cat(paste(xx,collapse=''),'\n')
               }
             }
             cat('=============================================================','\n')
             
           }
)
####----------

setMethod ('show' , '.varImportance',
           function (object) {
             cat('Relative Variable Importance \n')
             cat('=============================================================','\n')
             cat('method              : Permutation based on two metrics (Pearson Correlation and AUC)\n')
             cat('number of variables : ' , length(object@variables), '\n')
             cat('variable names      : ' , if (length(object@variables) > 3) paste(c(object@variables,'...'),collapse=', ') else paste(object@variables,collapse=', '), '\n')
             cat('=============================================================','\n')
             cat('Relative variable importance \n')
             cat('----------------------------------------------','\n')
             cat('Based on Correlation metric: \n')
             cat('----------------------------------------------','\n')
             for (n in object@variables) {
               .ns <- unlist(strsplit(n,''))
               w <- which(object@varImportance$variables == n)
               v <- round(object@varImportance$corTest[w]*50)
               vp <- round(object@varImportance$corTest[w]*100,1)
               
               if (length(.ns) >= 20) {
                 xx <- c(.ns[1:20],': ',rep('*',v))
                 xx <- c(xx,' (',vp,' %)')
                 cat(paste(xx,collapse=''),'\n')
               } else {
                 xx <- c(.ns,rep(' ',20 - length(.ns)), ': ',rep('*',v))
                 xx <- c(xx,' (',vp,' %)')
                 cat(paste(xx,collapse=''),'\n')
               }
             }
             
             cat('=============================================================','\n')
             cat('Based on AUC metric: \n')
             cat('----------------------------------------------','\n')
             for (n in object@variables) {
               .ns <- unlist(strsplit(n,''))
               w <- which(object@varImportance$variables == n)
               v <- round(object@varImportance$AUCtest[w]*50)
               vp <- round(object@varImportance$AUCtest[w]*100,1)
               
               if (length(.ns) >= 20) {
                 xx <- c(.ns[1:20],': ',rep('*',v))
                 xx <- c(xx,' (',vp,' %)')
                 cat(paste(xx,collapse=''),'\n')
               } else {
                 xx <- c(.ns,rep(' ',20 - length(.ns)), ': ',rep('*',v))
                 xx <- c(xx,' (',vp,' %)')
                 cat(paste(xx,collapse=''),'\n')
               }
             }
             cat('=============================================================','\n')
             
           }
)
#------

setMethod ('show' , '.pcaObject',
           function ( object ) {
             cat('class                                 :' , class(object), '\n')
             cat('===========================================================','\n')
             cat('number of components                  : ' , ncol(object@pcaObject$scores) , '\n')
             .vi <- object@pcaObject$sdev * object@pcaObject$sdev
             .vi <- .vi / sum(.vi)
             names(.vi) <- colnames(object@pcaObject$scores)
             
             .load <- object@pcaObject$loadings
             
             
             cat('----------------------------------------------','\n')
             cat('Proportion of Variance: \n\n')
             #cat('Proportion of Variance                : ' , if (length(.vi) > 4) paste(c(round(.vi[1:4],3),'...'),collapse=', ') else paste(round(.vi,3),collapse=', ') , '\n')
             for (n in names(.vi)) {
               .ns <- unlist(strsplit(n,''))
               v <- round(.vi[n]*50)
               vp <- round(.vi[n]*100,1) 
               xx <- c(.ns,rep(' ',10 - length(.ns)), ': ',rep('*',v))
               xx <- c(xx,' (',vp,' %)')
               cat(paste(xx,collapse=''),'\n')
             }
             cat('-----------------------------------------------','\n')
             
             .vic <- cumsum(.vi)
             .nc <- 1 - (1 / length(.vi))
             .nv <- which(.vic >= .nc)[1]
             .nc <- round(.nc,2)*100
             
             cat("nr of comp. explain more than one variable's worth of information (1 - [1 / N] =>",.nc,'%):' , .nv, '\n')
             cat('-----------------------------------------------','\n')
             cat('loadings:\n\n')
             if (.nv < 3) {
               .l1 <- abs(.load[,1])
               .l2 <- abs(.load[,2])
               cat(' -- The most important variables contributing to the first 2 components:\n')
               
               cat('            ---> ' , paste(c(names(.l1)[which.max(.l1)],names(.l2)[which.max(.l2)]),collapse=', ') , '\n')
             } else {
               .l1 <- abs(.load[,1])
               .l2 <- abs(.load[,2])
               .l3 <- abs(.load[,3])
               
               cat(' -- The most important variables contributing to the first 3 components:\n')
               cat('            ---> ' , paste(c(names(.l1)[which.max(.l1)],names(.l2)[which.max(.l2)],names(.l3)[which.max(.l3)]),collapse=', ') , '\n')
             }
             #----
             
             for (i in 1:.nv) {
               #cat('------------ \n')
               cat(paste0(' -- Important variables contributing to PC',i),':\n')
               .l1 <- abs(.load[,i])
               .l1 <- sort(.l1,decreasing = TRUE)[1]
               cat('            ---> ' , paste(paste0(names(.l1),' (',round(.l1,3),')'),collapse=', ') , '\n')
             }
             cat('----------------------------------------------','\n')
           }
)
#-----------



setMethod ('show' , 'sdmEvaluate',
           function ( object ) {
             cat('class                                 :' , class(object), '\n')
             cat('===========================================================','\n')
             cat('number of records                     : ' , length(object@observed) , '\n')
             cat('----------------------------------------------','\n')
             cat('Evaluation metrics (threshold_Independent)...                          \n')
             cat('----------------------------------------------','\n')
             print(object@statistics)
             cat('----------------------------------------------','\n')
             if (!is.null(object@threshold_based)) {
               cat('Threshold_based metrics (part of the table is printed here...)\n')
               cat('----------------------------------------------','\n')
               print(object@threshold_based[c(1,2,4,13,14),c(1,2,3,4,5,6,7)])
               #for ( i in c())
               cat('=============================================================','\n')
             }
             
           }
)
#-----------
# paste(colnames(e@threshold_based[1,c(2,3,4,5,6,7)]),collapse = ' | ')
# paste(e@threshold_based[2,1],'   |  ',paste(round(as.numeric(e@threshold_based[2,c(2,3,4,5,6,7)]),3),collapse = '   | '),'| ...')
# 
# .w <- length(strsplit(e@threshold_based[2,1],'')[[1]])
# paste(c(rep(' ',14-.w),'|'),collapse = '')
# c("Criteria  |  thr.    |  sen.   |  spe.   |   TSS   |   MCC   |   F1  | ...")

