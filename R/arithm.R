# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2018
# Last update: March 2018
# Version 1.0
# Licence GPL v3
##############################################

# compare two vectors to check if they have the same set of elements (can have different order)
.identicalContent <- function(x,y) {
  all(x %in% y) & all(y %in% x)
}

.combineSdmData <- function(d1,d2) {
  if (.identicalContent(d1@features.name,d2@features.name)) {
    stop('the datasets have different set of predictor variables')
  }
  
  if (.identicalContent(d1@species.names,d2@species.names)) {
    .sameSp <- TRUE
  } else .sameSp <- FALSE
  
  if (is.null(d1@info)) .d1Info <- NULL
  if (is.null(d2@info)) .d2Info <- NULL
  
  if (.identicalContent(names(d1@groups),names(d1@groups))) .sameGrp <- TRUE
  else .sameGrp <- FALSE
  
  if (nrow(d1@features) != nrow(d2@features))  .sameRecords <- FALSE
  else {
    .sameRecords <- all(sapply(rownames(d1@features),function(x) {
      identical(d1@features[x,d1@features.name],d2@features[x,d1@features.name])
    }))
  }
  
  
}

# 
# .combineSdmSettings <- function(d1,d2) {
#   
# }


.combineModels <- function(m1,m2) {
  smo <- new('sdmModels',setting=m1@setting,data=m1@data)
  .sameSp <- .identicalContent(names(m1@models),names(m2@models))
  .sameModels <- .identicalContent(names(m1@models[[1]]),names(m2@models[[1]]))
  if (!.sameSp && !.sameModels) stop('the sdm objects are not compatible (both different species and different models are TRUE)')
  
  .startID <- max(m1@run.info$modelID) + 1
  mi2 <- m2@run.info
  .newID <- .startID:(.startID+nrow(mi2)-1)
  for (i in 1:nrow(mi2)) {
    .m <- list(m2@models[[as.character(mi2$species[i])]][[as.character(mi2$method[i])]][[as.character(mi2$modelID[i])]])
    names(.m) <- as.character(.newID[i])
    m1@models[[as.character(mi2$species[i])]][[as.character(mi2$method[i])]] <- c(m1@models[[as.character(mi2$species[i])]][[as.character(mi2$method[i])]],.m)
    #m1@models[[as.character(mi2$species[i])]][[as.character(mi2$method[i])]][[as.character(.newID[i])]] <- .m
  }
  
  mi2$modelID <- .newID
  rownames(mi2) <-as.character(mi2$modelID)
  m1@run.info <- rbind(m1@run.info,mi2)
  
  # if (!.sameModels) {
  #   
  # } else {
  #   
  # }
  m1
}



setMethod("+", signature(e1='sdmModels', e2='sdmModels'),
          function(e1, e2){ 
            .combineModels(e1, e2)
          }
)
