# Author: Babak Naimi, naimi.b@gmail.com
# Date :  October 2024
# last update: October 2024
# Version 1.0
# Licence GPL v3
#-----------------



if (!isGeneric("getThreshold")) {
  setGeneric("getThreshold", function(x,id,opt,...)
    standardGeneric("getThreshold"))
}


if (!isGeneric("threshold")) {
  setGeneric("threshold", function(x,id,opt,...)
    standardGeneric("threshold"))
}

#----
setMethod('getThreshold', signature(x='sdmModels',id='numeric'), 
          function(x,id,opt,...) {
            if (missing(opt)) opt <- 2
            .id <- getModelId(x,success = TRUE)
            
            if (!any(id %in% .id)) stop('None of the specified modelIDs in "id" are available in the sdmModels object (x)!')
            
            if (!all(id %in% .id)) {
              id <- id[id %in% .id]
              warning('Some of the modelIDs specified in "id" are not available in the sdmModels object, and they are ignored!')
            }
            #-------
            
            th <- getEvaluation(x,id=id,stat = 'threshold',opt=opt)[,2]
            th
          }
)
#-----
setMethod('getThreshold', signature(x='sdmModels',id='character'), 
          function(x,id,opt,species=NULL,...) {
            if (missing(opt)) opt <- 2
            
            if (missing(species)) species <- NULL
            
            
            if (any(c('ens','en','ensemble','ensmble','ensmbl','ensembl','ensembel') %in% tolower(id))) id <- 'ensemble'
            
            if (id != 'ensemble') stop('id is not recognised: should be either a numeric vector with modelIDs or id = "ensemble"!')
            
            
            
            n <- x@data@species.names
            
            if (length(n) > 1) {
              if (is.null(species)) {
                warning('multiple species are available in the model and the species is not specified in the "species" argument. The first species is considered!')
                n <- n[1]
              } else {
                if (length(species) > 1) stop('only one species should be specified in the "species" argument!')
                else {
                  if (is.numeric(species)) {
                    if (length(n) < species) stop('The specified species is not available (species argument)!')
                    else n <- n[species]
                  } else if (is.character(species)) {
                    w <- which(n == species)
                    if (length(w) == 1) n <- n[w]
                    else stop('The specified species (in the "species" argument) is not available!')
                  } else stop('"species" should be a character or numeric!')
                }
              }
            }
            #---------
            if (!x@data@species[[n]]@type %in% c("Presence-Background","Presence-Absence")) stop('species data type is not Presence-Absence or Presence-Background!')
            
            en <- as.data.frame(x@data,sp=n)
            
            obs <- en[,n]
            
            en <- ensemble(x, en,...)
            if (length(x@data@species.names) > 1) {
              w <- which(grepl(n,colnames(en)))
              if (length(w) == 1) en <- en[,w]
              else {
                w <- which(x@data@species.names == n)
                en <- en[,w]
              }
            } else en <- en[,1]
            
            e <- evaluates(obs,en)
            
            
            e@threshold_based$threshold[opt]
          }
)
#-----------
####################
setMethod('threshold', signature(x='sdmModels',id='numeric'), 
          function(x,id,opt,...) {
            if (missing(opt)) opt <- 2
            
            getThreshold(x,id,opt,...)
            
          }
)
setMethod('threshold', signature(x='sdmModels',id='character'), 
          function(x,id,opt,species=NULL,...) {
            if (missing(opt)) opt <- 2
            
            if (missing(species)) species <- NULL
            
            
            if (any(c('ens','en','ensemble','ensmble','ensmbl','ensembl','ensembel') %in% tolower(id))) id <- 'ensemble'
            
            if (id != 'ensemble') stop('id is not recognised: should be either a numeric vector with modelIDs or id = "ensemble"!')
            
            getThreshold(x,id='ensemble',opt=opt,species=species,...)
            
          }
)