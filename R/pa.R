# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March 2024
# last update: March 2024
# Version 1.0
# Licence GPL v3
#-----------------



if (!isGeneric("pa")) {
  setGeneric("pa", function(x,y,id,opt,...)
    standardGeneric("pa"))
}



setMethod('pa', signature(x='SpatRaster',y='sdmModels'), 
          function(x,y,id,opt,...) {
            if (missing(opt)) opt <- 2
            
            if (missing(id)) {
              id <- getModelId(y,success = TRUE)
              if (length(id) != nlyr(x))  {
                if (nlyr(x) == 1) {
                  id <- 'ensemble'
                  cat('\nSince id is not specified, the threshold is obtained based on the ensemble of the models!\n')
                } else stop('The number of raster layers in x is not the same as the number of models y!')
              }
            } else if (is.character(id)) {
              if (any(c('ens','en','ensemble','ensmble','ensmbl') %in% tolower(id))) id <- 'ensemble'
            } else if (!is.numeric(id)) {
              .id <- getModelId(y,success = TRUE)
              if (is.logical(id) && length(id) == length(.id) && length(which(id)) == nlyr(x)) id <- .id[id]
              else if (length(.id) != nlyr(x)) stop('id is unknown; it should be either a numeric vector of modelIDs, or a single character: "ensemble"!')
              else id <- .id
            }
            #--------
            
            if (is.numeric(id)) {
              th <- getEvaluation(y,id=id,stat = 'threshold',opt=opt)[,2]
              
              if (length(th) == nlyr(x)) {
                pa <- x[[1]]
                pa <- ifel(pa >= th[1],1,0)
                if (nlyr(x) > 1) {
                  for (i in 2:length(th)) {
                    .pa <- x[[i]]
                    pa <- c(pa,ifel(.pa >= th[1],1,0))
                  }
                }
              } else stop('the length of thresholds extracted from the sdmModels object is not the same as the nlyr(x)...!')
            } else if (is.character(id)) {
              if (nlyr(x) == 1) {
                th <- evaluates(y@data,x)@threshold_based$threshold[opt]
                pa <- x
                pa <- ifel(pa >= th,1,0)
              } else stop('when id="ensemble", the nlyr(x) should be 1')
            }
            #-------------
            pa
          }
)

#--------------

setMethod('pa', signature(x='SpatRaster',y='numeric'), 
          function(x,y,id,opt,...) {
            if (missing(id)) id <- NULL
            
            if (missing(opt)) opt <- NULL
            
            if (is.matrix(y)) stop('y is matrix; should be a numeric vector!')
            
            if (!all(y >= 0 & y <= 1)) stop('threshold value(s) (provided in the y argument) should be within the range of [0, 1]')
            
            if (length(y) != nlyr(x)) stop('The length of y')
            #--------
            pa <- x[[1]]
            pa <- ifel(pa >= y[1],1,0)
            if (nlyr(x) > 1) {
              for (i in 2:length(y)) {
                .pa <- x[[i]]
                pa <- c(pa,ifel(.pa >= y[1],1,0))
              }
            }
            #-------------
            pa
          }
)







