# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Oct. 2024
# Version 1.0
# Licence GPL v3
#--------





# get Area of Applicability (assess which variables are outside of training data)
.getAOA <- function(d,x,vi=NULL) {
  vn <- d@sdmFormula@vars@numeric
  x2 <- x
  
  for (.n in vn$names) {
    .v <- vn[vn$names == .n,c('min','max')]
    x2[[.n]] <- ifel(x[[.n]] >= .v[1,1] & x[[.n]] <= .v[1,2],1,0)
  }
  
  if (!is.null(vi) && length(vi) == nlyr(x2)) {
    vi <- vi / sum(vi)
    weighted.mean(x2,vi)
  }
  else app(x2,'mean')
}
#-----------


if (!isGeneric("aoa")) {
  setGeneric("aoa", function(x,d,vi=NULL)
    standardGeneric("aoa"))
}


setMethod('aoa', signature(x='SpatRaster',d='sdmdata'), 
          function(x,d,vi=NULL)  {
            
            if (missing(vi)) vi <- NULL
            
            if (!is.null(vi) && length(vi) != nlyr(x)) stop('vi should have the length equal to the number of layers in x')
            
            if (!is.null(vi)) {
              if (!is.null(names(vi))) {
                if (all(names(vi) %in% names(x))) vi <- vi[names(x)]
                else {
                  warning('names of vi are not the same as the names layers in x; they are ignored and so it is assumed values in vi are in the same order as the layers in x!')
                }
              }
            }
            
            
            .getAOA(d,x,vi=vi)
          }
)


setMethod('aoa', signature(x='SpatRaster',d='sdmModels'), 
          function(x,d,vi=NULL)  {
            
            if (missing(vi)) vi <- NULL
            
            if (!is.null(vi) && length(vi) != nlyr(x)) stop('vi should have the length equal to the number of layers in x')
            
            if (!is.null(vi)) {
              if (!is.null(names(vi))) {
                if (all(names(vi) %in% names(x))) vi <- vi[names(x)]
                else {
                  warning('names of vi are not the same as the names layers in x; they are ignored and so it is assumed values in vi are in the same order as the layers in x!')
                }
              }
            }
              
              .getAOA(d@data,x,vi=vi)
          }
)

