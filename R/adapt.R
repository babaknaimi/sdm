# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  July 2024
# Version 1.1
# Licence GPL v3
#--------



if (!isGeneric("sdmAdapt")) {
  setGeneric("sdmAdapt", function(x)
    standardGeneric("sdmAdapt"))
}


setMethod('sdmAdapt', signature(x='sdmdata'),
          function(x) {
            .ty <- x@species[[1]]@type
            df <- x[1:nrow(x@features),][,-1]
            x <- sdmData(x@sdmFormula@formula,df)
            x@species[[1]]@type <- .ty
            x
          }
)



setMethod('sdmAdapt', signature(x='sdmModels'),
          function(x) {
            .ty <- x@data@species[[1]]@type
            df <- x@data[1:nrow(x@data@features),][,-1]
            d <- sdmData(x@data@sdmFormula@formula,df)
            d@species[[1]]@type <- .ty
            #-------------
            .s <- sdmSetting(x@setting@sdmFormula@formula,methods = x@setting@methods,data = d,
                             interaction.depth = x@setting@interaction.depth,n = x@setting@n.replicates,
                             replication = x@setting@replicate,cv.folds = x@setting@cv.folds,
                             test.percent = x@setting@test.percentage,bg = x@setting@pseudo.absence.methods,
                             bg.n = x@setting@n.pseudo.absence,var.importance = x@setting@varImportance.methods,
                             response.curve = x@setting@response.curve,var.selection = x@setting@var.selection,
                             modelSettings = x@setting@modelSettings,seed = x@setting@seed,
                             parallelSetting = x@setting@parallelSettings)
            x@data <- d
            x@setting <- .s
            #----------
            # adapt the evaluation object in the model:
            for (sp in names(x@models)) {
              for (m in names(x@models[[sp]])) {
                for (id in names(x@models[[sp]][[m]])) {
                  ev <- x@models[[sp]][[m]][[id]]@evaluation
                  if (!is.null(ev)) {
                    for (wt in names(ev)) {
                      .e <- evaluates(ev[[wt]]@observed,ev[[wt]]@predicted)
                      x@models[[sp]][[m]][[id]]@evaluation[[wt]] <- .e
                    }
                  }
                }
              }
            }
            #-----
            
            x
          }
)

