\name{predict}


\docType{methods}

\alias{predict}
\alias{predict,sdmModels-method}
\alias{predict,.maxlikeModel-method}

\title{sdm model prediction}

\description{
Make a Raster or matrix object (depending on input dataset) with predictions from one or several fitted models in \code{sdmModels} object.
}
\usage{
\S4method{predict}{sdmModels}(object, newdata, filename="", id=NULL,species=NULL
          ,method=NULL,replication=NULL,run=NULL,mean=FALSE,
          overwrite=TRUE,parallelSetting, ...)
}


\arguments{
  \item{object}{sdmModels object}
  \item{newdata}{SpatRaster object, or data.frame}
  \item{filename}{character, output filename, if missing, a name starts with sdm_prediction will be generated}
  \item{id}{numeric (optional), specifies which model(s) should be used if the object contains several models; with NULL all models are considered}
  \item{species}{character (optional), specifies which species should be used if the object contains models for multiple species; with NULL all species are used}
  \item{method}{character, names of fitted models, e.g., glm, brt, etc.}
  
  \item{replication}{character (optional), specifies the names of replication method,if NULL, all available replications are considered}
  \item{run}{numeric (optional), works if replication with multiple runs are used}
  \item{mean}{logical, works if replication with multiple runs are used to fit the models, and specifies whether a mean should be calculated over all predictions of a replication method (e.g., bootstrapping) for each modelling method.}
  
  \item{overwrite}{logical, whether the filename should be overwriten it it does exist}
  \item{parallelSetting}{default is NULL; a list contains setting items for parallel processing. The items in parallel setting include: ncore, method, type, hosts, doParallel, fork, and strategy; see details for more information.}
  \item{...}{ additional arguments, as for \code{\link[terra]{writeRaster}}}
  }
\details{
predict uses the fitted models in the sdmModels object to generate the predictions given newdata. A SpatRaster object (if the newdata is Raster) or a data.frame (if newdata is data.frame) will be returned.

The predictions can be generated for some of models in the \code{sdmModels} object by specifying \code{id} (modelIDs) or explicitely specifying the names of \code{species}, or \code{method}, \code{replication} or \code{run} (replications ID). 

For each prediction, a name is assigned which is an abbreviation representing the names of species, method, replication method, and run (replication ID). If the output is a SpatRaster object, \code{\link[terra]{metags}} function can be used to get full names of raster layers.


For parallel processing, a list of items can be passed to \code{parallelSetting}, including:

\code{ncore}: defines the number of cores (it can also be specified outside of this list

\code{method}: defines the parallelising engine. Currently, three options are available including 'parallel', 'foreach', and 'future'. default is 'parallel'

\code{doParallel}: Optional, definition to register for a backend for parallel processing (needed when method='foreach'). It should be provided as an R expression like the following example:

expression(registerDoParallel(parallelSetting@cl))

The above example is based on the function available in doParallel package. Other packages can also be used to provide and register backend technologies (e.g., doMC)

\code{cluster}: Optional; in case a cluster is created and available (e.g., using cl <- parallel::makeCluster(2)), the cluster object can be introduced here to be used as the parallel processing engine, otherwise, it is handled by the sdm package.

\code{hosts}: Optional; To use remote machines/clusters in the parallel processing, a character vector with the addresses (names or IPs) of the accessible (on the network) remote clusters can be provided here to be registered and used in parallel processing (still under development so it may not work appropriately!)

\code{fork}: Logical, Available for non-windows operating system and specifies whether a fork solution should be used for the parallelisation. Default is TRUE for non-windows OS and FALSE for windows.

\code{strategy}: character (default='auto'), specifies the parallelisation strategy that can be either 'data' (split data across multiple parallel cores) or 'model' (predict for different models in parallel); if 'auto' is selected, it is decided by the function depending on the size of dataset and number of models.


NOTE: Only use parallelSetting when you deal with a big dataset or large number of models otherwise, it make the procedure slower rather than faster if the procedure is quick without parallelising!

}

\value{
 a SpatRaster object or data.frame
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}

\seealso{# }

\examples{
\dontrun{

file <- system.file("external/species.shp", package="sdm") # get the location of the species data

species <- vect(file) # read the shapefile

path <- system.file("external", package="sdm") # path to the folder contains the data

lst <- list.files(path=path,pattern='asc$',full.names = T) # list the name of the raster files 


# stack is a function in the raster package, to read/create a multi-layers raster dataset
preds <- rast(lst) # making a raster object

d <- sdmData(formula=Occurrence~., train=species, predictors=preds)

d

# fit the models (5 methods, and 10 replications using bootstrapping procedure):
m <- sdm(Occurrence~.,data=d,methods=c('rf','tree','fda','mars','svm'),
          replicatin='boot',n=10)
    
# predict for all the methods and replications:    
p1 <- predict(m, newdata=preds, filename='preds.tif')
plot(p1)

# predict for all the methods but take the mean over all replications for each replication method:
p2 <- predict(m, newdata=preds, filename='preds.img',mean=T)
plot(p2)

# for parallel processing
p3 <- predict(m, newdata=preds, filename='preds.tif',parallelSetting=list(ncore=2))

}
}

\keyword{spatial}
\keyword{data}
\keyword{model}
