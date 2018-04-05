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
\S4method{predict}{sdmModels}(object, newdata, filename="", w=NULL,species=NULL
          ,method=NULL,replication=NULL,run=NULL,mean=FALSE,control=NULL,
          overwrite=TRUE,nc=1,obj.size=1, ...)
}


\arguments{
  \item{object}{sdmModels object}
  \item{newdata}{Raster* object, or data.frame}
  \item{filename}{character, output file name, if missing, a name starts with sdm_prediction will be generated}
  \item{w}{numeric, specifies which model(s) should be used if the object contains several models; with NULL all models are used}
  \item{species}{character, (optional), specifies which species should be used if the object contains models for multiple species; with NULL all species are used}
  \item{method}{character, names of fitted models, e.g., glm, brt, etc.}
  
  \item{replication}{character, specifies the names of replication method,if NULL, all available replications are considered}
  \item{run}{numeric, works if replication with multiple runs are used}
  \item{mean}{logical, works if replication with multiple runs are used to fit the models, and specifies whether a mean should be calculated over all predictions of a replication method (e.g., bootstrapping) for each modelling method.}
  
  \item{control}{not implemented yet!}
  
  \item{overwrite}{logical, whether the filename should be overwriten it it does exist}
  \item{nc}{number of cores for parallel running of the function}
  \item{obj.size}{the size of object can be kept in memory (default=1 Giga byte). Depending on the available memory, this value can be changed}
  \item{...}{ additional arguments, as for \code{\link[raster]{writeRaster}}}
  }
\details{
predict uses the fitted models in the sdmModel to generate the prediction given newdata. A raster (if newdata is Raster object) or data.frame (if newdata is data.frame) will be created.

The predictions can be generated for a specific set of models in the input \code{sdmModels} by determining either or a combination of the name of 

For each prediction, a name is assigned which is kind of abbreviation or codding that tells which species, which method, which replication method, and which run is the prediction for. If the output is a Raster object, \code{\link[raster]{setZ}} function can be used to get a full name of each layer.

}

\value{
 a Raster object or data.frame
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}

\seealso{# }

\examples{
\dontrun{

file <- system.file("external/species.shp", package="sdm") # get the location of the species data

species <- shapefile(file) # read the shapefile

path <- system.file("external", package="sdm") # path to the folder contains the data

lst <- list.files(path=path,pattern='asc$',full.names = T) # list the name of the raster files 


# stack is a function in the raster package, to read/create a multi-layers raster dataset
preds <- stack(lst) # making a raster object

d <- sdmData(formula=Occurrence~., train=species, predictors=preds)

d

# fit the models (5 methods, and 10 replications using bootstrapping procedure):
m <- sdm(Occurrence~.,data=d,methods=c('rf','tree','fda','mars','svm'),
          replicatin='boot',n=10)
    
# predict for all the methods and replications:    
p1 <- predict(m, newdata=preds, filename='preds.img')
plot(p1)

# predict for all the methods but take the mean over all replications for each replication method:
p2 <- predict(m, newdata=preds, filename='preds.img',mean=T)
plot(p2)

}
}

\keyword{spatial,species,data}
