\name{ensemble}


\docType{methods}

\alias{ensemble}
\alias{ensemble,sdmModels-method}
\alias{ensemble,sdmModels,Raster-method}

\title{Ensemble Forecasting of SDMs}

\description{
Make a Raster object with a weighted averaging over all predictions from several fitted model in a sdmModel object.
}
\usage{
\S4method{ensemble}{sdmModels}(x, newdata, filename="",setting,...)
}


\arguments{
  \item{x}{a sdmModels object}
  \item{newdata}{Raster* object or data.frame}
  \item{filename}{character, output file name}
  \item{setting}{list, contains the parameters that are used in the ensemble procedure; see details}
  \item{...}{additional arguments passed to the predict function}
  
  }
\details{
ensemble function uses the fitted models in an \code{sdmModels} object to generate an ensemble/consensus of predictions by individual models. Several methods do exist for this procedure, that are (or will be) implemented in this function, and can be defined in the method argument.
A list can be introduced in the \code{setting} argument in which several parameters can be set including:

- method: specify which ensemble method should be used. Currently, 'unweighted' (unweighted averaging), and 'weighted' (weighted averaging) are implemented, but more methods will be added.

- stat: if the method='weighted' is used, this specify which evaluation statistics can be used as weight in the weighted averaging procedure. Alternatively, one may directly introduce weights (see the next argument)

- weights: an optional numeric vector (with a length equal to the models that are successfully fitted) to specify the weights for weighted averaging procedure (if the method='weighted' is specified)

- id: specify the model IDs that should be considered in the ensemble procedure. If missing, all the models that are successfully fitted are considered.

- wtest: specify which test dataset ("training","test.dep","test.indep") should be used to extract the statistic (stat) values as weights (if a relevant method is specified)


}

\value{
 - a Raster object if \code{predictors} is a Raster object
 
 - a numeric vector if \code{predictors} is a data.frame object
}

\references{
#
}

\author{Babak Naimi \email{naimi.b@gmail.com}

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
    
# ensemble using weighted averaging based on AUC statistic:    
p1 <- ensemble(m, newdata=preds, filename='ens.img',setting=list(method='weighted',stat='AUC'))
plot(p1)

# ensemble using weighted averaging based on TSS statistic
# and optimum threshold critesion 2 (i.e., Max(spe+sen)) :    
p2 <- ensemble(m, newdata=preds, filename='ens2.img',setting=list(method='weighted',
                                                                  stat='TSS',opt=2))
plot(p2)

}


}

\keyword{spatial,species,data}
