\name{threshold}
\alias{getThreshold}
\alias{threshold,sdmModels,numeric-method}
\alias{threshold,sdmModels,character-method}
\alias{getThreshold,sdmModels,numeric-method}
\alias{getThreshold,sdmModels,character-method}

\title{Identifying/extracting the best threshold to convert probabilities to Presence/Absences}


\description{
To transform the values of predicted probability of occurrence (habitat suitability) to presence-absences, a threshold is needed. To identify the best threshold, several (15+) optimisation criteria are supported by the package which is calculated for each model in the sdmModels object. To extract the best threshold from each model (given its modelID specified in the \code{id} argument), these functions can be used. It is also possible to specify \code{id = "ensemble"} to identify the best threshold for the ensemble of models.
}

\usage{
threshold(x,id,opt,species=NULL,...)
getThreshold(x,id,opt,species=NULL,...)
}

\arguments{
  \item{x}{an sdmModels object (output of the sdmData function)}
  \item{id}{can be either a numeric vector specifying the modelIDs corresponding to SDMs in \code{x}, or a single character "ensemble" to ask the function to find the best threshold for the ensemble of models}
  \item{opt}{specifies the optimisation criterion based on which a threshold is identified; default is opt=2 (see \code{\link{evaluates}} for more details) }
  \item{species}{default=NULL; if the models for multiple species are available in \code{x}, \code{species} specifies the name of species or a single number specifies which species from the available species should be considered}
  \item{...}{if \code{id = "ensemble"}, additional arguments are passed to the ensemble function; for example, \code{setting} list can be provided for the ensemble function!}
  }

\details{
Both the \code{threshold} and \code{getThreshold} functions are the same

}

\value{
a numeric value
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881


}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org/}
}





\examples{
\dontrun{

# let's first fit a set of models and generate prediction and ensemble maps:

# get the path to the species data
file <- system.file("external/sp1.shp", package="sdm") 

sp <- vect(file) # read the species records

path <- system.file("external", package="sdm") # path to the folder contains the data

lst <- list.files(path=path,pattern='asc$',full.names = T) # list of predictor filenames

preds <- rast(lst) # making a SpatRaster object (predictors)

d <- sdmData(formula=Occurrence~., train=sp, predictors=preds)

d

# fit two models:
m <- sdm(Occurrence~.,data=d,methods=c('glmp','brt'),replication='boot',n=1)

m 

threshold(m,id=1,opt=1) # get threshold for the first model 

threshold(m,id=1:2,opt=1) # get thresholds for the first and second models

threshold(m,id=1:2,opt=2) # get thresholds for the first and second models

threshold(m,id="ensemble",opt=2) # get threshold based on the ensemble of models

}
}

\keyword{spatial}
\keyword{species}
\keyword{data}
