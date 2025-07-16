\name{aoa}

\alias{aoa}
\alias{aoa,SpatRaster,sdmdata-method}
\alias{aoa,SpatRaster,sdmModels-method}

\title{Area of Applicability}

\description{
aoa (Area of Applicability) measures whether the values of pixels in a SpatRaster object (x) with environmental layers used to train SDMs are within the range of variables contributed in the modelling.
}

\usage{
aoa(x, d, vi = NULL)

}

\arguments{
  \item{x}{a \code{SpatRaster} object with predictors used to train SDMs}
  \item{d}{a  \code{sdmdata} or \code{sdmModels} object}
  \item{vi}{optional; a numeric vector of variable importance values of the variables in \code{x}; vi should have the same length as the number of variables in \code{x}}
  }
\details{
The output of this function is a raster with values ranging between 0 and 1. It can be inferred as the degree of similarity between values of the environmental variables at each pixel and training range. When the values of all variables in \code{x} are within the range used to train the models (training range), the output of \code{aoa} is 1. Less than 1 refers to degree of dissimilarity to the range (it may be the proportion of variables outside of the training range; smaller value is more dissimilar) .

The \code{aoa} function can be used when the predict or ensemble function is used to predict or project distribution in a new area or a new time (e.g., future) where it is likely to have pixels with values outside of the training range.
By using the variable importance, higher weights are given to more important variables to assess similarity.

}

\value{
 a \code{SpatRaster} object 
}

\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://www.r-gis.net/}

\url{https://www.biogeoinformatics.org}
}



\examples{
\dontrun{


file <- system.file("external/species.shp", package = "sdm") # get the location of the species data

species <- vect(x = file) # read the shapefile

path <- system.file("external", package = "sdm") # path to the folder containing the data

lst <- list.files(path = path, pattern = 'asc$', full.names = T) # list the name of the raster files 


# rast is a function in the terra package, to read/create a multi-layers SpatRaster dataset
preds <- rast(x = lst) # making a SpatRaster object

modeling_extent <- c(350000, 600000, 4000000, 4200000) # extent to subset the species and preds for model building

species_cropped <- crop(x = species, y = modeling_extent) # crop species to the modeling extent
preds_cropped <- crop(x = preds, y = modeling_extent) # crop preds to the modeling extent

d <- sdmData(formula = Occurrence ~., train = species_cropped, predictors = preds_cropped)

aoa_layer <- aoa(x = preds, d = d) # provide whole extent predictor layers (preds) and the sdmdata object (d)

aoa_layer # a SpatRaster object with values ranging from 0 (highly dissimilar) to 1 (very similar)

plot(aoa_layer)

}


}

\keyword{spatial}
\keyword{data}