\name{niche}

\alias{niche}
\alias{niche,RasterStackBrick,RasterLayer-method}
\alias{niche,RasterStackBrick,SpatialPoints-method}
\alias{niche,RasterStackBrick,sdmdata-method}
\alias{niche,sdmdata-method}
\alias{plot,.nicheRaster-method}


\title{Generate and plot Ecological Niche}

\description{
This function maps the species data (either presence/absence or probability of occurrence/habitat suitability) into a two-dimensional environmental space (i.e., based on two environmental variables) to characterise ecological niche based on the specified environmental variables.
}

\usage{
niche(x,h,n,.size,plot,out,...)

}

\arguments{
  \item{x}{A \code{Raster*} object (or sdmdata) containing environmental variables}
  \item{h}{A \code{RasterLayer}, or \code{SpatialPoints}, or \code{sdmdata} object that represents species data either in the form of habitat suitability (e.g., probability of occurrence) or presence-absence (or even presence-only) data}
  \item{n}{A character vector specifying the names of environmental variables (two names) that should be used to map the ecological niche; if \code{h} is a \code{SpatialPoints} or \code{sdmdata} object, a third name may be added specifying the species name (e.g., the name of the column in \code{SpatialPointsDataFrame} contains species data)}
  \item{.size}{optional; a numeric value (default: 1e6) specifies the size of the maximum number of records should be used to generate the ecological niche map; would be useful when the \code{Raster*} object introduced in \code{x} is big, then a random sample with the specified .size will be drawn based on which the niche is generated}
  \item{plot}{logical, specifies whether the generated niche should be plotted}
  \item{out}{logical, specifies whether the niche should be returned by the function; it will be \code{TRUE} if \code{plot} is \code{FALSE}}
  \item{...}{additional arguments includingthe argument \code{gg} (see details) and other arguments that passed to the plot function}
  }
\details{
As an additional argument, a user may specify \code{gg} which is logical, specifies whether the plot should be generated using the ggplot2 package (if the package is installed), otherwise, the raster package is used to generate the plot.
  
- \code{...}: additional arguments for the \code{plot} function (e.g., xlab, ylab, main, col, ...) can be used with the function


}


\value{
 an object of class \code{.nicheRaster} that contains some information about the environmental variable, and a RasterLayer (100x100) that represents the two-dimensional ecological niche.
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}


\examples{
\dontrun{

file <- system.file("external/species.shp", package="sdm") # get the location of the species data

species <- shapefile(file) # read the shapefile

path <- system.file("external", package="sdm") # path to the folder contains the data

lst <- list.files(path=path,pattern='asc$',full.names = T) # list the name of the raster files 


# stack is a function in the raster package, to read/create a multi-layers raster dataset
preds <- stack(lst) # making a raster object

names(preds) # 4 environmental variables are used!

d <- sdmData(formula=Occurrence~., train=species, predictors=preds)

d

# fit models:
m <- sdm(Occurrence~.,data=d,methods=c('rf','glm','brt'))

# ensemble using weighted averaging based on AUC statistic:    
p1 <- ensemble(m, newdata=preds,filename='ens1.grd',setting=list(method='weighted',stat='AUC'))
plot(p1, main='Habitat Suitability in Geographic Space')

# Mapping Ecological Niche using selected two variables
niche(x=preds, h=p1, c('precipitation','temperature'))

niche(x=preds, h=p1, c('vegetation','temperature'))

# in case if you do not have the habitat suitability map but species data:

niche(x=preds, h=species, c('vegetation','temperature','Occurrence'))


niche(x=preds, h=d, n=c('vegetation','temperature','Occurrence'))



}
}

\keyword{spatial,species,data, niche}
