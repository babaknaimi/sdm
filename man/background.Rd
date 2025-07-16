\name{background}

\alias{background}
\alias{background,SpatRaster-method}
\alias{background,Raster-method}


\title{Generates background (pseudo-absence) records}

\description{
The function uses different methods to generate background or pseudo-absence records over a study area which is assumed to be the non-NA cells in the input Raster layer(s) in \code{x}. 
}

\usage{
background(x, n, method, bias, sp, setting)
}

\arguments{
  \item{x}{a spatRaster or RasterStack or RasterBrick object with explanatory (predictor) variables that will be used to fit SDMs}
  \item{n}{numeric, number of background records to sample}
  \item{method}{a character, specifies the method of background generation; can be either \code{gRandom}, \code{eRandom}, \code{gDist}, or \code{eDist}}
  \item{bias}{optional, a Raster object (SpatRaster or RasterLayer) with a single layer that specifies bias map which can ONLY be used by the method \code{gRandom} }
  \item{sp}{species presence locations (either as a SpatVector/SpatialPoints or a data.frame/matrix object); this argument is needed if the method is either \code{gDist} or \code{eDist} }
  \item{setting}{optional, a list containing additional settings required by different methods (see details)}
  
  }
\details{
The following methods are available:

- gRandom (random selection over geographical space): this method randomly selects the non-NA pixels over the study area. Same weights are given to each pixel through the random selection of points unless the \code{bias} layer is introduced by a user which is a single raster layer that specifies a weighting scheme for background generation. A pixel with a greater value in the \code{bias} layer would have a higher chance to be selected as a background record. It has been shown by some studies that if the same bias in collecting the presence records (e.g., locations that are close to roads and residential areas have higher chance to be visited for recording species presence) is used to generate background records, it can improve the performance of SDMs.

- eRandom (random selection over environmental space): this method tries to collect a uniform (i.e., evenly distributed) distribution of records over environmental gradients by sampling in the environmental space.

- gDist (random sampling weighted by geographic distance): This method uses a random selection of locations over geographical space but gives more weights to locations with larger distance to species presence locations.

- eDist (random sampling weighted by environmental distance): This method uses a random selection of locations over geographical space but gives more weights to locations with environmental conditions that are more dissimilar to the species presence locations.
}

\value{
 a data.frame with spatial coordinates of background locations and the values of predictor variables extracted over the locations.
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
  #########

# Let's read raster dataset containing predictor variables for this study area:

file <- system.file("external/predictors.tif", package = "sdm") # path to a raster object

r <- rast(file)

r # a SpatRaster object including 2 rasters (covariates)

plot(r)
#----

file <- system.file("external/po_spatial_points.shp", package = "sdm") # path to a shapefile

po <- vect(file) # spatial points with presence-only records


head(po) # it contains data for one species (sp4) and the dataset has only presence records!




b1 <- background(x = r, n = 20, method = 'gRandom') # you may specify the bias file (a raster object)

head(b1) # background records generated using gRandom method

b2 <- background(x = r, n = 20, method = 'eRandom')

head(b2) # background records generated using eRandom method

b3 <- background(x = r, n = 20, method = 'eDist', sp = po)

head(b3) # background records generated using eDist method

b4 <- background(x = r, n = 20, method = 'gDist')

head(b4) # background records generated using gDist method


}
}

\keyword{spatial}
\keyword{data}
\keyword{species}
