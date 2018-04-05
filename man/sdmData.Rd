\name{sdmData}

\alias{sdmData}
\alias{sdmData,ANY,data.frame,ANY,missing-method}
\alias{sdmData,data.frame,formula,ANY,missing-method}
\alias{sdmData,data.frame,missing,ANY,missing-method}

\alias{sdmData,ANY,data.frame,missing,missing-method}

\alias{sdmData,ANY,SpatialPoints,ANY,missing-method}
\alias{sdmData,ANY,SpatialPoints,ANY,Raster-method}

\title{creating sdm Data object}

\description{
Creates a \code{sdmdata} objects that holds species (single or multiple) and explanatory variates. In addition, more information such as spatial coordinates, time, grouping variables, and metadata (e.g., author, date, reference, etc.) can be included.
}

\usage{
sdmData(formula,train, test,predictors,bg,filename, crs,...)
}

\arguments{
  \item{formula}{Specifies which species and explanatory variables should be taken from the input data.  Other information (e.g., spatial coordinates, grouping variables, time, etc.) can be determined as well}
  \item{train}{Training data containing species observations as a \code{data.frame} or \code{SpatialPoints} or \code{SpatialPointsDataFrames}. It may contain predictor variables as well }
  \item{test}{Independent test data with the same structure as the train data}
  \item{predictors}{ explanatory variables (predictors), defined as a raster object (\code{RasterStack} or \code{RasterBrick}). Required if train data only contain species records, or background records (pseudo-absences) should be generated}
  
  
  \item{bg}{Background data (pseudo-absence), as a data.frame. It can also be a list contains the settings to generate background data (a Raster object is required in the predictors argument)}
  \item{filename}{filename of the sdm data object to store in the disk }
  \item{crs}{optional, coordinate reference system}
  \item{...}{Additional arguments (optional) that are used to create a metadata object. See details }
  }
\details{
sdmData creates a data object, for single or multiple species. It can automatically detect the variables containing species data (if a data.frame is provided in \code{train}), but it is recommended to use formula through which all species (in the left hand side, e.g., sp1+sp2+sp3 ~ .), and the explanatory variables (in the right hand side) can be determined. If there are additional information such as spatial coordinates, time, or some variables based on which the observation can be grouped, they can be determined in the right hand side of the formula in a flexsible way (e.g., \code{~ . + coords(x+y) + g(var)}; This right hand side formula, simply determines all variables (\code{.}) + \code{x} and \code{y} as spatial coordinates + grouping observations based on the variable \code{var}; for grouping, the variable (\code{var} in this example) should be categorical, i.e., \code{factor} ).

Additional arguments can be provided to determine metadata information including:
\code{author}, \code{website}, \code{citation}, \code{help}, \code{description}, \code{date}, and \code{license}

}

\value{
 an object of class \code{sdmdata}
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
# Example 1: a data.frame containing records for a species (sp) and two predictors (b15 & NDVI):

file <- system.file("external/pa_df.csv", package="sdm")

df <- read.csv(file)

head(df) 

d <- sdmData(sp~b15+NDVI,train=df)

d

# or simply:
d <- sdmData(sp~b15+NDVI,train=df)

d

#--------
# if formula is not specified, function tries to detect species and covariates, it works well only
# if dataset contains no additional columns but species and covariates!

d <- sdmData(train=df)

d

# # only right hand side of the formula is specified (one covariate), so function detects species:
d <- sdmData(~NDVI,train=df) 

d 

#----------
###########
# Example 2: a data.frame containing presence-absence records for 1 species, 4 covariates, and 
# x, y coordinates:

file <- system.file("external/pa_df_with_xy.csv", package="sdm")

df <- read.csv(file)

head(df) 

d <- sdmData(sp~b15+NDVI+categoric1+categoric2+coords(x+y),train=df) 

d
#----
# categoric1 and categoric2 are categorical variables (factors), if not sure the data.frame has 
# them as factor, it can be specified in the formula:
d <- sdmData(sp~b15+NDVI+f(categoric1)+f(categoric2)+coords(x+y),train=df) 

d
# more simple forms of the formula:
d <- sdmData(sp~.+coords(x+y),train=df) 

d

d <- sdmData(~.+coords(x+y),train=df)  # function detects the species

d
##############
# Example 3: a data.frame containing presence-absence records for 10 species:

file <- system.file("external/multi_pa_df.csv", package="sdm")

df <- read.csv(file)

head(df) 

# in the following formula, spatial coordinates columns are specified, and the rest is asked to
# be detected by the function:
d <- sdmData(~.+coords(x+y),train=df)  

d

#--- or it can be customized wich species and which covariates are needed:
d <- sdmData(sp1+sp2+sp3~b15+NDVI+f(categoric1) + coords(x+y),train=df) 

d # 3 species, 3 covariates, and coordinates
# just be careful that if you put "." in the right hand side, while not all species columns or 
# additional columns (e.g., coordinates, time) are specified in the formula, then it takes those
# columns as covariates which is NOT right!

#########
# Example 4: Spatial data:

file <- system.file("external/pa_spatial_points.shp", package="sdm") # path to a shapefile

# use a package like rgdal, or maptools, or shapefile function in package raster to read shapefile:
p <- shapefile(file)
class(p) # a "SpatialPointsDataFrame"

plot(p)

head(p) # it contains data for 3 species

# presence-absence plot for the first species (i.e., sp1)
plot(p[p@data$sp1 == 1,],col='blue',pch=16, main='Presence-Absence for sp1')

points(p[p@data$sp1 == 0,],col='red',pch=16)


# Let's read raster dataset containing predictor variables for this study area:

file <- system.file("external/predictors.grd", package="sdm") # path to a raster object

r <- brick(file)

r # a RasterBrick object including 2 rasters (covariates)

plot(r)

# now, we can use the species points and predictor rasters in sdmData function:
d <- sdmData(sp1+sp2+sp3~b15+NDVI,train=p,predictors = r)

d

##################
# Example 5: presence-only records:


file <- system.file("external/po_spatial_points.shp", package="sdm") # path to a shapefile

# use an appropriate function to read the shapefile (e.g., readOGR in rgdal, readShapeSpatial in 
# maptools, or shapefile in raster):

po <- shapefile(file)
class(po) # a "SpatialPointsDataFrame"


head(po) # it contains data for one species (sp4) and the column has only presence records!


d <- sdmData(sp4~b15+NDVI,train=po,predictors = r)

d # as you see in the type, the data is Presence-Only

### we can add another argument (i.e., bg) to generate background (pseudo-absence) records:

#------ in bg, we are going to provide a list containing the setting to generate background
#------ the setting includes n (number of background records), method (the method used for 
#------ background generation; gRandom refers to random in geographic space), and remove (whether 
#------ points located in presence sites should be removed).

d <- sdmData(sp4~b15+NDVI,train=po,predictors = r,bg=list(n=1000,method='gRandom',remove=TRUE))

d       # as you see in the type, the data is Presence-Background

# you can alternatively, put a data.frame including background records in bg!
}
}

\keyword{spatial,species,data}
