\name{pca}

\alias{pca}
\alias{pca,sdmdata-method}
\alias{pca,data.frame-method}
\alias{pca,RasterStackBrick-method}
\alias{pca,SpatRaster-method}
\alias{princomp-class}

  
\title{Principle Component Analysis}

\description{
pca performs a principal components analysis (using \code{princomp} function from \code{stats} package) on the given numeric data matrix and returns the results as an object of class princomp.
}

\usage{
\S4method{pca}{sdmdata}(x,scale,filename,...)

\S4method{pca}{data.frame}(x,scale,filename,...)

\S4method{pca}{RasterStackBrick}(x,scale,filename,...)

\S4method{pca}{SpatRaster}(x,scale,filename,...)

}

\arguments{
  \item{x}{sdmdata object, or a data.frame, or a Raster (either RasterStackBrick or SpatRaster) object}
  \item{scale}{logical; specifies whether the input data should be scaled (by subtracting the variable's mean, then dividing it by its standard deviation)}
  \item{filename}{optional character; specifies a filename that should be either a CSV file when \code{x} is \code{sdmdata} or \code{data.frame}, or a Raster file when \code{x} is a Raster object}
  \item{...}{additional arguments pass to \code{princomp} function}
}
\details{
pca analysis can be considered as a way to deal with multicollinearity problem and/or reduction of the data dimention. It returns two items in a list including data, and pca. The data contains the transoformed data into priciple components (the number of components is the same as the number of variables in the input data). You can check the pca item to see how many components (e.g., first 3) should be selected (e.g., by checking loadings). For more information on the calculation, see the \code{\link[stats]{princomp}} function.
}

\value{
a list including data (a data.frame or a RasterStack depending on the type of x), and pca results (output of the princomp function)
}


\author{Babak Naimi \email{naimi.b@gmail.com}

\url{https://r-gis.net/}

\url{https://biogeoinformatics.org/}

}


\examples{
filename <- system.file('external/predictors.tif',package='sdm')

r <- rast(filename)

p <- pca(r) # p is a .pcaObject

p

plot(p@pcaObject) # or biplot(p@pcaObject)

plot(p@data)
}

\keyword{spatial}
