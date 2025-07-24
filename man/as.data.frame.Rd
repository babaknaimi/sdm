\name{as.data.frame}

\alias{as.data.frame}
\alias{as.data.frame,sdmdata-method}


\title{Get a data.frame with record id values (rID)}

\description{
Converts a \code{sdmdata} object to a \code{data.frame}. By including additional arguments, it is possible to make a query on the dataset (see details).
}

\usage{
\S4method{as.data.frame}{sdmdata}(x, ...)
}

\arguments{
 \item{x}{sdmdata object}
 \item{...}{Additional arguments (optional, see details)} 
}
\details{
The following additional arguments can optionally be used to get a subset of data by specifying record IDs; or make a query by specifying the name of species, and/or the name of data groups, and/or a range of time period (if time is available in the training data):

\code{ind}: an integer vector with record Ids;

\code{sp}: a character vector with the name(s) of species;

\code{grp}: a character vector specifying groups of data (e.g., 'test', if independent test data is available)

\code{time}: a vector of times (an appropriate time class or a character that can be converted into a time format)
}
  
\value{
data.frame
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
file <- system.file("external/data.sdd", package="sdm")

d <- read.sdm(file)

d # a sdmdata object

df <- as.data.frame(d)

head(df)

# only the records with rID == c(1,2,3):

as.data.frame(d, ind=1:3)
}


}

\keyword{spatial}
\keyword{methods}
