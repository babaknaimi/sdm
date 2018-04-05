\name{as.data.frame}

\alias{as.data.frame}
\alias{as.data.frame,sdmdata-method}


\title{Get a data.frame with record id values (rID)}

\description{
Converts a \code{sdmdata} object to a \code{data.frame}. By additional arguments, it is possible to make a quary on the dataset (see details).
}

\usage{
\S4method{as.data.frame}{sdmdata}(x, ...)
}

\arguments{
 \item{x}{sdmdata object}
 \item{...}{Additional arguments (optional, see details)} 
}
\details{
Following additional arguments optionally can be used to get a subset of data by specifying the record IDs; or using a query  by specifying the name of species, and/or the name of data groups, and/or range of times (if time information are available):

\code{ind}: an intiger vector containing the record Ids;
\code{sp}: a character vector of the name of species
\code{grp}: a character vector of the group names (e.g., 'test'; if indipendent test is available)
\code{time}: a vector of times (an appropriate time class or a character that can be converted into a time format)
}
  
\value{
data.frame
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
