\name{read.sdm}
\alias{read.sdm}
\alias{read.sdm,character-method}
\alias{write.sdm}
\alias{write.sdm,sdmModels,character-method}
\alias{write.sdm,sdmdata,character-method}
\alias{write.sdm,.sdmCorSetting,character-method}


\title{read/write sdm* object from/to a file}

\description{
Read an sdm object from a file, or write it to a file.
}

\usage{
read.sdm(filename,...)

write.sdm(x,filename,overwrite,...)


}

\arguments{
  \item{filename}{ Filename (character)}
  \item{x}{a sdm object (e.g., sdmModels,sdmdata or sdmSetting)}
  \item{overwrite}{Logical. If \code{TRUE}, "filename" will be overwritten if it exists (default is \code{FALSE}) }
  \item{...}{additional arguments}
  }
\details{
 \code{read.sdm} function reads any files that has been writed by \code{write.sdm}.  These functions use \code{\link[base]{saveRDS}} and \code{\link[base]{readRDS}} to write and read the sdm objects. Additional arguments \code{...} pass to these functions. An \code{sdmModels} object is saved to a file with an extension of ".sdm". The file extensions for \code{sdmdata} and \code{sdmSetting} object are ".sdd", and "sds", respectively.
}


\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, 39:368-375, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi

\email{naimi@r-gis.net}

\url{http://r-gis.net}

\url{http://biogeoinformatics.org}
}


\examples{
\dontrun{

file <- system.file("external/data.sdd", package="sdm")

d <- read.sdm(file)

d
# can be used to read sdm models (sdmModels) and sdmSettings as well.

write.sdm(d,'dataset') 
# extension is created for data, model and settings as .sdd, .sds, and .sdm respectively.

list.files(pattern='dataset')
}
}

\keyword{ write, write sdm}
