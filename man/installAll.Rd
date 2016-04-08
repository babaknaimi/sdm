\name{installAll}

\alias{installAll}
\alias{installAll,ANY-method}


\title{install all packages that may be required by the package}

\description{
This function facilitates to install the required packages that some functions are dependent on in the sdm package. It first checks whether the package is already installed, and if not, it installs the packages. If update=TRUE is used, the packages re-installed if they were already installed.
}

\usage{
installAll(pkgs,update,...)

}

\arguments{
  \item{pkgs}{optional. the user provided list of packages (not required for the purpose of this function)}
  \item{update}{logical (default=FALSE), specifies whether the packages re-installed if they are already installed on the machine}
  \item{...}{Additional arguments passed to the \code{\link{install.packages}} function }
  }



\references{
Naimi, B., Araujo, M.B. (2016) sdm: a reproducible and extensible R platform for species distribution modelling, Ecography, DOI: 10.1111/ecog.01881
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://biogeoinformatics.org}
}

\seealso{# }

\examples{
\dontrun{

installAll()

}
}

\keyword{utilities, data}