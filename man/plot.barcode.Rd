\name{plot.barcode}
\alias{plot.barcode}
\title{Plot the Persistence Barcode}
\description{Plot the persistence barcode of the topology of a root system. More information can be found in \cite{Li et al (2017)}.}
\usage{
\method{plot}{barcode}(x, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...)
}
\arguments{
  \item{x}{A barcode object. Mandatory.}
  \item{xlim}{A numeric vector (length=2) giving the X limits of the plot.}
  \item{ylim}{A numeric vector (length=2) giving the Y limits of the plot.}
  \item{xlab}{A character string giving a label for the X axis.}
  \item{ylab}{A character string giving a label for the Y axis.}
  \item{...}{Additional graphical parameters.}
}
\value{
Returns a plot.
}
\references{
Li M., Duncan K., Topp C.N., Chitwood D.H. (2017) Persistent homology and the branching topologies of plants, \emph{American Journal of Botany}, DOI: 10.3732/ajb.1700046.
}
\author{Benjamin M. Delory, Guillaume Lobet
}
\seealso{\code{\link{perhomology}}, \code{\link{par}}.
}
\examples{
## Locate folder with RSML file
path <- system.file("extdata", package="archiDART")

table<-rsmlToTable(inputrsml=path, unitlength="cm", rsml.date="age")

#Persistent homology
ph<-perhomology(table, FUN="geodesic")

#Plot the barcode
plot(ph$'monocot-archisimple_1', las=1, xlab="Geodesic distance (cm)")
}