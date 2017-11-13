\name{bottleneckdist}
\alias{bottleneckdist}
\title{Computing a Pairwise Bottleneck Distance Matrix}
\description{Computing a pairwise distance matrix using a bottleneck distance method for the comparison of persistence barcodes. More information can be found in \cite{Li et al (2017)}.  
}
\usage{
bottleneckdist(x, show.progress=FALSE)
}
\arguments{
  \item{x}{A perhomology object created with the \code{perhomology} function of this package. Mandatory.}
  \item{show.progress}{A logical value. Should a progress bar be shown in the R console? Yes is TRUE, No is FALSE.}
}
\details{
The bottleneck distance between two persistence barcodes is calculated with the \code{bottleneck} function of the TDA package (\cite{Fasy et al (2017)}). 
}
\value{
Returns a pairwise distance matrix.
}
\references{
Fasy B.T., Kim J., Lecci F., Maria C., Rouvreau V. (2017). TDA: Statistical Tools for Topological Data Analysis. \url{https://CRAN.R-project.org/package=TDA}.

Li M., Duncan K., Topp C.N., Chitwood D.H. (2017) Persistent homology and the branching topologies of plants, \emph{American Journal of Botany}, DOI: 10.3732/ajb.1700046.
}
\author{Benjamin M. Delory, Guillaume Lobet
}
\seealso{\code{\link{perhomology}}.}
\examples{
## Locate folder with RSML file
path <- system.file("extdata", package="archiDART")

table<-rsmlToTable(inputrsml=path, unitlength="cm", rsml.date="age")

#Persistent homology
ph<-perhomology(table, FUN="geodesic")

#Compute pairwise distance matrix
distance<-bottleneckdist(ph)
}