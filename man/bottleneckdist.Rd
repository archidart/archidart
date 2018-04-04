\name{bottleneckdist}
\alias{bottleneckdist}
\title{Computing a Pairwise Bottleneck Distance Matrix}
\description{Computing a pairwise distance matrix using a bottleneck distance method for the comparison of persistence barcodes. More information can be found in \cite{Li et al (2017)} and \cite{Delory et al (2018)}.  
}
\usage{
bottleneckdist(x, show.progress=FALSE)
}
\arguments{
  \item{x}{A perhomology object created with the \code{perhomology} function of this package. Mandatory.}
  \item{show.progress}{A logical value. Should a progress bar be shown in the R console? Yes is TRUE, No is FALSE.}
}
\details{
The bottleneck distance between two persistence barcodes is calculated with the \code{bottleneck} function of the TDA package (\cite{Fasy et al (2017)}). \bold{This package is now listed in the \dQuote{Suggests} field of the description file and must be available in your R session (see \url{https://archidart.github.io/drat})}. An error message will be returned if TDA is not available.

\code{show.progress=TRUE} can only be used if \code{x} contains more than one element.
}
\value{
Returns a pairwise distance matrix.
}
\references{
Delory B.M., Li M., Topp C.N., Lobet G. (2018). archiDART v3.0: a new data analysis pipeline allowing the topological analysis of plant root systems, \emph{F1000Research}, 7:22, DOI: \url{http://dx.doi.org/10.12688/f1000research.13541.1}

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
ph<-perhomology(table)

#Compute pairwise distance matrix
if (requireNamespace("TDA", quietly=TRUE)) {distance<-bottleneckdist(ph)}
}