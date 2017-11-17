\name{perhomology}
\alias{perhomology}
\title{Topological Analysis using Persistent Homology}
\description{Computing the birth and death of homology group components (path-connected components) using root system architecture data encoded with the Root System Markup Language (RSML) or exported by Data Analysis of Root Tracings (DART). More information can be found in \cite{Le Bot et al (2010)}, \cite{Lobet et al (2015)} and \cite{Li et al (2017)}.  
}
\usage{
perhomology(x, show.progress=FALSE)
}
\arguments{
  \item{x}{A rsmlToTable or dartToTable object. Mandatory.}
  \item{show.progress}{A logical value. Should a progress bar be shown in the R console? Yes is TRUE, No is FALSE.}
}
\details{
In mathematics, homology refers to homology groups. Because roots are path-connected components, only zero order homology (H0) branches are considered here.

The persistence of each H0 branch is computed using a geodesic distance function. More information can be found in \cite{Li et al (2017)}. 
}
\value{
Returns a list of matrices. The list is a perhomology object. Each matrix of the list is a barcode object. Each element of the list is named as its corresponding RSML/DART file. Each matrix is associated with one root system only and possesses the following columns: 

\bold{dimension}: the homology group (0 for path-connected components).

\bold{birth}: the birth of an H0 branch along the distance function.

\bold{death}: the death of an H0 branch along the distance function.
}
\references{
Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.

Li M., Duncan K., Topp C.N., Chitwood D.H. (2017) Persistent homology and the branching topologies of plants, \emph{American Journal of Botany}, DOI: 10.3732/ajb.1700046.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625.
}
\author{Benjamin M. Delory, Guillaume Lobet
}
\seealso{\code{\link{dartToTable}}, \code{\link{rsmlToTable}}, \code{\link{plot.barcode}}.}
\examples{
## Locate folder with RSML file
path <- system.file("extdata", package="archiDART")

table<-rsmlToTable(inputrsml=path, unitlength="cm", rsml.date="age")

#Persistent homology
ph1<-perhomology(table)
plot(ph1$'monocot-archisimple_1', las=1, xlab="Geodesic distance (cm)",
ylab=expression(H[0]))
}