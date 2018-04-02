\name{root}
\alias{root}
\title{Computing Attributes for Each Root}
\description{Computing attributes for each individual root using root system architecture data encoded with the Root System Markup Language (RSML) or exported by Data Analysis of Root Tracings (DART). More information can be found in \cite{Le Bot et al (2010)}, \cite{Lobet et al (2015)}, and \cite{Delory et al (2018)}.  
}
\usage{
root(x, unitangle="d", vertical3d="y", last=TRUE, show.progress=FALSE)
}
\arguments{
  \item{x}{A rsmlToTable or dartToTable object. Mandatory.}
  \item{unitangle}{A character string specifying the unit that must be used to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{vertical3d}{The vertical axis for 3D RSA data (x, y, or z).}
  \item{last}{Should the calculations be performed for the last observation date only? Yes is TRUE, No is FALSE.}
  \item{show.progress}{A logical value. Should a progress bar be shown in the R console? Yes is TRUE, No is FALSE.}
}
\details{
If a rsmlToTable object was created with \code{rsml.connect=FALSE}, DBase values are approximated by the distance between the parent root base and the point located on the mother root that is closest to the point starting a lateral root.

\code{show.progress=TRUE} can only be used if more than one root system is stored in \code{x}.
}
\value{
Returns a data frame. Each line of the data frame gives the attributes of one root. The data frame includes the following columns: 

\bold{file}: the name of each DART or RSML file.

\bold{plant}: the plant identification number. Only if \code{x} is a rsmlToTable object.

\bold{root}: the root identification number.

\bold{order}: the root branching order.

\bold{parentroot}: the identification number of the parent root.

\bold{time}: the root system age.

\bold{DBase}: the distance between the branching point to the parent root base.

\bold{length}: the root length.

\bold{mean.diameter}: the mean root diameter. Only if \code{x} is a rsmlToTable object.

\bold{sd.diameter}: the standard deviation of the root diameter. Only if \code{x} is a rsmlToTable object.

\bold{nlat}: the number of lateral roots.

\bold{branching.angle}: the branching angle expressed in \code{unitangle}.

\bold{tortuosity}: the ratio of the final root length to the Euclidean distance between the branching point and the apex of the root.

\bold{growth}: the growth rate.

\bold{surface}: the total root surface area. Only if \code{x} is a rsmlToTable object.

\bold{volume}: the total root volume. Only if \code{x} is a rsmlToTable object.

\bold{lauz}: the length of the unbranched apical zone (see \cite{Pages et al (2010)}).
}
\references{
Delory B.M., Li M., Topp C.N., Lobet G. (2018). archiDART v3.0: a new data analysis pipeline allowing the topological analysis of plant root systems, \emph{F1000Research}, 7:22, DOI: \url{http://dx.doi.org/10.12688/f1000research.13541.1}

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625.

Pages L., Serra V., Draye X., Doussan C., Pierret A. (2010) Estimating root elongation rates from morphological measurements of the root tip, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0079-x. 
}
\author{Benjamin M. Delory, Guillaume Lobet
}
\seealso{\code{\link{dartToTable}}, \code{\link{rsmlToTable}}.}

\examples{
## Locate folder with RSML file
path <- system.file("extdata", package="archiDART")

table<-rsmlToTable(inputrsml=path, unitlength="cm", rsml.date="age")

#Compute attributes for all observation dates
roots<-root(table, last=FALSE)
}