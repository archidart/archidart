\name{dartToTable}
\alias{dartToTable}
\title{Import DART Files Into a Single Data Frame}
\description{Import DART files located in a folder into a single data frame. More information about DART files can be found in \cite{Le Bot et al (2010)}.  
}
\usage{
dartToTable(inputrac, inputlie, inputtps, res=NULL, unitlength="px",
            unitangle="d", fitter=FALSE)
}
\arguments{
  \item{inputrac}{A character string specifying the path to the folder containing the rac files created by DART. This argument could be the same as \code{inputlie}/\code{inputtps}.}
  \item{inputlie}{A character string specifying the path to the folder containing the lie files created by DART. This argument could be the same as \code{inputrac}/\code{inputtps}.}
  \item{inputtps}{A character string specifying the path to the folder containing the tps files created by DART. This argument could be the same as \code{inputrac}/\code{inputlie}.}
  \item{res}{\bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems with DART (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. This argument must be specified if \code{unitlength="mm"} or \code{unitlength="cm"}.}
  \item{unitlength}{A character string specifying the unit of length that must be used in the final data frame. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres. The function will use the metadata inside each RSML file (unit and resolution) to make the conversion.}
  \item{unitangle}{A character string specifying the unit that must be used to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{fitter}{A logical value. Should the Fitter topological indices be computed? Yes is TRUE, No is FALSE. See details.}
}
\details{
To run \code{dartToTable} efficiently, DART files must have been saved with their appropriate extension (.rac, .lie, and .tps).

If \code{fitter=TRUE}, topological indices described in \cite{Fitter (1986, 1987)} will be calculated. Here, we define each topological index as in \cite{Fitter (1986, 1987)}. We recommend the reading of these papers for more information.

The magnitude of a link is equal to the number of external links that feed into that link.

The path length of a link is equal to the number of links between a given link and the root system base.
}
\value{
Returns a dartToTable object. Each line of the data frame gives the attributes of one segment. The data frame includes the following columns:

\bold{file}: the name of each DART file located in \code{inputrac}.

\bold{root}: the root identification number.

\bold{order}: the root branching order.

\bold{parentroot}: the identification number of the parent root.

\bold{time}: the date of apparition.

\bold{bran}: the branching state (\sQuote{true} denotes a branching link while \sQuote{false} denotes a continuing link).

\bold{apic}: the apical state (\sQuote{true} denotes a segment ending a root).

\bold{x1}: the X coordinate of the node starting the segment (expressed in \code{unitlength}).

\bold{y1}: the Y coordinate of the node starting the segment (expressed in \code{unitlength}).

\bold{z1}: the Z coordinate of the node starting the segment (expressed in \code{unitlength}).

\bold{x2}: the X coordinate of the node ending the segment (expressed in \code{unitlength}).

\bold{y2}: the Y coordinate of the node ending the segment (expressed in \code{unitlength}).

\bold{z2}: the Z coordinate of the node ending the segment (expressed in \code{unitlength}).

\bold{length}: the length of the segment (expressed in \code{unitlength}).

\bold{blength}: the cumulative length measured along each root (expressed in \code{unitlength}).

\bold{orientation}: the angle between each segment and a vertical direction vector (expressed in \code{unitangle}).

\bold{growth}: the growth rate value of each segment (expressed in \code{unitlength}).

\bold{geodesic}: the geodesic distance measured along the roots (expressed in \code{unitlength}).

\bold{magnitude}: the magnitude of the segment (see \cite{Fitter et al (1987)}). Computed only if \code{fitter=TRUE}.

\bold{pathlength}: the external path length of the segment (see \cite{Fitter et al (1987)}). Computed only if \code{fitter=TRUE}.
}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2016) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Fitter A.H. (1986) The topology and geometry of plant root systems: influence of watering rate on root system topology in \emph{Trifolium pratense}, \emph{Annals of Botany}, 58, 91-101. 

Fitter A.H. (1987) An architectural approach to the comparative ecology of plant root systems, \emph{New Phytologist}, 106, 61-77. 

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.
}
\author{Benjamin M. Delory, Guillaume Lobet
}
\examples{
## Locate folder with DART files
path <- system.file("extdata", package="archiDART")

table<-dartToTable(inputrac=path, inputlie=path, inputtps=path, res=75, unitlength="cm")}