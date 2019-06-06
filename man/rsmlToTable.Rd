\name{rsmlToTable}
\alias{rsmlToTable}
\title{Import RSML Files Into a Single Data Frame}
\description{Import RSML files located in a folder into a single data frame. More information about RSML files can be found in \cite{Lobet et al (2015)} and \cite{Delory et al (2018)}.   
}
\usage{
rsmlToTable(inputrsml, unitlength="px", rsml.date=NULL, rsml.connect=TRUE,
            vertical3d="y", unitangle="d", fitter=FALSE, show.progress=FALSE)
}
\arguments{
  \item{inputrsml}{A character string specifying the path to the folder containing the RSML files. Mandatory.}
  \item{unitlength}{A character string specifying the unit of length that must be used in the final data frame. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres. The function will use the metadata inside each RSML file (unit and resolution) to make the conversion.}
  \item{rsml.date}{A numeric value specifying the root system age (the unit of time is defined by the user). If \code{rsml.date=NULL} (default), a default value of 1 is used by the function. If the root system age is stored as a continuous function along the root segments, a character string specifying the name/label of the function can be used (see details).}
  \item{rsml.connect}{A logical value. Mandatory. If \code{TRUE} (default), the lateral roots are connected to their parent root (see details). If \code{FALSE}, no connection is made between daughter and mother roots.}
  \item{vertical3d}{The vertical axis for 3D RSA data (x, y, or z).}
  \item{unitangle}{A character string specifying the unit that must be used to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{fitter}{A logical value. Should the Fitter topological indices be computed? Yes is TRUE, No is FALSE. See details.}
  \item{show.progress}{A logical value. Should a progress bar be shown in the R console? Yes is TRUE, No is FALSE.}
}
\details{
To run \code{rsmlToTable} efficiently, RSML (.rsml) files must have been saved with their appropriate extension.

The argument \code{rsml.date} is used to create a tps file for each root system located in \code{inputrsml}.  If the root system age is stored as a continuous function along the root segments, the \code{rsml.date} argument can be a character string specifying the name/label of the function.  For instance, the root system age at the first apparition of each node is stored as a continuous function called \dQuote{age} in RSML files containing RSA data simulated with the ArchiSimple model (\cite{Pages et al (2013)}).  The \dQuote{age} function attaches the root system age along the root segments, so that each node is associated with one age value.  An RSML library containing 3D RSA data simulated with ArchiSimple is available here: \url{https://doi.org/10.5281/zenodo.208214}.

The \code{rsml.connect} argument can be used to connect the lateral roots to their corresponding mother root.  If \code{rsml.connect=TRUE}, each point starting a lateral root is connected to the nearest point located on its mother root.  DBase is approximated by the distance between the parent root base and the point located on the mother root that is closest to the point starting a lateral root.  Using RSML files, only RSA data associated with roots possessing a branching order lower or equal to 7 are used for the computation of RSA traits.

\code{show.progress=TRUE} can only be used if more than one RSML file is stored in \code{inputrsml}.

If \code{fitter=TRUE}, topological indices described in \cite{Fitter (1986, 1987)} will be calculated. Here, we define each topological index as in \cite{Fitter (1986, 1987)}. We recommend the reading of these papers for more information.

The magnitude of a link is equal to the number of external links that feed into that link.

The path length of a link is equal to the number of links between a given link and the root system base.
}
\value{
Returns a rsmlToTable object. Each line of the data frame gives the attributes of one segment. The data frame includes the following columns:

\bold{file}: the name of each RSML file located in \code{inputrsml}.

\bold{plant}: the plant identification number.

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

\bold{diameter1}: the root diameter at the node starting the segment (expressed in \code{unitlength}).

\bold{diameter2}: the root diameter at the node ending the segment (expressed in \code{unitlength}).

\bold{length}: the length of the segment (expressed in \code{unitlength}).

\bold{blength}: the cumulative length measured along each root (expressed in \code{unitlength}).

\bold{surface}: the surface area of the root segment (expressed in \code{unitlength^2}).

\bold{volume}: the volume of the root segment (expressed in \code{unitlength^3}).

\bold{orientation}: the angle between each segment and a vertical direction vector (expressed in \code{unitangle}).

\bold{growth}: the growth rate value of each segment (expressed in \code{unitlength}).

\bold{geodesic}: the geodesic distance measured along the roots (expressed in \code{unitlength}).

\bold{magnitude}: the magnitude of the segment. See details for more information. Computed only if \code{fitter=TRUE}.

\bold{pathlength}: the path length of the segment. See details for more information. Computed only if \code{fitter=TRUE}.
}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2016) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Delory B.M., Li M., Topp C.N., Lobet G. (2018). archiDART v3.0: a new data analysis pipeline allowing the topological analysis of plant root systems, \emph{F1000Research}, 7:22, DOI: \url{http://dx.doi.org/10.12688/f1000research.13541.1}

Fitter A.H. (1986) The topology and geometry of plant root systems: influence of watering rate on root system topology in \emph{Trifolium pratense}, \emph{Annals of Botany}, 58, 91-101. 

Fitter A.H. (1987) An architectural approach to the comparative ecology of plant root systems, \emph{New Phytologist}, 106, 61-77. 

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625.

Pages L., Becel C., Boukcim H., Moreau D., Nguyen C., Voisin, A-S. (2013) Calibration and evaluation of ArchiSimple, a simple model of root system architecture, \emph{Ecological Modelling}, DOI: 10.1016/j.ecolmodel.2013.11.014.
}
\author{Benjamin M. Delory, Guillaume Lobet
}
\examples{
## Locate folder with RSML file
path <- system.file("extdata", package="archiDART")

table<-rsmlToTable(inputrsml=path, unitlength="cm", rsml.date="age")}