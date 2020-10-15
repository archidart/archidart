\name{architect}
\alias{architect}
\title{Computing Traits Describing The Global Root System Architecture}
\description{Calculation of common root system architecture (RSA) traits from Data Analysis of Root Tracings (DART) output files and RSA data encoded with the Root System Markup Language (RSML). More information can be found in \cite{Delory et al (2016)}, \cite{Le Bot et al (2010)}, and \cite{Lobet et al (2015)}.}
\usage{
architect(inputrac=NULL, inputtps=NULL, inputrsml=NULL, res=NULL, unitlength="px",
          rsml.date=NULL, rsml.connect=FALSE, vertical3d="y",
          fitter=FALSE)}
\arguments{
  \item{inputrac}{A dartToTable object or a character string specifying the path to the folder containing the rac files created by DART. If \code{inputrac} is a character string, it can be the same as \code{inputtps}/\code{inputrsml}.}
  \item{inputtps}{A character string specifying the path to the folder containing the tps files created by DART. This argument can be the same as \code{inputrac}/\code{inputrsml}. This argument will not be used if \code{inputrac} is a dartToTable object.}
  \item{inputrsml}{A rsmlToTable object or a character string specifying the path to the folder containing the RSML files. If \code{inputrsml} is a character string, it can be the same as \code{inputrac}/\code{inputtps}.}
  \item{res}{Mandatory for DART files only. \bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems with DART (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. For DART files, this argument must be specified if \code{inputrac} is not a dartToTable object and if \code{unitlength="mm"} or \code{unitlength="cm"}. For RSML files, the resolution is imported from the metadata and does not need to be specified as an argument.}
  \item{unitlength}{A character string specifying the unit of length that must be used to express the results. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres. This argument will not be used if \code{inputrac} is a dartToTable object and/or \code{inputrsml} is a rsmlToTable object.}
  \item{rsml.date}{Only used if \code{inputrsml} is not a rsmlToTable object. A numeric value specifying the root system age (the unit of time is defined by the user). If \code{rsml.date=NULL}, a default value of 1 is used by the function. If the root system age is stored as a continuous function along the root segments, a character string specifying the name/label of the function can be used (see details).}
  \item{rsml.connect}{Only used if \code{inputrsml} is not a rsmlToTable object. A logical value that must be specified when \code{is.null(inputrsml)=FALSE}. If \code{TRUE}, the lateral roots are connected to their parent root (see details). If \code{FALSE}, no connection is made between daughter and mother roots.}
  \item{vertical3d}{The vertical axis for 3D RSA data (x, y, or z).}
  \item{fitter}{A logical value. Should the Fitter topological indices be computed? Yes is TRUE, No is FALSE. See details.}
}
\details{

To run \code{architect} efficiently, DART (.rac or .tps) and RSML (.rsml) files must have been saved with their appropriate extension.

If \code{inputtps} contains a single tps file, it is used by \code{architect} for each rac file located in \code{inputrac}.  If \code{inputtps} contains more than one tps file, the number of tps files in \code{inputtps} must be equal to the number of rac files in \code{inputrac} and corresponding rac and tps files must have the same name.

After reading the RSML files located in \code{inputrsml}, the data associated with each root system are converted into a data frame possessing the same structure as the rac files created by DART.  The \code{architect} function then uses these data frames to calculate RSA traits describing the global architecture of each root system. 

The argument \code{rsml.date} is used to create a tps file for each root system located in \code{inputrsml}.  If the root system age is stored as a continuous function along the root segments, the \code{rsml.date} argument can be a character string specifying the name/label of the function.  For instance, the root system age at the first apparition of each node is stored as a continuous function called \dQuote{age} in RSML files containing RSA data simulated with the ArchiSimple model (\cite{Pages et al (2013)}).  The \dQuote{age} function attaches the root system age along the root segments, so that each node is associated with one age value.  An RSML library containing 3D RSA data simulated with ArchiSimple is available here: \url{https://doi.org/10.5281/zenodo.208214}.

The \code{rsml.connect} argument can be used to connect the lateral roots to their corresponding mother root.  If \code{rsml.connect=TRUE}, each point starting a lateral root is connected to the nearest point located on its mother root.  DBase is approximated by the distance between the parent root base and the point located on the mother root that is closest to the point starting a lateral root.  Using RSML files, only RSA data associated with roots possessing a branching order lower or equal to 7 are used for the computation of RSA traits.

Growth rates are calculated as follows:
\itemize{
\item \emph{First-order root growth rate}: for the first observation date, it is calculated as the ratio of the total first-order root length to the root system age. For other observation dates (t), it is calculated as the difference between the total first-order root length at time t and t-1 divided by the difference between the root system age at time t and t-1.
\item \emph{Lateral root growth rates}: for the first observation date, it is calculated as the ratio of the lateral root length to the root system age. For other observation dates (t), it is calculated as the difference between the lateral root length at time t and t-1 divided by the difference between the root system age at time t and t-1.
\item \emph{Growth rate of the root system}: for the first observation date, it is calculated as the ratio of the total root system length to the root system age. For other observation dates (t), it is calculated as the difference between the total root system length at time t and t-1 divided by the difference between the root system age at time t and t-1.}

If \code{inputrac} is a dartTotable object and/or \code{inputrsml} is a rsmlToTable object, topological indices described in \cite{Fitter (1986, 1987)} can be calculated (\code{fitter=TRUE}). Here, we define each topological index as in \cite{Fitter (1986, 1987)}. We recommend the reading of these papers for more information.

The magnitude of a root system is equal to the number of external links in the system. Therefore, the magnitude is also equal to the number of root apices.

The altitude of a root system is equal to the longest path length from one exterior link (or one root apex) to the root system base.

The total exterior path length is equal to the sum of all path lengths from all exterior links (or all root apices) to the root system base.
}
\value{
Returns a data frame summarizing all RSA traits calculated for each rac/RSML file located in \code{inputrac}/\code{inputrsml}. 

The data frame contains the following columns:

\bold{FileName}: (the name of the rac/RSML files).

\bold{Time}: the root system age expressed in the same unit (=unittime) as in the tps file(s) (for rac files) or in the same unit as \code{rsml.date} (for RSML files).

\bold{TRL}: the total root system length expressed in \code{unitlength}.

\bold{GRTR}: the growth rate of the root system expressed in \code{unitlength}/unittime.

\bold{L1R}: the first-order root length expressed in \code{unitlength}.

\bold{GR1R}: the first-order root growth rate expressed in \code{unitlength}/unittime.

\bold{TN1R}: the total number of first-order roots.

\bold{TNLR}: the total number of lateral roots.

\bold{TLRL}: the total lateral root length expressed in \code{unitlength}.

\bold{NxLR}: the total number of lateral roots with a branching order equal to x.

\bold{LxLR}: the total length of lateral roots with a branching order equal to x expressed in \code{unitlength}.

\bold{MLxLR}: the mean length of lateral roots with a branching order equal to x expressed in \code{unitlength}.

\bold{GRxL}: the growth rate of lateral roots with a branching order equal to x expressed in \code{unitlength}/unittime.

\bold{D2LR}: the density of secondary roots on the first-order root(s) expressed in root(s)/\code{unitlength}.

If \code{inputrac} is a dartTotable object and/or \code{inputrsml} is a rsmlToTable object, the following parameters will be calculated:

\bold{Height}: the root system height expressed in \code{unitlength}.

\bold{Width}: the root system width expressed in \code{unitlength}.

\bold{ConvexhullXY}: the area of the convex hull on the XY plane expressed in \code{unitlength^2}.

\bold{ConvexhullXZ}: the area of the convex hull on the XZ plane expressed in \code{unitlength^2}.

\bold{ConvexhullYZ}: the area of the convex hull on the YZ plane expressed in \code{unitlength^2}.

\bold{Convexhull3D}: the volume of the convex hull expressed in \code{unitlength^3}.

\bold{Magnitude}: the magnitude of the root system. See details for more information. Computed only if \code{fitter=TRUE}.

\bold{Altitude}: the altitude of the root system. See details for more information. Computed only if \code{fitter=TRUE}.

\bold{ExtPathLength}: the total exterior path length of the root system. See details for more information. Computed only if \code{fitter=TRUE}.

\bold{MDx}: the mean diameter of roots with a branching order equal to x expressed in \code{unitlength}.

\bold{MDLR}: the mean lateral root diameter expressed in \code{unitlength}.

\bold{Sx}: the total surface area of roots with a branching order equal to x expressed in \code{unitlength^2}.

\bold{Stot}: the total root surface area expressed in \code{unitlength^2}.

\bold{Vx}: the total volume of roots with a branching order equal to x expressed in \code{unitlength^3}.

\bold{Vtot}: the root system volume expressed in \code{unitlength^3}.

}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2016) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Fitter A.H. (1986) The topology and geometry of plant root systems: influence of watering rate on root system topology in \emph{Trifolium pratense}, \emph{Annals of Botany}, 58, 91-101. 

Fitter A.H. (1987) An architectural approach to the comparative ecology of plant root systems, \emph{New Phytologist}, 106, 61-77. 

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625.

Pages L., Becel C., Boukcim H., Moreau D., Nguyen C., Voisin, A-S. (2013) Calibration and evaluation of ArchiSimple, a simple model of root system architecture, \emph{Ecological Modelling}, DOI: 10.1016/j.ecolmodel.2013.11.014.
}
\author{Benjamin M. Delory, Caroline Baudson, Yves Brostaux, Guillaume Lobet, Loic Pages, Patrick du Jardin, Pierre Delaplace
}
\seealso{\code{\link{dartToTable}}, \code{\link{rsmlToTable}}.}
\examples{
## Locate folder with DART and RSML files
path <- system.file("extdata", package="archiDART")

## Compute RSA parameters for DART and RSML files
res2a <- architect(inputrac=path, inputtps=path, inputrsml=path, res=75, unitlength="cm",
rsml.connect=TRUE, rsml.date="age")

tablersml<-rsmlToTable(inputrsml=path, unitlength="cm", rsml.date="age", rsml.connect=TRUE)
res2b <- architect(inputrsml=tablersml, fitter=TRUE)
}