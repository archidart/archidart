\name{archidraw}
\alias{archidraw}
\title{Plotting Vectorized Root Systems}
\description{Plotting vectorized root systems for selected observation dates using Data Analysis of Root Tracings (DART) output files and RSA data encoded with the Root System Markup Language (RSML). More information can be found in \cite{Delory et al (2016)}, \cite{Le Bot et al (2010)}, and \cite{Lobet et al (2015)}.}
\usage{
archidraw(inputlie=NULL, inputrsml=NULL, res=NULL, unitlength="px", rsml.connect=FALSE,
          rsml.date=NULL, unitangle="d", rotation=0, numdate=NULL, finalscale=NULL,
          coldate=par("col"), twod=NULL, main=NULL, xlab=NULL, ylab=NULL, zlab=NULL, 
          xlim=NULL, ylim=NULL, zlim=NULL,...)
}
\arguments{
  \item{inputlie}{A character string specifying the path to the folder containing the lie files created by DART. This argument could be the same as \code{inputrsml}.}
  \item{inputrsml}{A character string specifying the path to the folder containing the RSML files. This argument could be the same as \code{inputlie}.}
  \item{res}{Mandatory for DART files only. \bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems with DART (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. For DART files, this argument must be specified if \code{unitlength="mm"} or \code{unitlength="cm"}. For RSML files, the resolution is imported from the metadata and does not need to be specified as an argument.}
  \item{unitlength}{A character string specifying the unit of length that must be used to plot the root systems. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres.}
  \item{rsml.connect}{Only used for RSML files. A logical value that must be specified when \code{is.null(inputrsml)=FALSE}. If \code{TRUE}, the lateral roots are connected to their parent root (see details). If \code{FALSE}, no connection is made between daughter and mother roots.}
  \item{rsml.date}{Only used for RSML files. A numeric value specifying the root system age (the unit of time is defined by the user). If \code{rsml.date=NULL}, a default value of 1 is used by the function. If the root system age is stored as a continuous function along the root segments, a character string specifying the name/label of the function can be used (see details).}
  \item{unitangle}{A character string specifying the unit that must be used to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{rotation}{A numeric value specifying the angle (expressed in \code{unitangle}) that must be used to perform a clockwise rotation of 2D root systems (see details).}
  \item{numdate}{A numeric value or a vector of numeric values specifying the identification number of all observation dates at which the root system(s) must be plotted (see details).}
  \item{finalscale}{A logical value that must be specified when \code{is.null(numdate)=FALSE}. If \code{TRUE}, the root system architecture at the selected observation dates is plotted with the same axis scales to that used to draw the root system architecture at the last observation date. If \code{FALSE}, the root system architecture at each observation date is plotted using different axis ranges of coordinates (\code{xlim} and \code{ylim}).}
  \item{coldate}{A vector specifying the color(s) that must be used to plot the root system(s) at the selected observation dates (see details).}
  \item{twod}{When RSML files contain 3D RSA data, should the root systems be plotted in 2D? This argument can be NULL (3D plots) or a character vector giving the name of the axes to plot (can be c(\dQuote{x}, \dQuote{y}), c(\dQuote{x}, \dQuote{z}), or c(\dQuote{y}, \dQuote{z})).}
  \item{main}{A character string giving a main title for the plot. When \code{main=NULL}, the default title for each plot is the name of the corresponding lie/RSML file.}
  \item{xlab}{A character string giving a label for the X axis. When \code{xlab=NULL}, the default label for the X axis is X (\code{unitlength}).}
  \item{ylab}{A character string giving a label for the Y axis. When \code{ylab=NULL}, the default label for the Y axis is Y (\code{unitlength}).}
  \item{zlab}{A character string giving a label for the Z axis. When \code{zlab=NULL}, the default label for the Z axis is Z (\code{unitlength}).}
  \item{xlim}{A numeric vector (length=2) giving the X limits of the plot using the same unit as \code{unitlength}.}
  \item{ylim}{A numeric vector (length=2) giving the Y limits of the plot using the same unit as \code{unitlength}.}
  \item{zlim}{A numeric vector (length=2) giving the Z limits of the plot using the same unit as \code{unitlength}.}
  \item{...}{Additional graphical parameters (see details).}
}
\details{
To run \code{archidraw} efficiently, DART and RSML files must have been saved with their appropriate extension (.lie or .rsml).  By default (when \code{numdate=NULL}), only the root system architecture at the last observation date is plotted.

If the number of elements specified in \code{coldate} is lower than the number of observation dates, \code{archidraw} replicates the elements stored in \code{coldate} to plot the root system(s).  When the number of elements specified in \code{coldate} is greater than the number of observation dates, only the first corresponding elements in \code{coldate} are used by \code{archidraw} to plot the root system(s).  

After reading the RSML files located in \code{inputrsml}, the data associated with each root system are converted into a data frame possessing the same structure as the lie files created by DART.  The \code{archidraw} function then uses these data frames to plot the root system(s). 

The argument \code{rsml.date} is used to create a tps file for each root system located in \code{inputrsml}.  If the root system age is stored as a continuous function along the root segments, the \code{rsml.date} argument can be a character string specifying the name/label of the function.  For instance, the root system age at the first apparition of each node is stored as a continuous function called \dQuote{age} in RSML files containing RSA data simulated with the ArchiSimple model (\cite{Pages et al (2013)}).  The \dQuote{age} function attaches the root system age along the root segments, so that each node is associated with one age value.  An RSML library containing 3D RSA data simulated with ArchiSimple is available here: \url{https://doi.org/10.5281/zenodo.208214}.

The \code{rsml.connect} argument can be used to connect the lateral roots to their corresponding mother root.  If \code{rsml.connect=TRUE}, each point starting a lateral root is connected to the nearest point located on its mother root.  DBase is approximated by the distance between the parent root base and the point located on the mother root that is closest to the point starting a lateral root.  Using RSML files, only RSA data associated with roots possessing a branching order lower or equal to 7 are used for the computation of RSA traits.

Due to technical limitations, it is sometimes easier to take a picture of a root system if it is rotated before image acquisition.  In this case, the vectorized root system depicts a rotation compared with the natural plant root system.  In order to make a correction, one can use \code{rotation} to specify an angle value expressed in \code{unitangle} to rotate the vectorized root system clockwise before plotting.  

Additional graphical parameters can be used to personalize the graphical outputs (see \code{\link{par}}).  
}
\value{Returns a plot for each vectorized root system located in \code{inputlie}/\code{inputrsml}.
}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2016) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625.

Pages L., Becel C., Boukcim H., Moreau D., Nguyen C., Voisin, A-S. (2013) Calibration and evaluation of ArchiSimple, a simple model of root system architecture, \emph{Ecological Modelling}, DOI: 10.1016/j.ecolmodel.2013.11.014.
}
\author{Benjamin M. Delory, Caroline Baudson, Yves Brostaux, Guillaume Lobet, Loic Pages, Patrick du Jardin, Pierre Delaplace
}
\seealso{\code{\link{par}}.
}
\examples{
## Locate folder with DART and RSML files
path <- system.file("extdata", package="archiDART")

##-----------------
## DART files only
##-----------------

## Plot RSA at several observation dates
## One color for each observation date
archidraw(inputlie=path, res=75, unitlength="cm", numdate=c(15,31), finalscale=TRUE,
coldate=rainbow(31), lwd=2,las=1, bty="l", asp=1, xaxp=c(0,30,3), yaxp=c(0,90,9))

##-----------------
## RSML files only
##-----------------

## Plot RSA at the last observation date
## One color for each observation date
##2D plot
archidraw(inputrsml=path, unitlength="cm", rsml.date="age", coldate=rainbow(15), lwd=2,
twod=c("x", "y"), asp=1, las=1, bty="l")
          
##---------------------
## DART and RSML files
##---------------------

archidraw(inputlie=path, inputrsml=path, res=75, unitlength="cm", rsml.date="age",
coldate="black", lwd=2, las=1, bty="l", asp=1)
}