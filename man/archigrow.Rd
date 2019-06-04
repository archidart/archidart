\name{archigrow}
\alias{archigrow}
\title{Computing Growth Rates and Plotting of Vectorized Root Systems}
\description{Calculation of growth rate matrices and plotting of vectorized root systems for selected observation dates using Data Analysis of Root Tracings (DART) output files and RSA data encoded with the Root System Markup Language (RSML). In the final plots, the color code used for each link constituting a vectorized root system depends on its corresponding growth rate value. More information can be found in \cite{Delory et al (2016)}, \cite{Le Bot et al (2010)}, and \cite{Lobet et al (2015)}.  
}
\usage{
archigrow(inputlie=NULL, inputtps=NULL, inputrsml=NULL, res=NULL, unitlength="px",
          rsml.date="age", rsml.connect=TRUE, plot=TRUE, twod=NULL, export.colors=FALSE,
          unittime=NULL, unitangle="d", rotation=0, numdate=NULL, finalscale=NULL,
          coldyn=NULL, GRscale=NULL, main=NULL, xlab=NULL, ylab=NULL, zlab=NULL,
          xlim=NULL, ylim=NULL, zlim=NULL, ...)}
\arguments{
  \item{inputlie}{A character string specifying the path to the folder containing the lie files created by DART. This argument could be the same as \code{inputtps}/\code{inputrsml}.}
  \item{inputtps}{A character string specifying the path to the folder containing the tps files created by DART. This argument could be the same as \code{inputlie}/\code{inputrsml}.}
  \item{inputrsml}{A character string specifying the path to the folder containing the RSML files. This argument could be the same as \code{inputlie}/\code{inputtps}.}
  \item{res}{Mandatory for DART files only. \bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems with DART (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. For DART files, this argument must be specified if \code{unitlength="mm"} or \code{unitlength="cm"}. For RSML files, the resolution is imported from the metadata and does not need to be specified as an argument.}
  \item{unitlength}{A character string specifying the unit of length that must be used by the function to plot the root systems. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres.}
  \item{rsml.date}{Only used for RSML files. A numeric value specifying the root system age (the unit of time is defined by the user). If \code{rsml.date=NULL}, a default value of 1 is used by the function. If the root system age is stored as a continuous function along the root segments, a character string specifying the name/label of the function can be used (see details).}
  \item{rsml.connect}{Only used for RSML files. A logical value that must be specified when \code{is.null(inputrsml)=FALSE}. If \code{TRUE}, the lateral roots are connected to their parent root (see details). If \code{FALSE}, no connection is made between daughter and mother roots.}
  \item{plot}{A logical value. Should the root system(s) be plotted? \code{TRUE} is yes, \code{FALSE} is no.}
  \item{twod}{When RSML files contain 3D RSA data, should the root systems be plotted in 2D? This argument can be NULL (3D plots) or a character vector giving the name of the axes to plot (can be c(\dQuote{x}, \dQuote{y}), c(\dQuote{x}, \dQuote{z}), or c(\dQuote{y}, \dQuote{z})).}
  \item{export.colors}{A logical value. Should the color values be exported? If \code{TRUE}, a matrix giving the color code of each root at each observation date is exported. If \code{FALSE}, no color matrix is exported.}
  \item{unittime}{A character string specifying the unit of time used to express the data (all character strings are valid for this argument). For RSML files, unittime is imported from the metadata and does not need to be specified as an argument.}
  \item{unitangle}{A character string specifying the unit that must be used by the function to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{rotation}{A numeric value specifying the angle (expressed in \code{unitangle}) that must be used to perform a clockwise rotation of 2D root systems (see details).}
  \item{numdate}{A numeric value or a vector of numeric values specifying the identification number of all observation dates at which the root system(s) must be plotted (see details).}
  \item{finalscale}{A logical value that must be specified when \code{is.null(numdate)=FALSE}. If \code{TRUE}, the root system architecture at the selected observation dates is plotted with the same axis scales to that used to draw the root system architecture at the last observation date. If \code{FALSE}, the root system architecture at each observation date is plotted using different axis ranges of coordinates (\code{xlim} and \code{ylim}).}
  \item{coldyn}{A vector specifying the color(s) that must be used to plot the root system(s) at the selected observation dates. The function colors each link constituting a vectorized root system depending on its corresponding growth rate value. To do so, \code{archigrow} interpolates the color(s) contained in \code{coldyn} between the minimum and the maximum values of the computed standardized growth rate matrix (if \code{is.null(GRscale)=TRUE}) or between \code{min(GRscale)} and \code{max(GRscale)} (if \code{is.null(GRscale)=FALSE}).}
  \item{GRscale}{A numeric vector (length=2) specifying the minimum and the maximum growth rate values (expressed in \code{unitlength}/\code{unittime}) that must be used to plot each vectorized root system.}
  \item{main}{A character string giving a main title for the plot. When \code{is.null(main)=TRUE}, the default title for each plot is the name of the corresponding lie/RSML file.}
  \item{xlab}{A character string giving a label for the X axis. When \code{is.null(xlab)=TRUE}, the default label for the X axis is X (\code{unitlength}).}
  \item{ylab}{A character string giving a label for the Y axis. When \code{is.null(ylab)=TRUE}, the default label for the Y axis is Y (\code{unitlength}).}
  \item{zlab}{A character string giving a label for the Z axis. When \code{is.null(zlab)=TRUE}, the default label for the Z axis is Z (\code{unitlength}).}
  \item{xlim}{A numeric vector (length=2) giving the X limits of the plot using the same unit as \code{unitlength}.}
  \item{ylim}{A numeric vector (length=2) giving the Y limits of the plot using the same unit as \code{unitlength}.}
  \item{zlim}{A numeric vector (length=2) giving the Z limits of the plot using the same unit as \code{unitlength}.}
  \item{...}{Additional graphical parameters (see details).}
}
\details{

To run \code{archigrow} efficiently, DART and RSML files must have been saved with their appropriate extension (.lie, .tps, or .rsml).  By default (when \code{numdate=NULL}), only the root system architecture at the last observation date is plotted.

If \code{inputtps} contains more than one tps file, the number of tps files in \code{inputtps} must be equal to the number of lie files in \code{inputlie} and corresponding lie and tps files must have the same name.

After reading the RSML files located in \code{inputrsml}, the data associated with each root system are converted into a data frame possessing the same structure as the lie files created by DART.  The \code{archigrow} function then uses these data frames to plot the root system(s). 

The argument \code{rsml.date} is used to create a tps file for each root system located in \code{inputrsml}.  If the root system age is stored as a continuous function along the root segments, the \code{rsml.date} argument can be a character string specifying the name/label of the function.  For instance, the root system age at the first apparition of each node is stored as a continuous function called \dQuote{age} in RSML files containing RSA data simulated with the ArchiSimple model (\cite{Pages et al (2013)}).  The \dQuote{age} function attaches the root system age along the root segments, so that each node is associated with one age value.  An RSML library containing 3D RSA data simulated with ArchiSimple is available here: \url{https://doi.org/10.5281/zenodo.208214}.

The \code{rsml.connect} argument can be used to connect the lateral roots to their corresponding mother root.  If \code{rsml.connect=TRUE}, each point starting a lateral root is connected to the nearest point located on its mother root.  DBase is approximated by the distance between the parent root base and the point located on the mother root that is closest to the point starting a lateral root.  Using RSML files, only RSA data associated with roots possessing a branching order lower or equal to 6 are used for the computation of RSA traits.

Due to technical limitations, it is sometimes easier to take a picture of a root system if it is rotated before image acquisition.  In this case, the vectorized root system depicts a rotation compared with the natural plant root system.  In order to make a correction, one can use \code{rotation} to specify an angle value expressed in \code{unitangle} to rotate the vectorized root system clockwise before plotting.  

Additional graphical parameters can be used to personalize the graphical outputs (see \code{\link{par}}).

Growth rates are calculated as follows:  
\itemize{
\item For the first observation date, it is calculated as the ratio of the root length to the root system age.  
\item For other observation dates (t), it is calculated as the difference between the root length at time t and t-1 divided by the difference between the root system age at time t and t-1.}
}

\value{
Returns a list of secondary lists. Each element of the primary list is named as its corresponding lie/RSML file. For each lie/RSML file, the secondary list contains the following items:

\item{unit}{The unit of length/time used to express the results.}
\item{time}{A data frame containing the following columns: \bold{Date} (the identification number of each observation date) and \bold{Age} (the root system age).}
\item{GR}{A matrix containing the growth rate values computed by the function. The matrix contains the following columns: \bold{Root} (the identification number of each root constituting a vectorized root system), \bold{Ord} (the branching order; only for DART files or for RSML files with \code{rsml.connect=TRUE}), \bold{DateX} (the growth rate value of each root at the observation date X).}
\item{colors}{Only if \code{export.colors=TRUE}. A matrix containing the colors used to plot the root system. The matrix contains the following columns: \bold{Root} (the identification number of each root constituting a vectorized root system), \bold{DateX} (the color used to plot each root at the observation date X).}

If \code{plot=TRUE}, returns a plot for each vectorized root system located in \code{inputlie}/\code{inputrsml}.
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
##---------------------------------------------
## Create a function (legendGR) to add a legend
##---------------------------------------------

legendGR<-function(list, coldyn, GRscale=NULL, raster=c(0.1,0.1,0.6,1), textx=0.8, l=10, ...){
start.col<-ncol(list$GR)-length(grep("Date", colnames(list$GR)))+1
end.col<-ncol(list$GR)
growthrate<-list$GR[,start.col:end.col]
colorlegend<-colorRampPalette(coldyn)(1000)
legendimage<-as.raster(matrix(rev(colorlegend), ncol=1))
par(mar=c(1,1,2,1))
plot(c(0,1),c(0,1),type="n", axes=FALSE, xlab="", ylab="", ...)
if (is.null(GRscale)==TRUE) {text(x=textx, y=seq(raster[2],raster[4],l=l),
labels=round(seq(round(min(growthrate),1),round(max(growthrate),1),l=l),1),...)}
else {text(x=textx, y=seq(raster[2],raster[4],l=l), labels=round(seq(round(min(GRscale),1),
round(max(GRscale),1),l=l),1),...)}
rasterImage(legendimage,xleft=raster[1],ybottom=raster[2],xright=raster[3],ytop=raster[4])}

## Locate folder with DART and RSML files
path <- system.file("extdata", package="archiDART")

##----------------
## DART files only
##----------------

layout(matrix(1:2,ncol=2),widths = c(3,1),heights = c(1,1))
par(mar=c(5,4,4,2)+0.1)
res1a <- archigrow(inputlie=path, inputtps=path, res=75, unittime="day",
unitlength="cm", plot=TRUE, export.colors=TRUE, coldyn=c("blue", "orange", "red"),
las=1, bty="l", asp=1, xaxp=c(0,30,3), lwd=2)
legendGR(list=res1a$ch7, coldyn=c("blue", "orange", "red"))

##--------------------
## DART and RSML files
##--------------------

layout(1)
par(mar=c(5,4,4,2)+0.1)
res1b <- archigrow(inputlie=path, inputtps=path, inputrsml=path, res=75, unittime="day",
unitlength="cm", rsml.connect=TRUE, rsml.date="age", plot=TRUE, export.colors=TRUE,
coldyn=c("black", "orange", "red"), las=1, bty="l", asp=1, xaxp=c(0,30,3), lwd=2)
unit<-res1b$`monocot-archisimple`$unit
legendGR(list=res1b$`monocot-archisimple`, coldyn=c("black", "orange", "red"),
main=paste("Growth rate monocot-archisimple (", unit, ")", sep=""), raster=c(0.2,0,0.7,1))
}