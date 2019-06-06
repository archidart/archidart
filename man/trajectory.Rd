\name{trajectory}
\alias{trajectory}
\title{Computing Root Growth Directions and Trajectories}
\description{Calculation of root system architecture (RSA) parameters describing the directions and the trajectories of the root system growth using Data Analysis of Root Tracings (DART) output files and RSA data encoded with the Root System Markup Language (RSML). This function allows the calculation of the branching angle of each daughter root on its corresponding mother root, statistical parameters describing the curvature of each root constituting a vectorized root system, and the root tip angles relative to a vertical direction vector at each observation date. More information can be found in \cite{Delory et al (2016)}, \cite{Le Bot et al (2010)}, and \cite{Lobet et al (2015)}.}
\usage{
trajectory(inputrac=NULL, inputlie=NULL, inputtps=NULL, inputrsml=NULL, res=NULL,
          unitlength="px", unitangle="d", rotation=0, l.brangle, l.curv, l.tipangle,
          rsml.date=NULL, vertical3d="y", plot=NULL, twod=NULL, colangle=NULL,
          export.colors=FALSE, BRscale=NULL, main=NULL, xlim=NULL, ylim=NULL,
          zlim=NULL, xlab=NULL, ylab=NULL, zlab=NULL, ...)}
\arguments{
  \item{inputrac}{A character string specifying the path to the folder containing the rac files created by DART. This argument could be the same as \code{inputlie}/\code{inputtps}/\code{inputrsml}.}
  \item{inputlie}{A character string specifying the path to the folder containing the lie files created by DART. This argument could be the same as \code{inputrac}/\code{inputtps}/\code{inputrsml}.}
  \item{inputtps}{A character string specifying the path to the folder containing the tps files created by DART. This argument could be the same as \code{inputrac}/\code{inputlie}/\code{inputrsml}.}
  \item{inputrsml}{A character string specifying the path to the folder containing the RSML files. This argument could be the same as \code{inputrac}/\code{inputlie}/\code{inputtps}.}
  \item{res}{Mandatory for DART files only. \bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems with DART (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. For DART files, this argument must be specified if \code{unitlength="mm"} or \code{unitlength="cm"}. For RSML files, the resolution is imported from the metadata and does not need to be specified as an argument.}
  \item{unitlength}{A character string specifying the unit of length that must be used to express the results. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres.}
  \item{unitangle}{A character string specifying the unit that must be used to express the calculated angles. The value acceptable for this argument could be either \dQuote{d} for degrees (default value) or \dQuote{r} for radians.}
  \item{rotation}{A numeric value specifying the angle (expressed in \code{unitangle}) that must be used to perform a clockwise rotation of 2D root systems (see details).}
  \item{l.brangle}{A numeric value specifying the minimum root length used to calculate the branching angle of a daughter root on its corresponding mother root (see details).}
  \item{l.curv}{A numeric value specifying the distance used to place equidistantly spaced interpolated points along each root before calculating the angles between the direction vectors of the successive links constituting a vectorized root (see details).}
  \item{l.tipangle}{A numeric value specifying the distance between the root tip and an interpolated point located along a root. Each interpolated point is used by the function to calculate a direction vector following the root growth direction near the apex prior to the calculation of the root tip angle relative to a vertical direction vector at each observation date.}
  \item{rsml.date}{Only used for RSML files. A numeric value specifying the root system age (the unit of time is defined by the user). If \code{rsml.date=NULL}, a default value of 1 is used by the function. If the root system age is stored as a continuous function along the root segments, a character string specifying the name/label of the function can be used (see details).}
  \item{vertical3d}{The vertical axis for 3D RSA data (x, y, or z).}
  \item{plot}{Either \dQuote{branching} or \dQuote{direction}. If \code{plot="branching"}, each root is plotted with a color code depending on the value of its branching angle. If \code{plot="direction"}, each root is plotted with a color code depending on the value of its growth direction angle.}
  \item{twod}{When RSML files contain 3D RSA data, should the root systems be plotted in 2D? This argument can be NULL (3D plots) or a character vector giving the name of the axes to plot (can be c(\dQuote{x}, \dQuote{y}), c(\dQuote{x}, \dQuote{z}), or c(\dQuote{y}, \dQuote{z})).}
    \item{colangle}{A vector specifying the color(s) that must be used to plot the root system(s). \code{trajectory} interpolates the color(s) contained in \code{colangle} between the minimum and the maximum values of the computed standardized angle matrix (if \code{is.null(BRscale)=TRUE}) or between \code{min(BRscale)} and \code{max(BRscale)} (if \code{is.null(BRscale)=FALSE}).}
  \item{export.colors}{A logical value. Should the color values be exported?}
  \item{BRscale}{A numeric vector (length=2) specifying the minimum and the maximum angle values (expressed in \code{unitangle}) that must be used to plot each vectorized root system.}
  \item{main}{A character string giving a main title for the plot. When \code{is.null(main)=TRUE}, the default title for each plot is the name of the corresponding lie/RSML file.}
  \item{xlim}{A numeric vector (length=2) giving the X limits of the plot using the same unit as \code{unitlength}.}
  \item{ylim}{A numeric vector (length=2) giving the Y limits of the plot using the same unit as \code{unitlength}.}
  \item{zlim}{A numeric vector (length=2) giving the Z limits of the plot using the same unit as \code{unitlength}.}
  \item{xlab}{A character string giving a label for the X axis. When \code{is.null(xlab)=TRUE}, the default label for the X axis is X (\code{unitlength}).}
  \item{ylab}{A character string giving a label for the Y axis. When \code{is.null(ylab)=TRUE}, the default label for the Y axis is Y (\code{unitlength}).}
  \item{zlab}{A character string giving a label for the Z axis. When \code{is.null(zlab)=TRUE}, the default label for the Z axis is Z (\code{unitlength}).}
  \item{...}{Additional graphical parameters (see details).}
}
\details{
To run \code{trajectory} efficiently, DART (.rac, .lie or .tps) and RSML (.rsml) files must have been saved with their appropriate extension.

If \code{inputtps} contains a single tps file, it is used by \code{trajectory} for each lie/rac file located in \code{inputlie}/\code{inputrac}.  The number of lie files in \code{inputlie} must be equal to the number of rac files in \code{inputrac} and corresponding rac and lie files must have the same name.  If \code{inputps} contains more than one tps file, the number of tps files in \code{inputtps} must be equal to the number of rac/lie files in \code{inputrac}/\code{inputlie} and corresponding rac/lie and tps files must have the same name.

After reading the RSML files located in \code{inputrsml}, the data associated with each root system are converted into data frames possessing the same structure as the lie/rac files created by DART.  The \code{trajectory} function then uses these data frames to compute RSA traits. 

The argument \code{rsml.date} is used to create a tps file for each root system located in \code{inputrsml}.  If the root system age is stored as a continuous function along the root segments, the \code{rsml.date} argument can be a character string specifying the name/label of the function.  For instance, the root system age at the first apparition of each node is stored as a continuous function called \dQuote{age} in RSML files containing RSA data simulated with the ArchiSimple model (\cite{Pages et al (2013)}).  The \dQuote{age} function attaches the root system age along the root segments, so that each node is associated with one age value.  An RSML library containing 3D RSA data simulated with ArchiSimple is available here: \url{https://doi.org/10.5281/zenodo.208214}.

The \code{rsml.connect} argument can be used to connect the lateral roots to their corresponding mother root.  If \code{rsml.connect=TRUE}, each point starting a lateral root is connected to the nearest point located on its mother root.  DBase is approximated by the distance between the parent root base and the point located on the mother root that is closest to the point starting a lateral root.  Using RSML files, only RSA data associated with roots possessing a branching order lower or equal to 7 are used for the computation of RSA traits.

Due to technical limitations, it is sometimes easier to take a picture of a root system if it is rotated before image acquisition.  In this case, the vectorized root system depicts a rotation compared with the natural plant root system.  In order to make a correction, one can use \code{rotation} to specify an angle value expressed in \code{unitangle} to rotate the vectorized root system clockwise before plotting.  
    
Hereafter, we will consider that the normal vector that is orthogonal to a direction vector (a,b) has the following coordinates: (b,-a).  
The direction of lateral root growth is only computed for 2D root systems. A daughter root is considered to grow at the left side of the mother root if the scalar product between a direction vector (=vector A) going from the branching point (X0) to the following point (Xd) on the daughter root and a vector (=vector B) normal to a direction vector (=vector C) going from the branching point (X0) to the following point on the mother root (Xm) is positive.  If the scalar product between A and B is negative, the daughter root is considered to grow at the right side of the mother root.  If the scalar product between A and B is equal to zero, the calculation of the scalar product between A and B will be performed again using a vector A going from the branching point (X0) to a point located on the daughter root at Xd<-Xd+1 until the calculated scalar product is different from zero.  If A is always normal to B, a random lateral root growth orientation (either left or right) is defined.  See table 2 and figure 1 of \cite{Delory et al (2016)} for more information.

The branching angle of a daughter root on its mother root is approximated by the angle between two direction vectors going from the branching point to an interpolated point located at a distance \code{l.brangle} from the branching point on the mother or on the daughter root, respectively.  The distance \code{l.brangle} between a branching point and an interpolated point is measured along the mother and the daughter roots.  A first-order root having no mother root, the reported angle is the angle measured between a vertical direction vector and a second direction vector going from the first node of the root to an interpolated point located at a distance \code{l.brangle} from the first node on the first-order root. If a root has a final length lower than \code{l.brangle}, no branching angle can be calculated and the function returns a \code{NA} value for Branching.Angle.  See table 2 and figure 1 of \cite{Delory et al (2016)} for more information.

For each root, the growth direction angle is calculated as the angle between a vertical direction vector and a second direction vector going from the branching point to the apex of the root.

The curvature of each root is evaluated by the mean and the standard deviation of the local angles formed by the direction vectors of the successive links constituting a vectorized root.  First, the function uses \code{l.curv} to interpolate equidistantly spaced points along each root constituting a vectorized root system.  Second, the angles between the direction vectors of successive links along each root are calculated. Then, the mean and the standard deviation of the calculated angles are determined for each root.  If a root has a final length lower than 2*\code{l.curv}, the function returns a \code{NA} value for Mean.Curv and SD.Curv. If a root has a final length lower than 3*\code{l.curv}, the function returns a \code{NA} value for SD.Curv as no standard deviation can be calculated on a single angle value.  See table 2 and figure 1 of \cite{Delory et al (2016)} for more information. 
}
\value{
Returns a list including:

\item{root}{A list of data frames. Each element of the list is named as its corresponding rac/lie/RSML file and contains the following columns: \bold{Root} (the identification number of each root constituting the vectorized root system), \bold{Mother} (the parent root identification number), \bold{Ord} (the branching order), \bold{DBase} (the distance between the branching point to the parent root base expressed in \code{unitlength}), \bold{FinalRootLength} (the root length at the last observation date expressed in \code{unitlength}), \bold{Tortuosity} (the ratio of the final root length to the Euclidean distance between the branching point and the apex of the root), \bold{Orientation} (the direction of the lateral root growth; only for 2D root systems), \bold{Branching.Angle} (the branching angle expressed in \code{unitangle}), \bold{Growth.Direction} (the growth direction angle expressed in \code{unitangle}), \bold{Mean.Curv} (the mean of the local angle values calculated between the direction vectors of the successive links constructed using equidistantly spaced interpolated points along each root of a vectorized root system; the calculated mean is expressed in \code{unitangle}), \bold{SD.Curv} (the standard deviation of the local angle values calculated between the direction vectors of the successive links constructed using equidistantly spaced interpolated points along each root of a vectorized root system; the calculated standard deviation is expressed in \code{unitangle}), \bold{Colors} (if \code{export.colors=TRUE}, the colors used to plot each root).}

\item{tip}{A list of data frames. Each element of the list is named as its corresponding rac/lie/RSML file and contains the following columns: \bold{Root} (the identification number of each root constituting the vectorized root system), \bold{Ang.DateX} (the calculated root tip angle relative to a vertical direction vector at the observation date X expressed in \code{unitangle}).}
}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2016) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625.

Pages L., Becel C., Boukcim H., Moreau D., Nguyen C., Voisin, A-S. (2013) Calibration and evaluation of ArchiSimple, a simple model of root system architecture, \emph{Ecological Modelling}, DOI: 10.1016/j.ecolmodel.2013.11.014.
}
\author{Benjamin M. Delory, Guillaume Lobet, Loic Pages
}
\examples{
## Locate folder with DART and RSML files
path <- system.file("extdata", package="archiDART")

##---------------------
## DART and RSML files
##---------------------

res4 <- trajectory(inputrac=path, inputlie=path, inputtps=path, inputrsml=path, res=75, 
unitlength="cm", rsml.date="age", l.brangle=1, l.curv=1, l.tipangle=0.5, vertical3d="y")

## Distribution of the calculated branching angles
hist(res4$root$ch7$Branching.Angle, breaks=seq(from=0, to=180, by=5),
main="ch7-Branching angle distribution", xlab="Branching angle (d)", las=1,
xaxp=c(0,150,15), xlim=c(0,150))

## Evolution of the root tip angle for the primary root (black) and the two longest lateral roots
## (green and red)
date<-c(1:31)
plot(date, res4$tip$ch7[1,2:ncol(res4$tip$ch7)], type="l", lwd=2, bty="l",
las=1, ylim=c(0,90), ylab="ch7 - Root tip angle (d)", xlab="Time (Num)", col="black",
main="Root tip angle", yaxp=c(0,90,9))
lines(date, res4$tip$ch7[206,2:ncol(res4$tip$ch7)], lwd=2, col="green")
lines(date, res4$tip$ch7[221,2:ncol(res4$tip$ch7)], lwd=2, col="red")

## Branching angles of lateral roots
monocot<-res4$root$'monocot-archisimple'
data<-data.frame(Factor=c(rep("Monocot", nrow(monocot))), monocot)
boxplot(Branching.Angle[Ord==2]~Factor[Ord==2], data=data, ylab="Branching angle (d)",
ylim=c(45,70), main="Branching angle laterals - Monocot", las=1, outline=FALSE,
range=0, notch=TRUE)

## Evolution of the root tip angle for the first-order roots of the monocot root system
tip<-res4$tip$'monocot-archisimple'
root<-res4$root$'monocot-archisimple'
date<-c(1:(ncol(tip)-1))
firstorderroots<-which(root$Ord==1)
colors<-colorRampPalette(c("red", "green"))(length(firstorderroots))

plot(date, tip[1,2:ncol(tip)], type="n", lwd=2, bty="l", las=1, ylim=c(0,90),
yaxp=c(0,90,9), xaxp=c(0,16,16), ylab="Root tip angle (d)", xlab="Time (Num)", col="black",
main="Monocot - Root tip angle first-order roots")

for (i in 1:length(firstorderroots)){lines(date, tip[firstorderroots[i], 2:ncol(tip)], lwd=2,
col=colors[i])}
}