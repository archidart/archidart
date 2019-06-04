\name{latdist}
\alias{latdist}
\title{Computing Lateral Root Length and Density Distribution}
\description{Analysing the lateral root length and density distribution on each mother root using Data Analysis of Root Tracings (DART) output files and RSA data encoded with the Root System Markup Language (RSML). More information can be found in \cite{Delory et al (2016)}, \cite{Le Bot et al (2010)}, and \cite{Lobet et al (2015)}.  
}
\usage{
latdist(inputrac=NULL, inputrsml=NULL, output=c("lrd","dtp"), res=NULL, unitlength="px",
        int.length=NULL, interpol=NULL, rsml.connect=TRUE)
}
\arguments{
  \item{inputrac}{A character string specifying the path to the folder containing the rac files created by DART. This argument could be the same as \code{inputrsml}.}
  \item{inputrsml}{A character string specifying the path to the folder containing the RSML files. This argument could be the same as \code{inputrac}.}
  \item{output}{A character string specifying the type of output that is returned by the function. Two values are acceptable for this argument: \dQuote{lrd} (default value; root length and density distribution) and \dQuote{dtp} (distance between neighbouring lateral roots). See details.}
  \item{res}{Mandatory for DART files only. \bold{If images were acquired with a flatbed scanner}: a numeric value specifying the resolution of the images used to vectorize the root systems with DART (resolution is expressed in dots/inch). \bold{If images were acquired with a digital camera}: the ratio between the length of a reference object located on the image expressed in pixels and the actual length of the same object expressed in inches. For DART files, this argument must be specified if \code{unitlength="mm"} or \code{unitlength="cm"}. For RSML files, the resolution is imported from the metadata and does not need to be specified as an argument.}
  \item{unitlength}{A character string specifying the unit of length that must be used to plot the root systems. The value acceptable for this argument could be either \dQuote{px} for pixels, \dQuote{cm} for centimetres or \dQuote{mm} for millimetres.}
  \item{int.length}{Mandatory when \code{output="lrd"}. A numeric value specifying the interval length (expressed in \code{unitlength}) used to calculate local lateral root length and density (see details).}
  \item{interpol}{A numeric value specifying the number of points used on each mother root to calculate local lateral root length and density. By default, \code{interpol=NULL} means that the function calculates local lateral root length and density at each branching point of each mother root. A numeric value can be specified when \code{output="lrd"} (see details).}
  \item{rsml.connect}{Only used for RSML files. A logical value that must be specified when \code{is.null(inputrsml)=FALSE}. If \code{TRUE}, the lateral roots are connected to their parent root (see details). If \code{FALSE}, no connection is made between daughter and mother roots.}
}
\details{
\bold{DBase}: distance between a branching point and the parent root base.

To run \code{latdist} efficiently, DART (.rac) and RSML (.rsml) files must have been saved with their appropriate extension.

After reading the RSML files located in \code{inputrsml}, the data associated with each root system are converted into a data frame possessing the same structure as the rac files created by DART. The \code{latdist} function then uses these data frames to compute lateral root length and density distribution.

The \code{rsml.connect} argument can be used to connect the lateral roots to their corresponding mother root.  If \code{rsml.connect=TRUE}, each point starting a lateral root is connected to the nearest point located on its mother root.  DBase is approximated by the distance between the parent root base and the point located on the mother root that is closest to the point starting a lateral root.  Using RSML files, only RSA data associated with roots possessing a branching order lower or equal to 6 are used for the computation of RSA traits.

If \code{output="lrd"}, the function uses \code{interpol} to select the appropriate DBase values for which the calculation of RSA parameters should be performed.  Then, a lateral root density and a total lateral root length are calculated on intervals possessing a length equal to \code{int.length} and centred on each DBase value.
}
\value{
Returns a list including:

\item{root}{A list of data frames. Each element of the list is named as its corresponding rac/RSML file and contains the following columns: \bold{Root} (the identification number of each root constituting a vectorized root system), \bold{Ord} (the branching order), \bold{LatRootNum} (the number of lateral roots), \bold{FinalRootLength} (the root length at the last observation date expressed in \code{unitlength}), \bold{LatRootDensity} (the lateral root density calculated as LatRootNum/FinalRootLength and expressed in root(s)/\code{unitlength}).}

\item{results}{A list of secondary lists. Each element of the primary list (results) is named as its corresponding rac/RSML file. The secondary lists contain as much elements as roots constituting a vectorized root system. Each element of the secondary lists could be either a \code{NULL} value if a root does not have lateral roots or if no DBase value is comprised in [\code{int.length}/2, FinalRootLength-(\code{int.length}/2)], or a data frame made of two or three columns depending on the \code{output} value. 

\bold{If \code{output="lrd"}}: \bold{DBase} (the distance between the branching point and the parent root base expressed in \code{unitlength}), \bold{LRD} (Local lateral root density expressed in root(s)/\code{unitlength}), \bold{LRL} (Local lateral root length expressed in \code{unitlength}/\code{unitlength}).

\bold{If \code{output="dtp"}}: \bold{DBase} (the distance between the branching point and the parent root base expressed in \code{unitlength}), \bold{DTP} (the distance to the previous lateral root expressed in \code{unitlength}).}
}
\references{
Delory B.M., Baudson C., Brostaux Y., Lobet G., du Jardin P., Pages L., Delaplace P. (2016) archiDART: an R package for the automated computation of plant root architectural traits, \emph{Plant and Soil}, DOI: 10.1007/s11104-015-2673-4.

Le Bot J., Serra V., Fabre J., Draye X., Adamowicz S., Pages L. (2010) DART: a software to analyse root system architecture and development from captured images, \emph{Plant and Soil}, DOI: 10.1007/s11104-009-0005-2.

Lobet G., Pound M.P., Diener J., Pradal C., Draye X., Godin C., Javaux M., Leitner D., Meunier F., Nacry P., Pridmore T.P., Schnepf A. (2015) Root System Markup Language: Toward a Unified Root Architecture Description Language, \emph{Plant Physiology}, DOI: 10.1104/pp.114.253625.
}
\author{Benjamin M. Delory, Guillaume Lobet, Loic Pages
}
\examples{
## Locate folder with DART and RSML files
path <- system.file("extdata", package="archiDART")

##-----------------
## DART files only
##-----------------

## Plotting the results for the primary root
## output="lrd", interpol=NULL, int.length=5
res3a <- latdist(inputrac=path, output="lrd", res=75, unitlength="cm",
int.length=5)

plot(res3a$results$ch7[[1]]$DBase, res3a$results$ch7[[1]]$LRD, pch=16, 
xlab="DBase (cm)",ylab="Lateral root density (roots/cm)", main="LRD-interpol=NULL-int.length=5",
las=1, bty="l", xaxp=c(0,90,9))

plot(res3a$results$ch7[[1]]$DBase, res3a$results$ch7[[1]]$LRL, pch=16, 
xlab="DBase (cm)",ylab="Lateral root length (cm/cm)", main="LRL-interpol=NULL-int.length=5",
las=1, bty="l", xaxp=c(0,90,9))

## output="lrd", interpol=1000, int.length=5
res3b <- latdist(inputrac=path, output="lrd", res=75, unitlength="cm",
int.length=5, interpol=1000)

plot(res3b$results$ch7[[1]]$DBase, res3b$results$ch7[[1]]$LRD, pch=16, 
xlab="DBase (cm)", ylab="Lateral root density (roots/cm)", main="LRD-interpol=1000-int.length=5", 
las=1, bty="l", xaxp=c(0,90,9))

plot(res3b$results$ch7[[1]]$DBase, res3b$results$ch7[[1]]$LRL, pch=16, 
xlab="DBase (cm)", ylab="Lateral root length (cm/cm)", main="LRL-interpol=1000-int.length=5", 
las=1, bty="l", xaxp=c(0,90,9))

## output="dtp"
res3c <- latdist(inputrac=path, output="dtp", res=75, unitlength="cm")

plot(res3c$results$ch7[[1]]$DBase, res3c$results$ch7[[1]]$DTP, pch=16, 
xlab="DBase (cm)", ylab="Distance to the previous root (cm)", main="DTP", las=1, bty="l",
xaxp=c(0,90,9))

##-----------------
## RSML files only
##-----------------

## output="lrd", interpol=200, int.length=1
res3d <- latdist(inputrsml=path, unitlength="cm", output="lrd", int.length=1, interpol=200,
rsml.connect=TRUE)

##output="dtp"
res3e <- latdist(inputrsml=path, output="dtp", unitlength="cm")

##---------------------
## DART and RSML files
##---------------------

## output="lrd", interpol=200, int.length=2
res3f <- latdist(inputrac=path, inputrsml=path, output="lrd", res=75, unitlength="cm",
int.length=2, interpol=200, rsml.connect=TRUE)

## output="dtp"
res3g <- latdist(inputrac=path, inputrsml=path, output="dtp", res=75, unitlength="cm")
}