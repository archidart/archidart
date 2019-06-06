archidraw<-function(inputlie=NULL, inputrsml=NULL, res=NULL, unitlength="px", rsml.connect=FALSE, rsml.date=NULL, unitangle="d", rotation=0, numdate=NULL, finalscale=NULL, coldate=par("col"), twod=NULL, main=NULL, xlab=NULL, ylab=NULL, zlab=NULL, xlim=NULL, ylim=NULL, zlim=NULL,...){
    
    # Errors interception
    
    if (is.null(inputlie)==TRUE & is.null(inputrsml)==TRUE){stop("inputlie and/or inputrsml must be provided")}
  
    if (is.null(inputlie)==FALSE) {if (mode(inputlie)!="character"){stop("mode(inputlie) must be character")}}
    
    if (is.null(inputrsml)==FALSE) {if (mode(inputrsml)!="character"){stop("mode(inputrsml) must be character")}}
    
    if (is.null(inputlie)==FALSE & is.null(res)==TRUE & unitlength!="px"){stop("If unitlength is not px, res must be specified")}
    if (is.null(res)==FALSE){
      if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
      if (res<=0){stop("res must be a positive value")}}
    
    if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
    if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
    
    if (mode(unitangle)!="character"){stop("mode(unitangle) must be character")}
    if(unitangle=="d"|unitangle=="r") {} else {stop("unitangle must be either d (degrees) or r (radians)")}
    
    if (mode(rotation)!="numeric"){stop("mode(rotation) must be numeric")}
    if (rotation<0){stop("rotation must be a positive value")}
    
    if (mode(rsml.connect)!="logical"){stop("mode(rsml.connect) must be logical")}
    
    if (is.null(numdate)==FALSE)
    {if (mode(numdate)!="numeric"){stop("mode(numdate) must be numeric")}
     for (i in 1:length(numdate)){if (numdate[i]<=0){stop("numdate must be either a positive value or a vector of positive values")}}
     numdate.sort<-sort(numdate)
     for (i in 1:length(numdate)) {if (numdate[i]!=numdate.sort[i]){stop("Numeric elements in numdate must be sorted by increasing values")}}}
    
    if (is.null(finalscale)==FALSE) {if (mode(finalscale)!="logical"){stop("mode(finalscale) must be logical")}}
  
    if (is.null(numdate)==FALSE & is.null(finalscale)==TRUE) {stop("If numdate is not NULL, finascale must be specified")}
  
    if (is.null(rsml.date)==FALSE){
      if (is.character(rsml.date)==TRUE|is.numeric(rsml.date)==TRUE){} else {stop("If rsml.date is not NULL, rsml.date must be a character string or a positive numeric value")}
      if (is.numeric(rsml.date)==TRUE){if (rsml.date<=0|length(rsml.date)>1){stop("If mode(rsml.date) is numeric, rsml.date must be a single positive value")}}}
  
    if (is.null(twod)==FALSE){
      if (mode(twod)!="character"){stop("mode(twod) must be character")}
      if (length(twod)!=2){stop("twod must be a vector of 2 character elements")}
      twod<-sort(twod)
      if (all.equal(twod, c("x", "y"))==TRUE|all.equal(twod, c("x", "z"))==TRUE|all.equal(twod, c("y", "z"))==TRUE) {} else {stop("twod must be c(x,y), c(x,z), or c(y,z)")}}
  
  # Reading of DART and rsml files
    
    if (is.null(inputlie)==FALSE){
      filenames.lie<-mixedsort(list.files(path=inputlie, pattern="\\.lie$"))
      path.lie<-rep(inputlie, length.out=length(filenames.lie))
      filenameslie<-sub(x=filenames.lie, pattern="\\.lie$", replacement="")
      message(paste("Number of DART lie files in inputlie:", length(filenames.lie), sep=" "))}
    
    if (is.null(inputrsml)==FALSE) {
      filenames.rsml<-mixedsort(list.files(path=inputrsml, pattern="\\.rsml$"))
      path.rsml<-rep(inputrsml, length.out=length(filenames.rsml))
      filenamesrsml<-sub(x=filenames.rsml, pattern="\\.rsml$", replacement="")
      message(paste("Number of rsml files in inputrsml:", length(filenames.rsml), sep=" "))}
    
    if (is.null(inputrsml)==TRUE){
      if (length(filenames.lie)==0){stop("There is no lie file in inputlie")}}
    else {
      if (is.null(inputlie)==TRUE){if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}
      else{
        if (length(filenames.lie)==0){stop("There is no lie file in inputlie")}
        if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}}
    
    if (is.null(inputrsml)==TRUE){ # Only DART files

    LIE<-lapply(paste(path.lie, "/", filenames.lie, sep=""), read.table, header=TRUE)
    
    #Unit conversion DART files
    if (unitlength=="mm") {cunit<-(10*cm(1)/res)}
    if (unitlength=="cm") {cunit<-(cm(1)/res)}
    if (unitlength=="px") {cunit<-1}}
    
    else {
      
      if (is.null(inputlie)==TRUE){ # Only rsml files
        
        LIE<-list()
        res1<-c()
        unitlength1<-c()
        filenameslie<-c()
        RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect)
        for (i in 1:length(RSML)){
          res1<-append(res1, rep(as.numeric(RSML[[i]]$resolution), length(RSML[[i]]$lie)))
          unitlength1<-append(unitlength1, rep(as.character(RSML[[i]]$length), length(RSML[[i]]$lie)))
          LIE<-append(LIE, RSML[[i]]$lie)
          length1<-length(RSML[[i]]$lie)
          if (length1>1){
            num<-c(1:length1)
            filenameslie[(length(filenameslie)+1):(length(filenameslie)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
          if (length1==1){
            filenameslie[(length(filenameslie)+1)]<-filenamesrsml[i]}}
        
        #Unit conversion rsml
        
        cunit1<-vector(length=length(res1))
        
        for (i in 1:length(res1)){
          
          if (unitlength=="cm"){
            
            if (unitlength1[i]=="pixel") {
              cunit1[i]<-1
              message(paste("Unit in ", filenameslie[i], " is pixel. Unitlength not used and results expressed in pixels", sep=""))}
            if (unitlength1[i]=="m") {cunit1[i]<-100/res1[i]}
            if (unitlength1[i]=="cm") {cunit1[i]<-1/res1[i]}
            if (unitlength1[i]=="mm") {cunit1[i]<-1/res1[i]/10}
            if (unitlength1[i]=="um") {cunit1[i]<-1/res1[i]/10000}
            if (unitlength1[i]=="nm") {cunit1[i]<-1/res1[i]/10000000}
            if (unitlength1[i]=="inch") {cunit1[i]<-1/res1[i]*cm(1)}}
          
          if (unitlength=="mm"){
            if (unitlength1[i]=="pixel") {
              cunit1[i]<-1
              message(paste("Unit in ", filenameslie[i], " is pixel. Unitlength not used and results expressed in pixels", sep=""))}
            if (unitlength1[i]=="m") {cunit1[i]<-1/res1[i]*1000}
            if (unitlength1[i]=="cm") {cunit1[i]<-1/res1[i]*10}
            if (unitlength1[i]=="mm") {cunit1[i]<-1/res1[i]}
            if (unitlength1[i]=="um") {cunit1[i]<-1/res1[i]/1000}
            if (unitlength1[i]=="nm") {cunit1[i]<-1/res1[i]/1000000}
            if (unitlength1[i]=="inch") {cunit1[i]<-1/res1[i]*cm(1)*10}}
          
          if (unitlength=="px"){cunit1[i]<-1}}}
      
      else { # DART and rsml files
        
        res1<-rep(res, length(filenames.lie))
        unitlength1<-rep(unitlength, length(filenames.lie))
        
        LIE<-lapply(paste(path.lie, "/", filenames.lie, sep=""), read.table, header=TRUE)
        RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect)
        for (i in 1:length(RSML)){
          res1<-append(res1, rep(as.numeric(RSML[[i]]$resolution), length(RSML[[i]]$lie)))
          unitlength1<-append(unitlength1, rep(as.character(RSML[[i]]$length), length(RSML[[i]]$lie)))
          LIE<-append(LIE, RSML[[i]]$lie)
          length1<-length(RSML[[i]]$lie)
          if (length1>1){
            num<-c(1:length1)
            filenameslie[(length(filenameslie)+1):(length(filenameslie)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
          if (length1==1){
            filenameslie[(length(filenameslie)+1)]<-filenamesrsml[i]}}
        
        #Unit conversion DART and RSML
        cunit1<-vector(length=length(res1))
        
        if (unitlength=="mm") {cunit1[1:length(filenames.lie)]<-(10*cm(1)/res)}
        if (unitlength=="cm") {cunit1[1:length(filenames.lie)]<-cm(1)/res}
        if (unitlength=="px") {cunit1[1:length(filenames.lie)]<-1}
        
        for (i in (length(filenames.lie)+1):length(res1)){
          
          if (unitlength=="cm"){
            
            if (unitlength1[i]=="pixel") {
              cunit[i]<-1
              message(paste("Unit in ", filenameslie[i], " is pixel. Unitlength not used and results expressed in pixels", sep=""))}
            if (unitlength1[i]=="m") {cunit1[i]<-100/res1[i]}
            if (unitlength1[i]=="cm") {cunit1[i]<-1/res1[i]}
            if (unitlength1[i]=="mm") {cunit1[i]<-1/res1[i]/10}
            if (unitlength1[i]=="um") {cunit1[i]<-1/res1[i]/10000}
            if (unitlength1[i]=="nm") {cunit1[i]<-1/res1[i]/10000000}
            if (unitlength1[i]=="inch") {cunit1[i]<-1/res1[i]*cm(1)}}
          
          if (unitlength=="mm"){
            if (unitlength1[i]=="pixel") {
              cunit[i]<-1
              message(paste("Unit in ", filenameslie[i], " is pixel. Unitlength not used and results expressed in pixels", sep=""))}
            if (unitlength1[i]=="m") {cunit1[i]<-1/res1[i]*1000}
            if (unitlength1[i]=="cm") {cunit1[i]<-1/res1[i]*10}
            if (unitlength1[i]=="mm") {cunit1[i]<-1/res1[i]}
            if (unitlength1[i]=="um") {cunit1[i]<-1/res1[i]/1000}
            if (unitlength1[i]=="nm") {cunit1[i]<-1/res1[i]/1000000}
            if (unitlength1[i]=="inch") {cunit1[i]<-1/res1[i]*cm(1)*10}}
          
          if (unitlength=="px"){cunit1[i]<-1}}}}
    
    # Unit conversion and rotation for 2D root systems
    
    if (unitangle=="r") {
      cunitangle<-1
      rotation<-rotation}
    if (unitangle=="d") {
      cunitangle<-180/pi
      rotation<-rotation*(1/cunitangle)}
    
    rot.matrix<-matrix(c(cos(rotation), sin(rotation), -sin(rotation), cos(rotation)), nrow=2, ncol=2)
    
    for (i in 1:length(LIE)){
      
      if (is.null(inputrsml)==FALSE){cunit<-cunit1[i]}
      
      if (!("Z" %in% colnames(LIE[[i]]))){
      LIE[[i]]$X<-LIE[[i]]$X*cunit
      LIE[[i]]$Y<-LIE[[i]]$Y*cunit
      newcoord<-rot.matrix%*%t(as.matrix(data.frame(LIE[[i]]$X, LIE[[i]]$Y)))
      LIE[[i]]$X<-newcoord[1,]
      LIE[[i]]$Y<-newcoord[2,]}
      else {
        LIE[[i]]$X<-LIE[[i]]$X*cunit
        LIE[[i]]$Y<-LIE[[i]]$Y*cunit
        LIE[[i]]$Z<-LIE[[i]]$Z*cunit}}
    
    # Drawing the root system architecture for each DART/RSML file
  
    for (i in 1:length(LIE)){
      
      a<-max(LIE[[i]]$Date)
      message(paste("Number of observation dates for ", filenameslie[i], ": ", a, sep=""))
      
    if (is.null(numdate)==TRUE){
      if (length(coldate)>max(LIE[[i]]$Date)){message(paste("Note: The number of colors in coldate is greater than the number of observation dates in ", filenameslie[i], sep=""))}
      if (length(coldate)<max(LIE[[i]]$Date)){message(paste("Note: The number of colors in coldate is lower than the number of observation dates in ", filenameslie[i], sep=""))}
      coldate1<-rep(coldate, length.out=max(LIE[[i]]$Date))
      LIE[[i]]$X<-LIE[[i]]$X-min(LIE[[i]]$X)
      LIE[[i]]$Y<-LIE[[i]]$Y-min(LIE[[i]]$Y)
      if ("Z" %in% colnames(LIE[[i]])) {LIE[[i]]$Z<-LIE[[i]]$Z-min(LIE[[i]]$Z)}
      minx<-min(LIE[[i]]$X)
      maxx<-max(LIE[[i]]$X)
      miny<-min(LIE[[i]]$Y)
      maxy<-max(LIE[[i]]$Y)
      if ("Z" %in% colnames(LIE[[i]])){
        minz<-min(LIE[[i]]$Z)
        maxz<-max(LIE[[i]]$Z)}
      if (is.null(main)==TRUE){main1<-filenameslie[i]} else {main1<-main}
      if (is.null(xlim)==TRUE){xlim1<-c(minx,maxx)} else {xlim1<-xlim}
      if (is.null(ylim)==TRUE){ylim1<-c(maxy,miny)} else {ylim1<-ylim}
      if ("Z" %in% colnames(LIE[[i]]) & is.null(zlim)==TRUE){zlim1<-c(minz,maxz)} else {zlim1<-zlim}
      if (is.null(xlab)==TRUE){xlab1<-paste("X (", unitlength, ")", sep="")} else {xlab1<-xlab}
      if (is.null(ylab)==TRUE){ylab1<-paste("Y (", unitlength, ")", sep="")} else {ylab1<-ylab}
      if ("Z" %in% colnames(LIE[[i]]) & is.null(zlab)==TRUE){zlab1<-paste("Z (", unitlength, ")", sep="")} else {zlab1<-zlab}
      
      if ("Z" %in% colnames(LIE[[i]])){
        
        if (is.null(twod)==TRUE){
        
        LIE[[i]]$col<-LIE[[i]]$Date
        for (l in 1:length(coldate1)){LIE[[i]]$col<-replace(LIE[[i]]$col, LIE[[i]]$col==l, coldate1[l])}
        root<-sum(LIE[[i]]$Suiv==0)
        end<-which(LIE[[i]]$Suiv==0)
        open3d()
        plot3d(x=LIE[[i]]$X[1], y=LIE[[i]]$Y[1], z=LIE[[i]]$Z[1], type="n", xlim=xlim1, ylim=ylim1, zlim=zlim1, main=main1, ylab=ylab1, xlab=xlab1, zlab=zlab1,...)
        
        for (k in 1:root){
          if (k==1) {
            dataroot<-as.matrix(LIE[[i]][k:end[k],c(7:9,14)])} 
          else {
            if (LIE[[i]]$Prec[end[k-1]+1]!=0){
              dataroot<-as.matrix(LIE[[i]][(end[k-1]+1):end[k],c(7:9,14)])
              dataroot<-rbind(LIE[[i]][LIE[[i]]$Num==LIE[[i]]$Prec[end[k-1]+1],c(7:9,14)], dataroot)}
            else{
              dataroot<-as.matrix(LIE[[i]][(end[k-1]+1):end[k],c(7:9,14)])}}
  
          lines3d(dataroot[,1:3], col=dataroot[,4], smooth=FALSE, ...)}}
        
        else {
          
          if (all.equal(twod, c("x", "y"))==TRUE){
          
            plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
            for (k in 1:nrow(LIE[[i]])){
              a<-LIE[[i]]$Prec[k]
              b<-LIE[[i]]$Date[k]
              if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=coldate1[b],...)}}}
          
          if (all.equal(twod, c("x", "z"))==TRUE){
            
            plot(LIE[[i]]$X[1], LIE[[i]]$Z[1], type="n", xlim=xlim1, ylim=zlim1, main=main1, ylab=zlab1, xlab=xlab1,...)
            for (k in 1:nrow(LIE[[i]])){
              a<-LIE[[i]]$Prec[k]
              b<-LIE[[i]]$Date[k]
              if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Z[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Z[k], col=coldate1[b],...)}}}
          
          if (all.equal(twod, c("y", "z"))==TRUE){
            
            plot(LIE[[i]]$Z[1], LIE[[i]]$Y[1], type="n", xlim=zlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=zlab1,...)
            for (k in 1:nrow(LIE[[i]])){
              a<-LIE[[i]]$Prec[k]
              b<-LIE[[i]]$Date[k]
              if (a!=0){segments(x0=LIE[[i]]$Z[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$Z[k],y1=LIE[[i]]$Y[k], col=coldate1[b],...)}}}}}
      
      else {
            plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
            for (k in 1:nrow(LIE[[i]])){
              a<-LIE[[i]]$Prec[k]
              b<-LIE[[i]]$Date[k]
              if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=coldate1[b],...)}}}}
    
    if (is.null(numdate)==FALSE) {
      
      for (l in 1:length(numdate)){if (numdate[l]>max(LIE[[i]]$Date)){message(paste("Note: numdate contains a numerical value greater than the number of observation dates in ", filenameslie[i], ". The root system at the last observation date will be plotted.", sep=""))}}
      if (length(coldate)>max(numdate)){message("Note: The number of colors in coldate is greater than max(numdate)")}
      if (length(coldate)<max(numdate)){message("Note: The number of colors in coldate is lower than max(numdate)")}
      coldate1<-rep(coldate, length.out=max(numdate))
      
      if (finalscale==TRUE){
        LIE[[i]]$X<-LIE[[i]]$X-min(LIE[[i]]$X)
        LIE[[i]]$Y<-LIE[[i]]$Y-min(LIE[[i]]$Y)
        if ("Z" %in% colnames(LIE[[i]])) {LIE[[i]]$Z<-LIE[[i]]$Z-min(LIE[[i]]$Z)}
        minx<-min(LIE[[i]]$X)
        maxx<-max(LIE[[i]]$X)
        miny<-min(LIE[[i]]$Y)
        maxy<-max(LIE[[i]]$Y)
        if ("Z" %in% colnames(LIE[[i]])){
          minz<-min(LIE[[i]]$Z)
          maxz<-max(LIE[[i]]$Z)}
        for (j in 1:length(numdate)){
          if (is.null(main)==TRUE){main1<-paste(filenameslie[i], "-numdate=", numdate[j], sep="")} else {main1<-main}
          if (is.null(xlim)==TRUE){xlim1<-c(minx,maxx)} else {xlim1<-xlim}
          if (is.null(ylim)==TRUE){ylim1<-c(maxy,miny)} else {ylim1<-ylim}
          if ("Z" %in% colnames(LIE[[i]]) & is.null(zlim)==TRUE){zlim1<-c(minz,maxz)} else {zlim1<-zlim}
          if (is.null(xlab)==TRUE){xlab1<-paste("X (", unitlength, ")", sep="")} else {xlab1<-xlab}
          if (is.null(ylab)==TRUE){ylab1<-paste("Y (", unitlength, ")", sep="")} else {ylab1<-ylab}
          if ("Z" %in% colnames(LIE[[i]]) & is.null(zlab)==TRUE){zlab1<-paste("Z (", unitlength, ")", sep="")} else {zlab1<-zlab}
          
          if ("Z" %in% colnames(LIE[[i]])){
            
            if (is.null(twod)==TRUE){
            
            LIE[[i]]$col<-LIE[[i]]$Date
            for (l in 1:length(coldate1)){LIE[[i]]$col<-replace(LIE[[i]]$col, LIE[[i]]$col==l, coldate1[l])}
            dataroot1<-LIE[[i]][LIE[[i]]$Date<=numdate[j],]
            root<-sum(dataroot1$Bran=="true")
            if (root>1) {
              end<-which(dataroot1$Bran=="true")-1
              end<-end[-1]
              end[length(end)+1]<-nrow(dataroot1)}
            else {end<-nrow(dataroot1)}
            open3d()
            plot3d(x=LIE[[i]]$X[1], y=LIE[[i]]$Y[1], z=LIE[[i]]$Z[1], type="n", xlim=xlim1, ylim=ylim1, zlim=zlim1, main=main1, ylab=ylab1, xlab=xlab1, zlab=zlab1,...)
            for (k in 1:root){
                if (k==1) {
                  dataroot<-as.matrix(dataroot1[k:end[k],c(7:9,14)])} 
                else {
                  if (dataroot1$Prec[end[k-1]+1]!=0){
                    dataroot<-as.matrix(dataroot1[(end[k-1]+1):end[k],c(7:9,14)])
                    dataroot<-rbind(dataroot1[dataroot1$Num==dataroot1$Prec[end[k-1]+1],c(7:9,14)], dataroot)}
                  else{
                    dataroot<-as.matrix(dataroot1[(end[k-1]+1):end[k],c(7:9,14)])}}
                
              lines3d(dataroot[,1:3], col=dataroot[,4], smooth=FALSE, ...)}}
            
            else {
              
              if (all.equal(twod, c("x", "y"))==TRUE){
                
                plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
                for (k in 1:nrow(LIE[[i]])){
                  if (LIE[[i]]$Date[k]<=numdate[j]){
                  a<-LIE[[i]]$Prec[k]
                  b<-LIE[[i]]$Date[k]
                  if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=coldate1[b],...)}}}}
              
              if (all.equal(twod, c("x", "z"))==TRUE){
                
                plot(LIE[[i]]$X[1], LIE[[i]]$Z[1], type="n", xlim=xlim1, ylim=zlim1, main=main1, ylab=zlab1, xlab=xlab1,...)
                for (k in 1:nrow(LIE[[i]])){
                  if (LIE[[i]]$Date[k]<=numdate[j]){
                  a<-LIE[[i]]$Prec[k]
                  b<-LIE[[i]]$Date[k]
                  if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Z[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Z[k], col=coldate1[b],...)}}}}
              
              if (all.equal(twod, c("y", "z"))==TRUE){
                
                plot(LIE[[i]]$Z[1], LIE[[i]]$Y[1], type="n", xlim=zlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=zlab1,...)
                for (k in 1:nrow(LIE[[i]])){
                  if (LIE[[i]]$Date[k]<=numdate[j]){
                  a<-LIE[[i]]$Prec[k]
                  b<-LIE[[i]]$Date[k]
                  if (a!=0){segments(x0=LIE[[i]]$Z[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$Z[k],y1=LIE[[i]]$Y[k], col=coldate1[b],...)}}}}}} 
              
              else {
                plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
                for (k in 1:nrow(LIE[[i]])){
                  if (LIE[[i]]$Date[k]<=numdate[j]){
                  a<-LIE[[i]]$Prec[k]
                  b<-LIE[[i]]$Date[k]
                  if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=coldate1[b],...)}}}}}}
      
      if (finalscale==FALSE){
        for (j in 1:length(numdate)){
          Date=NULL
          Num=NULL
          Suiv=NULL
          Prec=NULL
          X=NULL
          Y=NULL
          Z=NULL
          if ("Z" %in% colnames(LIE[[i]])) {subset.LIE<-subset(LIE[[i]], Date<=numdate[j], select=c(Num, Date, Suiv, Prec, X, Y, Z))} else {subset.LIE<-subset(LIE[[i]], Date<=numdate[j], select=c(Num, Date, Suiv, Prec, X, Y))}
          LIE[[i]]$X<-LIE[[i]]$X-min(subset.LIE$X)
          LIE[[i]]$Y<-LIE[[i]]$Y-min(subset.LIE$Y)
          if ("Z" %in% colnames(LIE[[i]])) {LIE[[i]]$Z<-LIE[[i]]$Z-min(subset.LIE$Z)}
          subset.LIE$X<-subset.LIE$X-min(subset.LIE$X)
          subset.LIE$Y<-subset.LIE$Y-min(subset.LIE$Y)
          if ("Z" %in% colnames(LIE[[i]])) {subset.LIE$Z<-subset.LIE$Z-min(subset.LIE$Z)}
          minx<-min(subset.LIE$X)
          maxx<-max(subset.LIE$X)
          miny<-min(subset.LIE$Y)
          maxy<-max(subset.LIE$Y)
          if ("Z" %in% colnames(LIE[[i]])){
            minz<-min(subset.LIE$Z)
            maxz<-max(subset.LIE$Z)}
          if (is.null(main)==TRUE){main1<-paste(filenameslie[i], "-numdate=", numdate[j], sep="")} else {main1<-main}
          if (is.null(xlim)==TRUE){xlim1<-c(minx,maxx)} else {xlim1<-xlim}
          if (is.null(ylim)==TRUE){ylim1<-c(maxy,miny)} else {ylim1<-ylim}
          if ("Z" %in% colnames(LIE[[i]]) & is.null(zlim)==TRUE){zlim1<-c(minz,maxz)} else {zlim1<-zlim}
          if (is.null(xlab)==TRUE){xlab1<-paste("X (", unitlength, ")", sep="")} else {xlab1<-xlab}
          if (is.null(ylab)==TRUE){ylab1<-paste("Y (", unitlength, ")", sep="")} else {ylab1<-ylab}
          if ("Z" %in% colnames(LIE[[i]]) & is.null(zlab)==TRUE){zlab1<-paste("Z (", unitlength, ")", sep="")} else {zlab1<-zlab}
          
          if ("Z" %in% colnames(LIE[[i]])){
            
            if (is.null(twod)==TRUE){
            
            LIE[[i]]$col<-LIE[[i]]$Date
            for (l in 1:length(coldate1)){LIE[[i]]$col<-replace(LIE[[i]]$col, LIE[[i]]$col==l, coldate1[l])}
            dataroot1<-LIE[[i]][LIE[[i]]$Date<=numdate[j],]
            root<-sum(dataroot1$Bran=="true")
            if (root>1) {
              end<-which(dataroot1$Bran=="true")-1
              end<-end[-1]
              end[length(end)+1]<-nrow(dataroot1)}
            else {end<-nrow(dataroot1)}
            open3d()
            plot3d(x=LIE[[i]]$X[1], y=LIE[[i]]$Y[1], z=LIE[[i]]$Z[1], type="n", xlim=xlim1, ylim=ylim1, zlim=zlim1, main=main1, ylab=ylab1, xlab=xlab1, zlab=zlab1,...)
            
            for (k in 1:root){
              if (k==1) {
                dataroot<-as.matrix(dataroot1[k:end[k],c(7:9,14)])} 
              else {
                if (dataroot1$Prec[end[k-1]+1]!=0){
                  dataroot<-as.matrix(dataroot1[(end[k-1]+1):end[k],c(7:9,14)])
                  dataroot<-rbind(dataroot1[dataroot1$Num==dataroot1$Prec[end[k-1]+1],c(7:9,14)], dataroot)}
                else{
                  dataroot<-as.matrix(dataroot1[(end[k-1]+1):end[k],c(7:9,14)])}}
              
              lines3d(dataroot[,1:3], col=dataroot[,4], smooth=FALSE, ...)}}
            
            else {
              
              if (all.equal(twod, c("x", "y"))==TRUE){
                
                plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
                for (k in 1:nrow(LIE[[i]])){
                  if (LIE[[i]]$Date[k]<=numdate[j]){
                    a<-LIE[[i]]$Prec[k]
                    b<-LIE[[i]]$Date[k]
                    if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=coldate1[b],...)}}}}
              
              if (all.equal(twod, c("x", "z"))==TRUE){
                
                plot(LIE[[i]]$X[1], LIE[[i]]$Z[1], type="n", xlim=xlim1, ylim=zlim1, main=main1, ylab=zlab1, xlab=xlab1,...)
                for (k in 1:nrow(LIE[[i]])){
                  if (LIE[[i]]$Date[k]<=numdate[j]){
                    a<-LIE[[i]]$Prec[k]
                    b<-LIE[[i]]$Date[k]
                    if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Z[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Z[k], col=coldate1[b],...)}}}}
              
              if (all.equal(twod, c("y", "z"))==TRUE){
                
                plot(LIE[[i]]$Z[1], LIE[[i]]$Y[1], type="n", xlim=zlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=zlab1,...)
                for (k in 1:nrow(LIE[[i]])){
                  if (LIE[[i]]$Date[k]<=numdate[j]){
                    a<-LIE[[i]]$Prec[k]
                    b<-LIE[[i]]$Date[k]
                    if (a!=0){segments(x0=LIE[[i]]$Z[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$Z[k],y1=LIE[[i]]$Y[k], col=coldate1[b],...)}}}}}} 
          
          else {
                  plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
                  for (k in 1:nrow(LIE[[i]])){
                    if (LIE[[i]]$Date[k]<=numdate[j]){
                      a<-LIE[[i]]$Prec[k]
                      b<-LIE[[i]]$Date[k]
                      if (a!=0){segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=coldate1[b],...)}}}}}}}}}