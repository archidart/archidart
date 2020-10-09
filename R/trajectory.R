trajectory<-function(inputrac=NULL, inputlie=NULL, inputtps=NULL, inputrsml=NULL, res=NULL, unitlength="px", unitangle="d", rotation=0, l.brangle, l.curv, l.tipangle, rsml.date=NULL, vertical3d="y", plot=NULL, twod=NULL, colangle=NULL, export.colors=FALSE, BRscale=NULL, main=NULL, xlim=NULL, ylim=NULL, zlim=NULL, xlab=NULL, ylab=NULL, zlab=NULL, ...){
  
  # Errors interception
  
  if (is.null(inputrac)==TRUE & is.null(inputlie)==TRUE & is.null(inputtps)==TRUE & is.null(inputrsml)==TRUE){stop("inputrac/inpulie/inputtps and/or inputrsml must be provided")}
  
  if (is.null(inputrac)==FALSE) {if (mode(inputrac)!="character"){stop("mode(inputrac) must be character")}}
  
  if (is.null(inputlie)==FALSE) {if (mode(inputlie)!="character"){stop("mode(inputlie) must be character")}}
  
  if (is.null(inputtps)==FALSE) {if (mode(inputtps)!="character"){stop("mode(inputtps) must be character")}}
  
  if (is.null(inputrsml)==FALSE) {if (mode(inputrsml)!="character"){stop("mode(inputrsml) must be character")}}
  
  if (is.null(inputrac)==FALSE|is.null(inputtps)==FALSE|is.null(inputlie)==FALSE){
    if (is.null(inputrac)==TRUE|is.null(inputtps)==TRUE|is.null(inputlie)==TRUE){stop("If inputrac/inputlie/inputtps is not NULL, inputrac/inputlie/inputtps must be provided")}}
  
  if (is.null(inputrac)==FALSE & is.null(res)==TRUE & unitlength!="px"){stop("If unitlength is not px, res must be specified")}
  if (is.null(res)==FALSE){
    if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
    if (res<=0){stop("res must be a positive value")}}
  
  if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
  if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
  
  if (mode(unitangle)!="character"){stop("mode(unitangle) must be character")}
  if(unitangle=="d"|unitangle=="r") {} else {stop("unitangle must be either d (degrees) or r (radians)")}
  
  if (mode(rotation)!="numeric"){stop("mode(rotation) must be numeric")}
  if (rotation<0){stop("rotation must be a positive value")}
  
  if (mode(l.brangle)!="numeric"){stop("mode(l.brangle) must be numeric")}
  if (l.brangle<=0) {stop("l.brangle must be a positive value")}
  
  if (mode(l.curv)!="numeric"){stop("mode(l.curv) must be numeric")}
  if (l.curv<=0) {stop("l.curv must be a positive value")}
  
  if (mode(l.tipangle)!="numeric"){stop("mode(l.tipangle) must be numeric")}
  if (l.tipangle<=0) {stop("l.tipangle must be a positive value")}
  
  if (is.null(rsml.date)==FALSE){
    if (is.character(rsml.date)==TRUE|is.numeric(rsml.date)==TRUE){} else {stop("If rsml.date is not NULL, rsml.date must be a character string or a positive numeric value")}
    if (is.numeric(rsml.date)==TRUE){if (rsml.date<=0|length(rsml.date)>1){stop("If mode(rsml.date) is numeric, rsml.date must be a single positive value")}}}
  
  if (vertical3d=="x"|vertical3d=="y"|vertical3d=="z") {} else {stop("vertical3d must be x, y, or z")}
  
  if (is.null(plot)==FALSE){if (plot=="branching"|plot=="direction") {} else {stop("If plot is not NULL, plot most be branching or direction")}}
  
  if (is.null(export.colors)==FALSE) {if (mode(export.colors)!="logical"){stop("mode(export.colors) must be logical")}}
  
  if (is.null(plot)==FALSE & is.null(colangle)==TRUE){stop("If plot is not NULL, colangle must be specified")}
  
  if (is.null(BRscale)==FALSE){
    if (mode(BRscale)!="numeric") {stop("mode(BRscale) must be numeric")}
    if (length(BRscale)!=2|diff(BRscale)==0){stop("length(BRscale) must be equal to 2: c(min, max)")}
    if (BRscale[1]<0|BRscale[2]<0){stop("BRscale must be a vector of positive values")}}
  
  if (is.null(twod)==FALSE){
    if (mode(twod)!="character"){stop("mode(twod) must be character")}
    if (length(twod)!=2){stop("twod must be a vector of 2 character elements")}
    twod<-sort(twod)
    if (all.equal(twod, c("x", "y"))==TRUE|all.equal(twod, c("x", "z"))==TRUE|all.equal(twod, c("y", "z"))==TRUE) {} else {stop("twod must be c(x,y), c(x,z), or c(y,z)")}}
  
  # Reading of DART and rsml files
  
  if (is.null(inputrac)==FALSE){
    filenames.rac<-mixedsort(list.files(path=inputrac, pattern="\\.rac$"))
    path.rac<-rep(inputrac, length.out=length(filenames.rac))
    filenamesrac<-sub(x=filenames.rac, pattern="\\.rac$", replacement="")
    message(paste("Number of DART rac files in inputrac:", length(filenames.rac), sep=" "))}
  
  if (is.null(inputtps)==FALSE){
    filenames.tps<-mixedsort(list.files(path=inputtps, pattern="\\.tps$"))
    path.tps<-rep(inputtps, length.out=length(filenames.tps))
    filenamestps<-sub(x=filenames.tps, pattern="\\.tps$", replacement="")
    message(paste("Number of DART tps files in inputtps:", length(filenames.tps), sep=" "))}
  
  if (is.null(inputlie)==FALSE){
    filenames.lie<-mixedsort(list.files(path=inputlie, pattern="\\.lie$"))
    path.lie<-rep(inputlie, length.out=length(filenames.lie))
    filenameslie<-sub(x=filenames.lie, pattern="\\.lie$", replacement="")
    message(paste("Number of DART lie files in inputlie:", length(filenames.lie), sep=" "))}
  
  if (is.null(inputrsml)==FALSE){
    filenames.rsml<-mixedsort(list.files(path=inputrsml, pattern="\\.rsml$"))
    path.rsml<-rep(inputrsml, length.out=length(filenames.rsml))
    filenamesrsml<-sub(x=filenames.rsml, pattern="\\.rsml$", replacement="")
    message(paste("Number of rsml files in inputrsml:", length(filenames.rsml), sep=" "))}
  
  if (is.null(inputrsml)==TRUE){
    if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
    if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}
    if (length(filenames.lie)==0){stop("There is no lie file in inputlie")}}
  else {
    if (is.null(inputrac)==TRUE){if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}
    else{
      if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
      if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}
      if (length(filenames.lie)==0){stop("There is no lie file in inputlie")}
      if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}}
  
  if (is.null(inputrsml)==TRUE){ # Only DART files
    
    LIE<-lapply(paste(path.lie, "/", filenames.lie, sep=""), read.table, header=TRUE)
    
    DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
    for (i in 1:length(DATA)) {
      colnames(DATA[[i]])<-c()
      colnames(DATA[[i]])[1]<-"Root"
      colnames(DATA[[i]])[2]<-"Mother"
      colnames(DATA[[i]])[3]<-"Ord"
      colnames(DATA[[i]])[4]<-"DBase"
      colnames(DATA[[i]])[5]<-"DApp"
      for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
    
    TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE)
    
    #Unit conversion DART files
    if (unitlength=="mm") {cunit<-(10*cm(1)/res)}
    if (unitlength=="cm") {cunit<-(cm(1)/res)}
    if (unitlength=="px") {cunit<-1}
  
    if (length(LIE)!=length(DATA)) {stop("The number of rac files in inputrac and lie files in inputlie must be equal")}
    else {
      for (i in 1:length(DATA)) {if(filenamesrac[i]!=filenameslie[i]) {stop("Input rac files and their corresponding lie files must have the same name")}}}  	
    
    if (length(TIME)==1) {
      for (i in 1:length(DATA)) {if(length(TIME[[1]]$Date)!=(ncol(DATA[[i]])-5)){stop("The number of observation dates between corresponding rac et tps files must be equal")}}}
    else {
      if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac/lie files in inputrac/inputlie and tps files in inputtps must be equal")}
      else {
        for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac/lie files and their corresponding tps files must have the same name")}}
        for (i in 1:length(DATA)) {if (length(TIME[[i]]$Date)!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}} 
  
  else {
    
    if (is.null(inputrac)==TRUE){ # Only rsml files
      
      LIE<-list()
      DATA<-list()
      TIME<-list()
      res1<-c()
      unitlength1<-c()
      filenameslie<-c()
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=TRUE)
      for (i in 1:length(RSML)){
        res1<-append(res1, rep(as.numeric(RSML[[i]]$resolution, length(RSML[[i]]$lie))))
        unitlength1<-append(unitlength1, rep(as.character(RSML[[i]]$length), length(RSML[[i]]$lie)))
        DATA<-append(DATA, RSML[[i]]$rac)
        LIE<-append(LIE, RSML[[i]]$lie)
        TIME<-append(TIME, RSML[[i]]$tps)
        length1<-length(RSML[[i]]$rac)
        if (length1>1){
          num<-c(1:length1)
          filenameslie[(length(filenameslie)+1):(length(filenameslie)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenameslie[(length(filenameslie)+1)]<-filenamesrsml[i]}}
      
      #Unit conversion RSML
      cunit1<-vector(length=length(res1))
      
      for (i in 1:length(res1)){
        
        if (unitlength=="cm"){
          
          if (unitlength1[i]=="pixel") {
            cunit1[i]<-1
            message(paste("Unit in ", filenamesrac[i], " is pixel. Unitlength not used and results expressed in pixels", sep=""))}
          if (unitlength1[i]=="m") {cunit1[i]<-100/res1[i]}
          if (unitlength1[i]=="cm") {cunit1[i]<-1/res1[i]}
          if (unitlength1[i]=="mm") {cunit1[i]<-1/res1[i]/10}
          if (unitlength1[i]=="um") {cunit1[i]<-1/res1[i]/10000}
          if (unitlength1[i]=="nm") {cunit1[i]<-1/res1[i]/10000000}
          if (unitlength1[i]=="inch") {cunit1[i]<-1/res1[i]*cm(1)}}
        
        if (unitlength=="mm"){
          if (unitlength1[i]=="pixel") {
            cunit1[i]<-1
            message(paste("Unit in ", filenamesrac[i], " is pixel. Unitlength not used and results expressed in pixels", sep=""))}
          if (unitlength1[i]=="m") {cunit1[i]<-1/res1[i]*1000}
          if (unitlength1[i]=="cm") {cunit1[i]<-1/res1[i]*10}
          if (unitlength1[i]=="mm") {cunit1[i]<-1/res1[i]}
          if (unitlength1[i]=="um") {cunit1[i]<-1/res1[i]/1000}
          if (unitlength1[i]=="nm") {cunit1[i]<-1/res1[i]/1000000}
          if (unitlength1[i]=="inch") {cunit1[i]<-1/res1[i]*cm(1)*10}}
        
        if (unitlength=="px"){cunit1[i]<-1}}}
    
    else { # DART and rsml files
      
      LIE<-lapply(paste(path.lie, "/", filenames.lie, sep=""), read.table, header=TRUE)
      
      DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
      for (i in 1:length(DATA)) {
        colnames(DATA[[i]])<-c()
        colnames(DATA[[i]])[1]<-"Root"
        colnames(DATA[[i]])[2]<-"Mother"
        colnames(DATA[[i]])[3]<-"Ord"
        colnames(DATA[[i]])[4]<-"DBase"
        colnames(DATA[[i]])[5]<-"DApp"
        for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
      
      TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE)
      
      if (length(LIE)!=length(DATA)) {stop("The number of rac files in inputrac and lie files in inputlie must be equal")}
      else {
        for (i in 1:length(DATA)) {if(filenamesrac[i]!=filenameslie[i]) {stop("Input rac files and their corresponding lie files must have the same name")}}}  	
      
      if (length(TIME)==1) {
        for (i in 1:length(DATA)) {if(length(TIME[[1]]$Date)!=(ncol(DATA[[i]])-5)){stop("The number of observation dates between corresponding rac et tps files must be equal")}}
        TIME[1:length(DATA)]<-TIME[1]}
      else {
        if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac/lie files in inputrac/inputlie and tps files in inputtps must be equal")}
        else {
          for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac/lie files and their corresponding tps files must have the same name")}}
          for (i in 1:length(DATA)) {if (length(TIME[[i]]$Date)!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}
      
      res1<-rep(res, length(filenames.rac))
      unitlength1<-rep(unitlength, length(filenames.rac))
      
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=TRUE)
      
      for (i in 1:length(RSML)){
        res1<-append(res1, rep(as.numeric(RSML[[i]]$resolution, length(RSML[[i]]$lie))))
        unitlength1<-append(unitlength1, rep(as.character(RSML[[i]]$length), length(RSML[[i]]$lie)))
        DATA<-append(DATA, RSML[[i]]$rac)
        LIE<-append(LIE, RSML[[i]]$lie)
        TIME<-append(TIME, RSML[[i]]$tps)
        length1<-length(RSML[[i]]$rac)
        if (length1>1){
          num<-c(1:length1)
          filenameslie[(length(filenameslie)+1):(length(filenameslie)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenameslie[(length(filenameslie)+1)]<-filenamesrsml[i]}}
      
      #Unit conversion DART and RSML
      cunit1<-vector(length=length(res1))
      
      if (unitlength=="mm") {cunit1[1:length(filenames.rac)]<-(10*cm(1)/res)}
      if (unitlength=="cm") {cunit1[1:length(filenames.rac)]<-cm(1)/res}
      if (unitlength=="px") {cunit1[1:length(filenames.rac)]<-1}
      
      for (i in (length(filenames.rac)+1):length(res1)){
        
        if (unitlength=="cm"){
          
          if (unitlength1[i]=="pixel") {
            cunit[i]<-1
            message(paste("Unit in ", filenamesrac[i], " is pixel. Unitlength not used and results expressed in pixels", sep=""))}
          if (unitlength1[i]=="m") {cunit1[i]<-100/res1[i]}
          if (unitlength1[i]=="cm") {cunit1[i]<-1/res1[i]}
          if (unitlength1[i]=="mm") {cunit1[i]<-1/res1[i]/10}
          if (unitlength1[i]=="um") {cunit1[i]<-1/res1[i]/10000}
          if (unitlength1[i]=="nm") {cunit1[i]<-1/res1[i]/10000000}
          if (unitlength1[i]=="inch") {cunit1[i]<-1/res1[i]*cm(1)}}
        
        if (unitlength=="mm"){
          if (unitlength1[i]=="pixel") {
            cunit[i]<-1
            message(paste("Unit in ", filenamesrac[i], " is pixel. Unitlength not used and results expressed in pixels", sep=""))}
          if (unitlength1[i]=="m") {cunit1[i]<-1/res1[i]*1000}
          if (unitlength1[i]=="cm") {cunit1[i]<-1/res1[i]*10}
          if (unitlength1[i]=="mm") {cunit1[i]<-1/res1[i]}
          if (unitlength1[i]=="um") {cunit1[i]<-1/res1[i]/1000}
          if (unitlength1[i]=="nm") {cunit1[i]<-1/res1[i]/1000000}
          if (unitlength1[i]=="inch") {cunit1[i]<-1/res1[i]*cm(1)*10}}
        
        if (unitlength=="px"){cunit1[i]<-1}}}}
  
  # Unit conversion and rotation
  
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
    
  # Vertical direction vector
  
  if ("Z" %in% colnames(LIE[[i]])) {
    
    if (vertical3d=="x") {
      if (max(LIE[[i]]$X)+min(LIE[[i]]$X)>0) {dirvert<-c(1,0,0)} else {dirvert<-c(-1,0,0)}}
    if (vertical3d=="y") {
      if (max(LIE[[i]]$Y)+min(LIE[[i]]$Y)>0) {dirvert<-c(0,1,0)} else {dirvert<-c(0,-1,0)}}
    if (vertical3d=="z") {
      if (max(LIE[[i]]$Z)+min(LIE[[i]]$Z)>0) {dirvert<-c(0,0,1)} else {dirvert<-c(0,0,-1)}}}
  
  # Creating vectors and matrices for root architecture parameters calculation
  
    filenames<-c()
    finallength.lie<-list()
    Root<-c()
    Mother<-c()
    Ord<-c()
    DBase<-c()
    DApp<-c()
    Finalrootlength<-c()
    Orientation<-c()
    Br.Angle<-c()
    MeanAngleVar<-c()
    SDAngleVar<-c()
    rac<-list()
    tip<-list()
  
  # Storing the final root length of each root constituting the vectorized root systems
  
  for (i in 1:length(LIE)){
    
    finallength<-c()
    k<-0
    
    for (j in 1:nrow(LIE[[i]])){
      
      if (LIE[[i]]$Bran[j]=="true"){
        k<-k+1
        prec<-LIE[[i]]$Prec[j]
        
        if (prec==0){
          finallength[k]<-0
          m<-LIE[[i]]$Suiv[j]}
        
        if (prec!=0){
          prec<-LIE[[i]]$Prec[j]
          if (!("Z" %in% colnames(LIE[[i]]))) {finallength[k]<-sqrt((LIE[[i]]$X[j]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[j]-LIE[[i]]$Y[prec])^2)} else {finallength[k]<-sqrt((LIE[[i]]$X[j]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[j]-LIE[[i]]$Y[prec])^2+(LIE[[i]]$Z[j]-LIE[[i]]$Z[prec])^2)}
          m<-LIE[[i]]$Suiv[j]}
        
          while (m!=0){
            prec<-LIE[[i]]$Prec[m]
            if (!("Z" %in% colnames(LIE[[i]]))) {finallength[k]<-finallength[k]+sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2)} else {finallength[k]<-finallength[k]+sqrt((LIE[[i]]$X[m]-LIE[[i]]$X[prec])^2+(LIE[[i]]$Y[m]-LIE[[i]]$Y[prec])^2+(LIE[[i]]$Z[m]-LIE[[i]]$Z[prec])^2)}
            m<-LIE[[i]]$Suiv[m]}}}
    
    finallength.lie[[i]]<-finallength}
    
  # Calculating the coordinates of interpolated points
  
  t<-1
  for (i in 1:length(LIE)){
    
    LIE[[i]]$X<-LIE[[i]]$X-min(LIE[[i]]$X)
    LIE[[i]]$Y<-LIE[[i]]$Y-min(LIE[[i]]$Y)
    if ("Z" %in% colnames(LIE[[i]])) {LIE[[i]]$Z<-LIE[[i]]$Z-min(LIE[[i]]$Z)}
    
    k<-0
    XYcurv<-list()
    XYangle<-list()
    XYdir<-list()
    if (!("Z" %in% colnames(LIE[[i]]))) {orientation<-c()}
    tortuosity<-c()
    if (length(TIME)==1){
      num<-TIME[[1]]$Num
      date<-TIME[[1]]$Date}
    else {
      num<-TIME[[i]]$Num
      date<-TIME[[i]]$Date}
    tipangle<-matrix(nrow=nrow(DATA[[i]]), ncol=length(date))
    
    for (j in 1:nrow(LIE[[i]])){
      
      if (LIE[[i]]$Bran[j]=="true"){
        
        k<-k+1
        l<-1
        distangle<-0
        prec1<-LIE[[i]]$Prec[j]
        if (finallength.lie[[i]][k]<l.tipangle) {tipangle[k,]<-NA}
        
        # Points used for root curvature calculation
        
        if (finallength.lie[[i]][k]<l.curv){XYcurv[[k]]<-NA}
        else {
          if (!("Z" %in% colnames(LIE[[i]]))) {XYcurv[[k]]<-matrix(ncol=2, nrow=floor(finallength.lie[[i]][k]/l.curv)+1)} else {XYcurv[[k]]<-matrix(ncol=3, nrow=floor(finallength.lie[[i]][k]/l.curv)+1)}
        if (prec1==0){
          if (!("Z" %in% colnames(LIE[[i]]))){
              XYcurv[[k]][l,1]<-LIE[[i]]$X[j]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[j]}
          else {
              XYcurv[[k]][l,1]<-LIE[[i]]$X[j]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[j]
              XYcurv[[k]][l,3]<-LIE[[i]]$Z[j]}
          
              suiv<-LIE[[i]]$Suiv[j]
              l<-l+1}
        
        if (prec1!=0){
          if (!("Z" %in% colnames(LIE[[i]]))){
              XYcurv[[k]][l,1]<-LIE[[i]]$X[LIE[[i]]$Prec[j]]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[LIE[[i]]$Prec[j]]}
          else {
              XYcurv[[k]][l,1]<-LIE[[i]]$X[LIE[[i]]$Prec[j]]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[LIE[[i]]$Prec[j]]
              XYcurv[[k]][l,3]<-LIE[[i]]$Z[LIE[[i]]$Prec[j]]}
              
              suiv<-j
              l<-l+1}
                
        if (!("Z" %in% colnames(LIE[[i]]))) {D<-distance2D(x1=XYcurv[[k]][l-1,1], y1=XYcurv[[k]][l-1,2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])} else {D<-distance3D(x1=XYcurv[[k]][l-1,1], y1=XYcurv[[k]][l-1,2], z1=XYcurv[[k]][l-1,3], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv])}
        distcurv<-D
        
        while (l<=floor(finallength.lie[[i]][k]/l.curv)+1){
          
          # First situation
          if (distcurv<l.curv){
            while (distcurv<l.curv){
              suiv<-LIE[[i]]$Suiv[suiv]
              if (!("Z" %in% colnames(LIE[[i]]))) {D<-distance2D(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])} else {D<-distance3D(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], z1=LIE[[i]]$Z[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv])}
              distcurv<-distcurv+D}
            
            if (distcurv==l.curv){
              if (!("Z" %in% colnames(LIE[[i]]))){
                XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]}
              else {
                XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
                XYcurv[[k]][l,3]<-LIE[[i]]$Z[suiv]}
              l<-l+1
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
              distcurv<-0}
            
            if (distcurv>l.curv){
              if (!("Z" %in% colnames(LIE[[i]]))){
                xycoord<-XYcoord2D(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], d=l.curv-distcurv+D)
                XYcurv[[k]][l,1]<-xycoord[1]
                XYcurv[[k]][l,2]<-xycoord[2]}
              else {
                xycoord<-XYcoord3D(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], z1=LIE[[i]]$Z[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv], d=l.curv-distcurv+D)
                XYcurv[[k]][l,1]<-xycoord[1]
                XYcurv[[k]][l,2]<-xycoord[2]
                XYcurv[[k]][l,3]<-xycoord[3]}
              l<-l+1
              
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
              if (!("Z" %in% colnames(LIE[[i]]))) {distcurv<-distance2D(x1=xycoord[1], y1=xycoord[2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])} else {distcurv<-distance3D(x1=xycoord[1], y1=xycoord[2], z1=xycoord[3], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv])}
              
              if (distcurv==l.curv){
                if (!("Z" %in% colnames(LIE[[i]]))){
                  XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                  XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]}
                else {
                  XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                  XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
                  XYcurv[[k]][l,3]<-LIE[[i]]$Z[suiv]}
                
                l<-l+1
                if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
                distcurv<-0}
              
              while (distcurv>l.curv){
                if (!("Z" %in% colnames(LIE[[i]]))){
                  xycoord<-XYcoord2D(x1=XYcurv[[k]][l-1,1], y1=XYcurv[[k]][l-1,2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], d=l.curv)
                  XYcurv[[k]][l,1]<-xycoord[1]
                  XYcurv[[k]][l,2]<-xycoord[2]}
                else {
                  xycoord<-XYcoord3D(x1=XYcurv[[k]][l-1,1], y1=XYcurv[[k]][l-1,2], z1=XYcurv[[k]][l-1,3], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv], d=l.curv)
                  XYcurv[[k]][l,1]<-xycoord[1]
                  XYcurv[[k]][l,2]<-xycoord[2]
                  XYcurv[[k]][l,3]<-xycoord[3]}
                
                l<-l+1
                if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
                
                if (!("Z" %in% colnames(LIE[[i]]))) {distcurv<-distance2D(x1=xycoord[1], y1=xycoord[2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])} else {distcurv<-distance3D(x1=xycoord[1], y1=xycoord[2], z1=xycoord[3], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv])}
                
                if (distcurv==l.curv){
                  if (!("Z" %in% colnames(LIE[[i]]))){
                    XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                    XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]}
                  else {
                    XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                    XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
                    XYcurv[[k]][l,3]<-LIE[[i]]$Z[suiv]}
                  
                  l<-l+1
                  if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
                  distcurv<-0}}
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}}
            suiv<-LIE[[i]]$Suiv[suiv]}
          
          # Second situation
          if (distcurv>l.curv){
            if (!("Z" %in% colnames(LIE[[i]]))){
              xycoord<-XYcoord2D(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], d=l.curv-distcurv+D)
              XYcurv[[k]][l,1]<-xycoord[1]
              XYcurv[[k]][l,2]<-xycoord[2]}
            else {
              xycoord<-XYcoord3D(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], z1=LIE[[i]]$Z[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv], d=l.curv-distcurv+D)
              XYcurv[[k]][l,1]<-xycoord[1]
              XYcurv[[k]][l,2]<-xycoord[2]
              XYcurv[[k]][l,3]<-xycoord[3]}
            
            l<-l+1
            if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
            if (!("Z" %in% colnames(LIE[[i]]))) {distcurv<-distance2D(x1=xycoord[1], y1=xycoord[2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])} else {distcurv<-distance3D(x1=xycoord[1], y1=xycoord[2], z1=xycoord[3], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv])}
            
            if (distcurv==l.curv){
              if (!("Z" %in% colnames(LIE[[i]]))){
                XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]}
              else {
                XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
                XYcurv[[k]][l,3]<-LIE[[i]]$Z[suiv]}
              
              l<-l+1
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
              distcurv<-0}
            
            while (distcurv>l.curv){
              if (!("Z" %in% colnames(LIE[[i]]))){
                xycoord<-XYcoord2D(x1=XYcurv[[k]][l-1,1], y1=XYcurv[[k]][l-1,2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], d=l.curv)
                XYcurv[[k]][l,1]<-xycoord[1]
                XYcurv[[k]][l,2]<-xycoord[2]}
              else {
                xycoord<-XYcoord3D(x1=XYcurv[[k]][l-1,1], y1=XYcurv[[k]][l-1,2], z1=XYcurv[[k]][l-1,3], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv], d=l.curv)
                XYcurv[[k]][l,1]<-xycoord[1]
                XYcurv[[k]][l,2]<-xycoord[2]
                XYcurv[[k]][l,3]<-xycoord[3]}
              
              l<-l+1
              if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
              
              if (!("Z" %in% colnames(LIE[[i]]))) {distcurv<-distance2D(x1=xycoord[1], y1=xycoord[2], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])} else {distcurv<-distance3D(x1=xycoord[1], y1=xycoord[2], z1=xycoord[3], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv])}
              
              if (distcurv==l.curv){
                if (!("Z" %in% colnames(LIE[[i]]))){
                  XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                  XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]}
                else {
                  XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
                  XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
                  XYcurv[[k]][l,3]<-LIE[[i]]$Z[suiv]}
                
                l<-l+1
                if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
                distcurv<-0}}
            
            if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
            
          suiv<-LIE[[i]]$Suiv[suiv]}
          
          # Third situation
          if (distcurv==l.curv){
            if (!("Z" %in% colnames(LIE[[i]]))){
              XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]}
            else {
              XYcurv[[k]][l,1]<-LIE[[i]]$X[suiv]
              XYcurv[[k]][l,2]<-LIE[[i]]$Y[suiv]
              XYcurv[[k]][l,3]<-LIE[[i]]$Z[suiv]}
            
            l<-l+1
            if (l>floor(finallength.lie[[i]][k]/l.curv)+1) {break}
            distcurv<-0
            suiv<-LIE[[i]]$Suiv[suiv]}
        
          if (!("Z" %in% colnames(LIE[[i]]))) {D<-distance2D(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv])} else {D<-distance3D(x1=LIE[[i]]$X[LIE[[i]]$Prec[suiv]], y1=LIE[[i]]$Y[LIE[[i]]$Prec[suiv]], z1=LIE[[i]]$Z[LIE[[i]]$Prec[suiv]], x2=LIE[[i]]$X[suiv], y2=LIE[[i]]$Y[suiv], z2=LIE[[i]]$Z[suiv])}
          distcurv<-distcurv+D}}
      
      # Points used to calculate the branching angles of daughter roots on mother roots
        
        if (prec1==0) {if (!("Z" %in% colnames(LIE[[i]]))) {orientation[k]<-NA}}
        
        if (finallength.lie[[i]][k]<l.brangle){XYangle[[k]]<-NA}
        else{
          if (!("Z" %in% colnames(LIE[[i]]))){
            if (prec1!=0){
              XYangle[[k]]<-matrix(nrow=3, ncol=2)
              XYangle[[k]][1,1]<-LIE[[i]]$X[LIE[[i]]$Prec[j]]
              XYangle[[k]][1,2]<-LIE[[i]]$Y[LIE[[i]]$Prec[j]]}
            else {
              XYangle[[k]]<-matrix(nrow=3, ncol=2)
              XYangle[[k]][1,1]<-LIE[[i]]$X[j]
              XYangle[[k]][1,2]<-LIE[[i]]$Y[j]}}
          
          else {
            if (prec1!=0){
              XYangle[[k]]<-matrix(nrow=3, ncol=3)
              XYangle[[k]][1,1]<-LIE[[i]]$X[LIE[[i]]$Prec[j]]
              XYangle[[k]][1,2]<-LIE[[i]]$Y[LIE[[i]]$Prec[j]]
              XYangle[[k]][1,3]<-LIE[[i]]$Z[LIE[[i]]$Prec[j]]}
            else{
              XYangle[[k]]<-matrix(nrow=3, ncol=3)
              XYangle[[k]][1,1]<-LIE[[i]]$X[j]
              XYangle[[k]][1,2]<-LIE[[i]]$Y[j]
              XYangle[[k]][1,3]<-LIE[[i]]$Z[j]}}
            
            # For daughter roots
            m<-j
            
            while(distangle<l.brangle){
              prec<-LIE[[i]]$Prec[m]
              
              if (!("Z" %in% colnames(LIE[[i]]))) {
                
                if (prec!=0) {D<-distance2D(x1=LIE[[i]]$X[prec], x2=LIE[[i]]$X[m], y1=LIE[[i]]$Y[prec], y2=LIE[[i]]$Y[m])}
                else {D<-0}} 
              
              else {
                
                if (prec!=0) {D<-distance3D(x1=LIE[[i]]$X[prec], x2=LIE[[i]]$X[m], y1=LIE[[i]]$Y[prec], y2=LIE[[i]]$Y[m], z1=LIE[[i]]$Z[prec], z2=LIE[[i]]$Z[m])}
                else {D<-0}}
              
              distangle<-distangle+D
              if (distangle<l.brangle) {m<-LIE[[i]]$Suiv[m]}}
            
            if (distangle==l.brangle){
              if (!("Z" %in% colnames(LIE[[i]]))){
                XYangle[[k]][3,1]<-LIE[[i]]$X[m]
                XYangle[[k]][3,2]<-LIE[[i]]$Y[m]}
              else {
                XYangle[[k]][3,1]<-LIE[[i]]$X[m]
                XYangle[[k]][3,2]<-LIE[[i]]$Y[m]
                XYangle[[k]][3,3]<-LIE[[i]]$Z[m]}}
            
            if (distangle>l.brangle){
              if (!("Z" %in% colnames(LIE[[i]]))) {
                xycoord<-XYcoord2D(x1=LIE[[i]]$X[prec], x2=LIE[[i]]$X[m], y1=LIE[[i]]$Y[prec], y2=LIE[[i]]$Y[m], d=l.brangle-distangle+D)
                XYangle[[k]][3,1]<-xycoord[1]
                XYangle[[k]][3,2]<-xycoord[2]}
              else {
                xycoord<-XYcoord3D(x1=LIE[[i]]$X[prec], x2=LIE[[i]]$X[m], y1=LIE[[i]]$Y[prec], y2=LIE[[i]]$Y[m], z1=LIE[[i]]$Z[prec], z2=LIE[[i]]$Z[m], d=l.brangle-distangle+D)
                XYangle[[k]][3,1]<-xycoord[1]
                XYangle[[k]][3,2]<-xycoord[2]
                XYangle[[k]][3,3]<-xycoord[3]}}
            
            # For mother roots
            
            if (prec1==0){
              
              if (!("Z" %in% colnames(LIE[[i]]))){
                XYangle[[k]][2,1]<-XYangle[[k]][1,1]+0
                XYangle[[k]][2,2]<-XYangle[[k]][1,2]+1}
              else {
                XYangle[[k]][2,1]<-XYangle[[k]][1,1]+dirvert[1]
                XYangle[[k]][2,2]<-XYangle[[k]][1,2]+dirvert[2]
                XYangle[[k]][2,3]<-XYangle[[k]][1,3]+dirvert[3]}}
            
            else{
              
              if((finallength.lie[[i]][DATA[[i]]$Mother[k]+1]-(cunit*DATA[[i]]$DBase[k]))<l.brangle) {XYangle[[k]]<-NA}
              
              else{
            
            distangle<-0
            m<-LIE[[i]]$Prec[j]
            
            while(distangle<l.brangle){
              suiv<-LIE[[i]]$Suiv[m]
              if (!("Z" %in% colnames(LIE[[i]]))) {D<-distance2D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[suiv], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[suiv])} else {D<-distance3D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[suiv], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[suiv], z1=LIE[[i]]$Z[m], z2=LIE[[i]]$Z[suiv])}
              distangle<-distangle+D
              if (distangle<l.brangle) {m<-suiv}}
            
            if (distangle==l.brangle){
              if (!("Z" %in% colnames(LIE[[i]]))){
                XYangle[[k]][2,1]<-LIE[[i]]$X[suiv]
                XYangle[[k]][2,2]<-LIE[[i]]$Y[suiv]}
              else {
                XYangle[[k]][2,1]<-LIE[[i]]$X[suiv]
                XYangle[[k]][2,2]<-LIE[[i]]$Y[suiv]
                XYangle[[k]][2,3]<-LIE[[i]]$Z[suiv]}}
            
            if (distangle>l.brangle){
              if (!("Z" %in% colnames(LIE[[i]]))){
                xycoord<-XYcoord2D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[suiv], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[suiv], d=l.brangle-distangle+D)
                XYangle[[k]][2,1]<-xycoord[1]
                XYangle[[k]][2,2]<-xycoord[2]}
              else{
                xycoord<-XYcoord3D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[suiv], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[suiv], z1=LIE[[i]]$Z[m], z2=LIE[[i]]$Z[suiv], d=l.brangle-distangle+D)
                XYangle[[k]][2,1]<-xycoord[1]
                XYangle[[k]][2,2]<-xycoord[2]
                XYangle[[k]][2,3]<-xycoord[3]}}}}}
          
          # Calculating the Orientation of lateral roots
          
          if (prec1!=0 & !("Z" %in% colnames(LIE[[i]]))) {
            prec<-LIE[[i]]$Prec[j]
            suivMR<-LIE[[i]]$Suiv[prec]
            suivLR<-LIE[[i]]$Suiv[j]
            u<-c(LIE[[i]]$X[suivMR]-LIE[[i]]$X[prec], LIE[[i]]$Y[suivMR]-LIE[[i]]$Y[prec])
            while (u[1]==0 & u[2]==0){
              suivMR<-LIE[[i]]$Suiv[suivMR]
              u<-c(LIE[[i]]$X[suivMR]-LIE[[i]]$X[prec], LIE[[i]]$Y[suivMR]-LIE[[i]]$Y[prec])}
            n<-normal(u)
            lateral<-c(LIE[[i]]$X[j]-LIE[[i]]$X[prec], LIE[[i]]$Y[j]-LIE[[i]]$Y[prec])
            while (lateral[1]==0 & lateral[2]==0){
              lateral<-c(LIE[[i]]$X[suivLR]-LIE[[i]]$X[prec], LIE[[i]]$Y[suivLR]-LIE[[i]]$Y[prec])
              suivLR<-LIE[[i]]$Suiv[suivLR]}
          
            if (n%*%lateral>0){orientation[k]<-"Left"}
            if (n%*%lateral<0){orientation[k]<-"Right"}
          
            if (n%*%lateral==0 & suivLR==0){
              if (runif(1, min=-1, max=1)>0) {orientation[k]<-"Left"} else {orientation[k]<-"Right"}}
          
            if (n%*%lateral==0 & suivLR!=0){
            
            while(n%*%lateral==0){
                lateral<-c(LIE[[i]]$X[suivLR]-LIE[[i]]$X[prec], LIE[[i]]$Y[suivLR]-LIE[[i]]$Y[prec])
                suivLR<-LIE[[i]]$Suiv[suivLR]
                if (suivLR==0 & n%*%lateral==0){
                  if (runif(1, min=-1, max=1)>0) {orientation[k]<-"Left"} else {orientation[k]<-"Right"}
                  break}}
            
              if (n%*%lateral>0){orientation[k]<-"Left"}
              if (n%*%lateral<0){orientation[k]<-"Right"}}}
      
      # Calculating the points used for the calculation of tip angles
        
      if (finallength.lie[[i]][k]>=l.tipangle){
       
        m<-j
        
        while (LIE[[i]]$Apic[m]!="true"){
          
          a<-LIE[[i]]$Date[m]
          b<-LIE[[i]]$Date[LIE[[i]]$Suiv[m]]
          
          if (a!=0 & a<b){
            
            if (DATA[[i]][k,5+a]*cunit>=l.tipangle) {
            
            if (!("Z" %in% colnames(LIE[[i]]))) {disttip<-distance2D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]])} else {disttip<-distance3D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]], z1=LIE[[i]]$Z[m], z2=LIE[[i]]$Z[LIE[[i]]$Prec[m]])}
            if (disttip>l.tipangle){if (!("Z" %in% colnames(LIE[[i]]))) {xycoord<-XYcoord2D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]], d=l.tipangle)} else {xycoord<-XYcoord3D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]], z1=LIE[[i]]$Z[m], z2=LIE[[i]]$Z[LIE[[i]]$Prec[m]], d=l.tipangle)}}
            if (disttip==l.tipangle){if (!("Z" %in% colnames(LIE[[i]]))) {xycoord<-c(LIE[[i]]$X[LIE[[i]]$Prec[m]], LIE[[i]]$Y[LIE[[i]]$Prec[m]])} else {xycoord<-c(LIE[[i]]$X[LIE[[i]]$Prec[m]], LIE[[i]]$Y[LIE[[i]]$Prec[m]], LIE[[i]]$Z[LIE[[i]]$Prec[m]])}}
            if (disttip<l.tipangle){
              prec1<-LIE[[i]]$Prec[m]
              while(disttip<l.tipangle){
                prec2<-LIE[[i]]$Prec[prec1]
                if (!("Z" %in% colnames(LIE[[i]]))) {D<-distance2D(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2])} else {D<-distance3D(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2], z1=LIE[[i]]$Z[prec1], z2=LIE[[i]]$Z[prec2])}
                disttip<-disttip+D
                if (disttip<l.tipangle) {prec1<-prec2}}
              if (disttip==l.tipangle){if (!("Z" %in% colnames(LIE[[i]]))) {xycoord<-c(LIE[[i]]$X[prec2], LIE[[i]]$Y[prec2])} else {xycoord<-c(LIE[[i]]$X[prec2], LIE[[i]]$Y[prec2], LIE[[i]]$Z[prec2])}}
              if (disttip>l.tipangle) {if (!("Z" %in% colnames(LIE[[i]]))) {xycoord<-XYcoord2D(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2], d=l.tipangle-disttip+D)} else {xycoord<-XYcoord3D(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2], z1=LIE[[i]]$Z[prec1], z2=LIE[[i]]$Z[prec2], d=l.tipangle-disttip+D)}}}
            
              if (!("Z" %in% colnames(LIE[[i]]))) {tipangle[k,a:(b-1)]<-acos((c(0,1)%*%c(LIE[[i]]$X[m]-xycoord[1], LIE[[i]]$Y[m]-xycoord[2]))/(sqrt((LIE[[i]]$X[m]-xycoord[1])^2+(LIE[[i]]$Y[m]-xycoord[2])^2)))*cunitangle} else {tipangle[k,a:(b-1)]<-acos((dirvert%*%c(LIE[[i]]$X[m]-xycoord[1], LIE[[i]]$Y[m]-xycoord[2], LIE[[i]]$Z[m]-xycoord[3]))/(sqrt((LIE[[i]]$X[m]-xycoord[1])^2+(LIE[[i]]$Y[m]-xycoord[2])^2+(LIE[[i]]$Z[m]-xycoord[3])^2)))*cunitangle}}}
          
          m<-LIE[[i]]$Suiv[m]}
      
      if (LIE[[i]]$Apic[m]=="true"){
                
          a<-LIE[[i]]$Date[m]
        
          if (!("Z" %in% colnames(LIE[[i]]))) {disttip<-distance2D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]])} else {disttip<-distance3D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]], z1=LIE[[i]]$Z[m], z2=LIE[[i]]$Z[LIE[[i]]$Prec[m]])}
          if (disttip>l.tipangle){if (!("Z" %in% colnames(LIE[[i]]))) {xycoord<-XYcoord2D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]], d=l.tipangle)} else {xycoord<-XYcoord3D(x1=LIE[[i]]$X[m], x2=LIE[[i]]$X[LIE[[i]]$Prec[m]], y1=LIE[[i]]$Y[m], y2=LIE[[i]]$Y[LIE[[i]]$Prec[m]], z1=LIE[[i]]$Z[m], z2=LIE[[i]]$Z[LIE[[i]]$Prec[m]], d=l.tipangle)}}
          if (disttip==l.tipangle){if (!("Z" %in% colnames(LIE[[i]]))) {xycoord<-c(LIE[[i]]$X[LIE[[i]]$Prec[m]], LIE[[i]]$Y[LIE[[i]]$Prec[m]])} else {xycoord<-c(LIE[[i]]$X[LIE[[i]]$Prec[m]], LIE[[i]]$Y[LIE[[i]]$Prec[m]], LIE[[i]]$Z[LIE[[i]]$Prec[m]])}}
          if (disttip<l.tipangle){
            prec1<-LIE[[i]]$Prec[m]
            while(disttip<l.tipangle){
              prec2<-LIE[[i]]$Prec[prec1]
              if (!("Z" %in% colnames(LIE[[i]]))) {D<-distance2D(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2])} else {D<-distance3D(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2], z1=LIE[[i]]$Z[prec1], z2=LIE[[i]]$Z[prec2])}
              disttip<-disttip+D
              if (disttip<l.tipangle) {prec1<-prec2}}
            if (disttip==l.tipangle){if (!("Z" %in% colnames(LIE[[i]]))) {xycoord<-c(LIE[[i]]$X[prec2], LIE[[i]]$Y[prec2])} else {xycoord<-c(LIE[[i]]$X[prec2], LIE[[i]]$Y[prec2], LIE[[i]]$Z[prec2])}}
            if (disttip>l.tipangle) {if (!("Z" %in% colnames(LIE[[i]]))) {xycoord<-XYcoord2D(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2], d=l.tipangle-disttip+D)} else {xycoord<-XYcoord3D(x1=LIE[[i]]$X[prec1], x2=LIE[[i]]$X[prec2], y1=LIE[[i]]$Y[prec1], y2=LIE[[i]]$Y[prec2], z1=LIE[[i]]$Z[prec1], z2=LIE[[i]]$Z[prec2], d=l.tipangle-disttip+D)}}}
          
          if (!("Z" %in% colnames(LIE[[i]]))) {tipangle[k,(a:ncol(tipangle))]<-acos((c(0,1)%*%c(LIE[[i]]$X[m]-xycoord[1], LIE[[i]]$Y[m]-xycoord[2]))/(sqrt((LIE[[i]]$X[m]-xycoord[1])^2+(LIE[[i]]$Y[m]-xycoord[2])^2)))*cunitangle} else {tipangle[k,(a:ncol(tipangle))]<-acos((dirvert%*%c(LIE[[i]]$X[m]-xycoord[1], LIE[[i]]$Y[m]-xycoord[2], LIE[[i]]$Z[m]-xycoord[3]))/(sqrt((LIE[[i]]$X[m]-xycoord[1])^2+(LIE[[i]]$Y[m]-xycoord[2])^2+(LIE[[i]]$Z[m]-xycoord[3])^2)))*cunitangle}}}
      
      # Tortuosity
      
      m<-j
      prec<-LIE[[i]]$Prec[j]
      
      if (!("Z" %in% colnames(LIE[[i]]))) {XYdir[[k]]<-matrix(nrow=2, ncol=2)} else {XYdir[[k]]<-matrix(nrow=2, ncol=3)}
      
      while (LIE[[i]]$Apic[m]!="true"){m<-LIE[[i]]$Suiv[m]}
      
      if (prec==0) {
        if (!("Z" %in% colnames(LIE[[i]]))) {
          tortuosity[k]<-finallength.lie[[i]][k]/sqrt((LIE[[i]]$X[j]-LIE[[i]]$X[m])^2+(LIE[[i]]$Y[j]-LIE[[i]]$Y[m])^2)
          XYdir[[k]][1,]<-c(LIE[[i]]$X[j], LIE[[i]]$Y[j])
          XYdir[[k]][2,]<-c(LIE[[i]]$X[m], LIE[[i]]$Y[m])}
        
        else {
          tortuosity[k]<-finallength.lie[[i]][k]/sqrt((LIE[[i]]$X[j]-LIE[[i]]$X[m])^2+(LIE[[i]]$Y[j]-LIE[[i]]$Y[m])^2+(LIE[[i]]$Z[j]-LIE[[i]]$Z[m])^2)
          XYdir[[k]][1,]<-c(LIE[[i]]$X[j], LIE[[i]]$Y[j], LIE[[i]]$Z[j])
          XYdir[[k]][2,]<-c(LIE[[i]]$X[m], LIE[[i]]$Y[m], LIE[[i]]$Z[j])}} 
      
      else {
        if (!("Z" %in% colnames(LIE[[i]]))) {
          tortuosity[k]<-finallength.lie[[i]][k]/sqrt((LIE[[i]]$X[prec]-LIE[[i]]$X[m])^2+(LIE[[i]]$Y[prec]-LIE[[i]]$Y[m])^2)
          XYdir[[k]][1,]<-c(LIE[[i]]$X[prec], LIE[[i]]$Y[prec])
          XYdir[[k]][2,]<-c(LIE[[i]]$X[m], LIE[[i]]$Y[m])}
        
        else {
          tortuosity[k]<-finallength.lie[[i]][k]/sqrt((LIE[[i]]$X[prec]-LIE[[i]]$X[m])^2+(LIE[[i]]$Y[prec]-LIE[[i]]$Y[m])^2+(LIE[[i]]$Z[prec]-LIE[[i]]$Z[m])^2)
          XYdir[[k]][1,]<-c(LIE[[i]]$X[prec], LIE[[i]]$Y[prec], LIE[[i]]$Z[prec])
          XYdir[[k]][2,]<-c(LIE[[i]]$X[m], LIE[[i]]$Y[m], LIE[[i]]$Z[j])}}}}
    
    tip[[i]]<-data.frame(DATA[[i]]$Root, tipangle)
    colnames(tip[[i]])<-c("Root", paste("Ang.Date", t(num), sep=""))
    
  # Calculating branching angles and curvatures + root growth direction relative to the vertical
  
    br.angle<-c()
    gr.dir<-c()
    meananglevar<-c()
    sdanglevar<-c()
    
    for (j in 1:length(XYangle)){
      
      if ("matrix" %in% class(XYdir[[j]])){
        
        VECTangle<-XYdir[[j]][2,]-XYdir[[j]][1,]
        normVECTangle<-sqrt(VECTangle[1]^2+VECTangle[2]^2)
        if (length(VECTangle)==2) {gr.dir[j]<-acos((VECTangle%*%c(0,1))/(normVECTangle*1))*cunitangle}
        if (length(VECTangle)==3) {gr.dir[j]<-acos((VECTangle%*%dirvert)/(normVECTangle*1))*cunitangle}}
      
      if ("matrix" %in% class(XYangle[[j]])){
        
        if (!("Z" %in% colnames(LIE[[i]]))){
          VECTangle<-matrix(nrow=2, ncol=2)
          VECTangle[1,]<-XYangle[[j]][2,]-XYangle[[j]][1,] # For mother roots
          VECTangle[2,]<-XYangle[[j]][3,]-XYangle[[j]][1,] # For daugther roots
          normVECTangle<-sqrt(VECTangle[,1]^2+VECTangle[,2]^2)
          br.angle[j]<-acos((VECTangle[1,]%*%VECTangle[2,])/(normVECTangle[1]*normVECTangle[2]))*cunitangle}
        else {
          VECTangle<-matrix(nrow=2, ncol=3)
          VECTangle[1,]<-XYangle[[j]][2,]-XYangle[[j]][1,] # For mother roots
          VECTangle[2,]<-XYangle[[j]][3,]-XYangle[[j]][1,] # For daugther roots
          normVECTangle<-sqrt(VECTangle[,1]^2+VECTangle[,2]^2+VECTangle[,3]^2)
          br.angle[j]<-acos((VECTangle[1,]%*%VECTangle[2,])/(normVECTangle[1]*normVECTangle[2]))*cunitangle}} 
      
      else {br.angle[j]<-NA}
      
      if ("matrix" %in% class(XYcurv[[j]])) {
        
        if (nrow(XYcurv[[j]])>2){
          
          VECTcurv<-diff(XYcurv[[j]])
          angle<-c()
          for (k in 1:(nrow(XYcurv[[j]])-2)){
            if (!("Z" %in% colnames(LIE[[i]]))) {ratio<-(VECTcurv[k,]%*%VECTcurv[k+1,])/(sqrt(VECTcurv[k,1]^2+VECTcurv[k,2]^2)*sqrt(VECTcurv[k+1,1]^2+VECTcurv[k+1,2]^2))} else {ratio<-(VECTcurv[k,]%*%VECTcurv[k+1,])/(sqrt(VECTcurv[k,1]^2+VECTcurv[k,2]^2+VECTcurv[k,3]^2)*sqrt(VECTcurv[k+1,1]^2+VECTcurv[k+1,2]^2+VECTcurv[k+1,3]^2))}
            if (ratio>1) {ratio<-1} # The value of a cosinus must be < or equal to 1
            if (ratio<(1*(-1))) {ratio<-1*(-1)} # The value of a cosinus must be > or equal to -1
            angle[k]<-acos(ratio)*cunitangle}
          meananglevar[j]<-mean(angle)
          sdanglevar[j]<-sd(angle)}
      
        else {
          
          meananglevar[j]<-NA
          sdanglevar[j]<-NA}} 
      
      else {
        
        meananglevar[j]<-NA
        sdanglevar[j]<-NA}}
  
  if (!("Z" %in% colnames(LIE[[i]]))) {rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Mother=DATA[[i]]$Mother, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase*cunit, FinalRootLength=finallength.lie[[i]], Tortuosity=tortuosity, Orientation=orientation, Branching.Angle=br.angle, Growth.Direction=gr.dir, Mean.Curv=meananglevar, SD.Curv=sdanglevar)}
  else {rac[[i]]<-data.frame(Root=DATA[[i]]$Root, Mother=DATA[[i]]$Mother, Ord=DATA[[i]]$Ord, DBase=DATA[[i]]$DBase*cunit, FinalRootLength=finallength.lie[[i]], Tortuosity=tortuosity, Branching.Angle=br.angle, Growth.Direction=gr.dir, Mean.Curv=meananglevar, SD.Curv=sdanglevar)}
  
  #Plotting root system with a color code for the branching angle
  
  if (is.null(plot)==FALSE){
    
  if (plot=="branching"){
  
    if (is.null(BRscale)==TRUE) {
      maxi<-max(rac[[i]]$Branching.Angle[2:length(rac[[i]]$Branching.Angle)], na.rm=TRUE)
      mini<-min(rac[[i]]$Branching.Angle[2:length(rac[[i]]$Branching.Angle)], na.rm=TRUE)}
    else {
      maxi<-max(BRscale)
      mini<-min(BRscale)}
    branchingangled<-(rac[[i]]$Branching.Angle[2:length(rac[[i]]$Branching.Angle)]-mini)/(maxi-mini)
    posNEG<-which(branchingangled<0)
    pos1plus<-which(branchingangled>1)
    posNA<-which(is.na(branchingangled))
    branchingangled[is.na(branchingangled)|branchingangled<0|branchingangled>1]<-0
    pal<-colorRamp(colangle)
    colors<-rgb(pal(branchingangled), maxColorValue=255)
    colors[posNA]<-"black"
    colors[posNEG]<-"black"
    colors[pos1plus]<-"black"
    colors<-c("black", colors)
    if (export.colors==TRUE){rac[[i]]$Colors<-colors}}
    
  if (plot=="direction"){
    
    if (is.null(BRscale)==TRUE) {
      maxi<-max(rac[[i]]$Growth.Direction[2:length(rac[[i]]$Growth.Direction)], na.rm=TRUE)
      mini<-min(rac[[i]]$Growth.Direction[2:length(rac[[i]]$Growth.Direction)], na.rm=TRUE)}
    else {
      maxi<-max(BRscale)
      mini<-min(BRscale)}
    growthdirectiond<-(rac[[i]]$Growth.Direction[2:length(rac[[i]]$Growth.Direction)]-mini)/(maxi-mini)
    posNEG<-which(growthdirectiond<0)
    pos1plus<-which(growthdirectiond>1)
    posNA<-which(is.na(growthdirectiond))
    growthdirectiond[is.na(growthdirectiond)|growthdirectiond<0|growthdirectiond>1]<-0
    pal<-colorRamp(colangle)
    colors<-rgb(pal(growthdirectiond), maxColorValue=255)
    colors[posNA]<-"black"
    colors[posNEG]<-"black"
    colors[pos1plus]<-"black"
    colors<-c("black", colors)
    if (export.colors==TRUE){rac[[i]]$Colors<-colors}}
  
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
  
  if (!("Z" %in% colnames(LIE[[i]]))){
    
    plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
    r<-0
    for (k in 1:nrow(LIE[[i]])){
      if(LIE[[i]]$Bran[k]=="true"){
        r<-r+1
        a<-LIE[[i]]$Prec[k]
        b<-LIE[[i]]$Date[k]
        if (a!=0) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=colors[r],...)}
        m<-LIE[[i]]$Suiv[k]
        if (m!=0){
          while (LIE[[i]]$Apic[m]=="false"){
            a<-LIE[[i]]$Prec[m]
            b<-LIE[[i]]$Date[m]
            segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r],...)
            m<-LIE[[i]]$Suiv[m]}
          if (LIE[[i]]$Apic[m]=="true"){
            a<-LIE[[i]]$Prec[m]
            b<-LIE[[i]]$Date[m]
            segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r],...)}}}}}
  
  else {
    
    if (is.null(twod)==TRUE){
    
    root<-sum(LIE[[i]]$Suiv==0)
    end<-which(LIE[[i]]$Suiv==0)
    open3d()
    plot3d(x=LIE[[i]]$X[1], y=LIE[[i]]$Y[1], z=LIE[[i]]$Z[1], type="n", xlim=xlim1, ylim=ylim1, zlim=zlim1, main=main1, ylab=ylab1, xlab=xlab1, zlab=zlab1,...)
    for (k in 1:root){
      if (k==1) {
        dataroot<-as.matrix(LIE[[i]][k:end[k],7:9])} 
      else {
        if (LIE[[i]]$Prec[end[k-1]+1]!=0){
          dataroot<-as.matrix(LIE[[i]][(end[k-1]+1):end[k],7:9])
          dataroot<-rbind(LIE[[i]][LIE[[i]]$Num==LIE[[i]]$Prec[end[k-1]+1],7:9], dataroot)}
        else{
          dataroot<-as.matrix(LIE[[i]][(end[k-1]+1):end[k],7:9])}}
  
      lines3d(dataroot, col=colors[k], smooth=FALSE, ...)}}
    
    else {
      
      if (all.equal(twod, c("x", "y"))==TRUE){
        
        plot(LIE[[i]]$X[1], LIE[[i]]$Y[1], type="n", xlim=xlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=xlab1,...)
        r<-0
        for (k in 1:nrow(LIE[[i]])){
          if(LIE[[i]]$Bran[k]=="true"){
            r<-r+1
            a<-LIE[[i]]$Prec[k]
            b<-LIE[[i]]$Date[k]
            if (a!=0) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Y[k], col=colors[r],...)}
            m<-LIE[[i]]$Suiv[k]
            if (m!=0){
              while (LIE[[i]]$Apic[m]=="false"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r],...)
                m<-LIE[[i]]$Suiv[m]}
              if (LIE[[i]]$Apic[m]=="true"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Y[m], col=colors[r],...)}}}}}
      
      if (all.equal(twod, c("x", "z"))==TRUE){
        
        plot(LIE[[i]]$X[1], LIE[[i]]$Z[1], type="n", xlim=xlim1, ylim=zlim1, main=main1, ylab=zlab1, xlab=xlab1,...)
        r<-0
        for (k in 1:nrow(LIE[[i]])){
          if(LIE[[i]]$Bran[k]=="true"){
            r<-r+1
            a<-LIE[[i]]$Prec[k]
            b<-LIE[[i]]$Date[k]
            if (a!=0) {segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Z[a],x1=LIE[[i]]$X[k],y1=LIE[[i]]$Z[k], col=colors[r],...)}
            m<-LIE[[i]]$Suiv[k]
            if (m!=0){
              while (LIE[[i]]$Apic[m]=="false"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Z[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Z[m], col=colors[r],...)
                m<-LIE[[i]]$Suiv[m]}
              if (LIE[[i]]$Apic[m]=="true"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                segments(x0=LIE[[i]]$X[a],y0=LIE[[i]]$Z[a],x1=LIE[[i]]$X[m],y1=LIE[[i]]$Z[m], col=colors[r],...)}}}}}
      
      if (all.equal(twod, c("y", "z"))==TRUE){
        
        plot(LIE[[i]]$Z[1], LIE[[i]]$Y[1], type="n", xlim=zlim1, ylim=ylim1, main=main1, ylab=ylab1, xlab=zlab1,...)
        r<-0
        for (k in 1:nrow(LIE[[i]])){
          if(LIE[[i]]$Bran[k]=="true"){
            r<-r+1
            a<-LIE[[i]]$Prec[k]
            b<-LIE[[i]]$Date[k]
            if (a!=0) {segments(x0=LIE[[i]]$Z[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$Z[k],y1=LIE[[i]]$Y[k], col=colors[r],...)}
            m<-LIE[[i]]$Suiv[k]
            if (m!=0){
              while (LIE[[i]]$Apic[m]=="false"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                segments(x0=LIE[[i]]$Z[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$Z[m],y1=LIE[[i]]$Y[m], col=colors[r],...)
                m<-LIE[[i]]$Suiv[m]}
              if (LIE[[i]]$Apic[m]=="true"){
                a<-LIE[[i]]$Prec[m]
                b<-LIE[[i]]$Date[m]
                segments(x0=LIE[[i]]$Z[a],y0=LIE[[i]]$Y[a],x1=LIE[[i]]$Z[m],y1=LIE[[i]]$Y[m], col=colors[r],...)}}}}}}}}}
  
  names(rac)<-filenameslie
  names(tip)<-filenameslie
  outputresults<-list(root=rac, tip=tip)
  return(outputresults)}