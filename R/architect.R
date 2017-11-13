architect<-function(inputrac=NULL, inputtps=NULL, inputrsml=NULL, res=NULL, unitlength="px", rsml.date=NULL, rsml.connect=FALSE, vertical3d="y", fitter=TRUE){
  
  # Find out what kind of data input we have
  
  if (is.null(inputrac)==FALSE){
    if ("character" %in% class(inputrac) | "dartToTable" %in% class(inputrac)){} else {stop("inputrac must be a character string or a dartToTable object")}}
  
  if (is.null(inputrsml)==FALSE){
    if ("character" %in% class(inputrsml) | "rsmlToTable" %in% class(inputrsml)){} else {stop("inputrsml must be a character string or a rsmlToTable object")}}
  
  if (is.null(inputrac)==TRUE & is.null(inputrsml)==TRUE){algo<-"rawfiles"}
  if (is.null(inputrac)==TRUE & "rsmlToTable" %in% class(inputrsml)){algo<-"tables"}
  if (is.null(inputrac)==TRUE & mode(inputrsml)=="character"){algo<-"rawfiles"}
  
  if ("dartToTable" %in% class(inputrac) & is.null(inputrsml)==TRUE){algo<-"tables"}
  if ("dartToTable" %in% class(inputrac) & "rsmlToTable" %in% class(inputrsml)){algo<-"tables"}
  if ("dartToTable" %in% class(inputrac) & mode(inputrsml)=="character"){stop("inputrsml must be a rsmlToTable object")}
  
  if (mode(inputrac)=="character" & is.null(inputrsml)==TRUE){algo<-"rawfiles"}
  if (mode(inputrac)=="character" & "rsmlToTable" %in% class(inputrsml)){stop("inputrac must be a dartToTable object")}
  if (mode(inputrac)=="character" & mode(inputrsml)=="character"){algo<-"rawfiles"}
  
##Two possible algorithms for architect
##rawfiles if raw DART and RSML files have to be loaded (the raw files are stored in a folder)
##tables if input data are stored in rsmlToTable or dartToTable objects
  
#################
#First algorithm
#################
  
if (algo=="rawfiles"){
  
  # Errors interception
  
  if (is.null(inputrac)==TRUE & is.null(inputtps)==TRUE & is.null(inputrsml)==TRUE){stop("inputrac/inputps and/or inputrsml must be provided")}
  
  if (is.null(inputrac)==FALSE) {if (mode(inputrac)!="character"){stop("mode(inputrac) must be character")}}
  
  if (is.null(inputtps)==FALSE) {if (mode(inputtps)!="character"){stop("mode(inputtps) must be character")}}
  
  if (is.null(inputrac)==FALSE|is.null(inputtps)==FALSE){
    if (is.null(inputrac)==TRUE|is.null(inputtps)==TRUE){stop("If inputrac/inputtps is not NULL, inputtps/inputrac must be provided")}}
  
  if (is.null(inputrsml)==FALSE) {if (mode(inputrsml)!="character"){stop("mode(inputrsml) must be character")}}
  
  if (is.null(inputrac)==FALSE & is.null(res)==TRUE & unitlength!="px"){stop("If unitlength is not px, res must be specified")}
  if (is.null(res)==FALSE){
    if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
    if (res<=0|length(res)>1){stop("res must be a single positive value")}}
  
  if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
  if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
  
  if (is.null(rsml.date)==FALSE){
    if (is.character(rsml.date)==TRUE|is.numeric(rsml.date)==TRUE){} else {stop("If rsml.date is not NULL, rsml.date must be a character string or a positive numeric value")}
    if (is.numeric(rsml.date)==TRUE){if (rsml.date<=0|length(rsml.date)>1){stop("If mode(rsml.date) is numeric, rsml.date must be a single positive value")}}}
  
  if (mode(rsml.connect)!="logical"){stop("mode(rsml.connect) must be logical")}
  
  # Reading of DART and rsml files
  
  if (is.null(inputtps)==FALSE){
    filenames.tps<-list.files(path=inputtps, pattern="\\.tps$")
    path.tps<-rep(inputtps, length.out=length(filenames.tps))
    filenamestps<-sub(x=filenames.tps, pattern="\\.tps$", replacement="")
    message(paste("Number of DART tps files in inputtps:", length(filenames.tps), sep=" "))}
  
  if (is.null(inputrac)==FALSE){
    filenames.rac<-list.files(path=inputrac, pattern="\\.rac$")
    path.rac<-rep(inputrac, length.out=length(filenames.rac))
    filenamesrac<-sub(x=filenames.rac, pattern="\\.rac$", replacement="")
    message(paste("Number of DART rac files in inputrac:", length(filenames.rac), sep=" "))}
  
  if (is.null(inputrsml)==FALSE) {
    filenames.rsml<-list.files(path=inputrsml, pattern="\\.rsml$")
    path.rsml<-rep(inputrsml, length.out=length(filenames.rsml))
    filenamesrsml<-sub(x=filenames.rsml, pattern="\\.rsml$", replacement="")
    message(paste("Number of rsml files in inputrsml:", length(filenames.rsml), sep=" "))}
  
  if (is.null(inputrsml)==TRUE){
    if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
    if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}}
  else {
    if (is.null(inputrac)==TRUE){if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}
    else{
      if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
      if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}
      if (length(filenames.rsml)==0){stop("There is no rsml file in inputrsml")}}}

  if (is.null(inputrsml)==TRUE){ # Only DART files
  
      TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE)
      
      DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1)
      for (i in 1:length(DATA)) {
        colnames(DATA[[i]])<-c()
        colnames(DATA[[i]])[1]<-"Root"
        colnames(DATA[[i]])[2]<-"Mother"
        colnames(DATA[[i]])[3]<-"Ord"
        colnames(DATA[[i]])[4]<-"DBase"
        colnames(DATA[[i]])[5]<-"DApp"
        for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
      
      #Unit conversion DART files
      if (unitlength=="mm") {cunit<-(10*cm(1)/res)}
      if (unitlength=="cm") {cunit<-(cm(1)/res)}
      if (unitlength=="px") {cunit<-1}
      
      if (length(TIME)==1) {
        age<-list()
        obstot<-0
        for (i in 1:length(DATA)) {
          age[[i]]<-TIME[[1]]$Date
          obstot<-obstot+length(age[[i]])}
        for (i in 1:length(DATA)) {if(length(age[[i]])!=(ncol(DATA[[i]])-5)){stop("The number of observation dates between corresponding rac et tps files must be equal")}}}
      
      else {
        if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac files in inputrac and tps files in inputtps must be equal")}
        else {
          for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac files and their corresponding tps files must have the same name")}}
          age<-list()
          obstot<-0
          for (i in 1:length(TIME)) {
            age[[i]]<-TIME[[i]]$Date
            obstot<-obstot+length(age[[i]])}
          for (i in 1:length(DATA)) {if (length(age[[i]])!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}}
  
  else {
    if (is.null(inputrac)==TRUE){ # Only RSML files
      
      DATA<-list()
      age<-list()
      res1<-c()
      unitlength1<-c()
      filenamesrac<-c()
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect) # Read RSML files
      obstot<-0
      for (i in 1:length(RSML)){
        res1<-append(res1, rep(as.numeric(RSML[[i]]$resolution), length(RSML[[i]]$lie)))
        unitlength1<-append(unitlength1, rep(as.character(RSML[[i]]$length), length(RSML[[i]]$lie)))
        DATA<-append(DATA, RSML[[i]]$rac)
        length1<-length(RSML[[i]]$rac)
        
        l<-length(age)
        for (j in 1:length1) {
          age[[l+j]]<-RSML[[i]]$tps[[j]]$Date
          obstot<-obstot+length(age[[l+j]])}

        if (length1>1){
          num<-c(1:length1)
          filenamesrac[(length(filenamesrac)+1):(length(filenamesrac)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenamesrac[(length(filenamesrac)+1)]<-filenamesrsml[i]}}
      
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
    
    else { # DART and RSML files
      
      TIME<-lapply(paste(path.tps, "/", filenames.tps, sep=""), read.table, header=TRUE) # Read TPS files
      
      DATA<-lapply(paste(path.rac, "/", filenames.rac, sep=""), read.table, skip=1) # Read RAC files
      for (i in 1:length(DATA)) {
        colnames(DATA[[i]])<-c()
        colnames(DATA[[i]])[1]<-"Root"
        colnames(DATA[[i]])[2]<-"Mother"
        colnames(DATA[[i]])[3]<-"Ord"
        colnames(DATA[[i]])[4]<-"DBase"
        colnames(DATA[[i]])[5]<-"DApp"
        for (j in 6:ncol(DATA[[i]])-5) {colnames(DATA[[i]])[j+5]<-paste("Lengths", j, sep="")}}
      
      res1<-rep(res, length(filenames.rac))
      unitlength1<-rep(unitlength, length(filenames.rac))
      
      if (length(TIME)==1) {
        age<-list()
        obstot<-0
        for (i in 1:length(DATA)) {
          age[[i]]<-TIME[[1]]$Date
          obstot<-obstot+length(age[[i]])}
        for (i in 1:length(DATA)) {if (length(age[[i]])!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}
      else {
        if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac files in inputrac and tps files in inputtps must be equal")}
        else {
          for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac files and their corresponding tps files must have the same name")}}
          age<-list()
          obstot<-0
          for (i in 1:length(DATA)) {
            age[[i]]<-TIME[[i]]$Date
            obstot<-obstot+length(age[[i]])}
          for (i in 1:length(DATA)) {if (length(age[[i]])!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}
      
      RSML <- lapply(paste(path.rsml, "/", filenames.rsml, sep=""), rsmlToDART, final.date=rsml.date, connect=rsml.connect) # Read RSML files
      
      for (i in 1:length(RSML)){
        res1<-append(res1, rep(as.numeric(RSML[[i]]$resolution), length(RSML[[i]]$lie)))
        unitlength1<-append(unitlength1, rep(as.character(RSML[[i]]$length), length(RSML[[i]]$lie)))
        DATA<-append(DATA, RSML[[i]]$rac)
        length1<-length(RSML[[i]]$rac)
        
        l<-length(age)
        for (j in 1:length1) {
          age[[l+j]]<-RSML[[i]]$tps[[j]]$Date
          obstot<-obstot+length(age[[l+j]])}
        
        if (length1>1){
          num<-c(1:length1)
          filenamesrac[(length(filenamesrac)+1):(length(filenamesrac)+length1)]<-paste(rep(filenamesrsml[i], length.out=length1), num, sep="")}
        if (length1==1){
          filenamesrac[(length(filenamesrac)+1)]<-filenamesrsml[i]}}
      
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
  
  # Creating vectors and matrices for root architecture parameters calculation
  
  FileNames<-c()
  Time<-c()
  FirstOrderRootLength<-c()
  FirstOrderRootNumber<-c()
  TotalRootLength<-c()
  TotalLateralRootNumber<-c()
  TotalLateralRootLength<-c()
  LateralRootDensity<-c()
  GrowthRateFirstOrderRoot<-c()
  GrowthRateTotalRoot<-c()
  
  # Calculation of root architecture parameters
  
  k<-0
  maxord<-max(sapply(DATA, function(x){max(x[[3]])}))
  if (maxord>1) {latroot<-matrix(ncol=4*(maxord-1), nrow=obstot)}

  for(i in 1:length(DATA)){
    
    LastFirstOrderRootLength<-DATA[[i]][[ncol(DATA[[i]])]][DATA[[i]]$Ord==1]
    
    for (t in 1:(ncol(DATA[[i]])-5)){
      
      k<-k+1
      
      # Time
      Time[k]<-age[[i]][t]
      
      #Split data per root order
      data.split<-split(DATA[[i]][[paste("Lengths",t,sep="")]],as.factor(DATA[[i]]$Ord))
      
      # File names
      
      FileNames[k]<-filenamesrac[i]
      
      # Unit conversion
      
      if (is.null(inputrsml)==FALSE){cunit<-cunit1[i]}
      
      # Total root length
      
      TotalRootLength[k]<-sum(DATA[[i]][[paste("Lengths",t,sep="")]])*cunit
      
      # Growth rate of the root system
      
      if (t==1) {GrowthRateTotalRoot[k]<-TotalRootLength[k]/age[[i]][t]}
      if (t>1) {GrowthRateTotalRoot[k]<-(TotalRootLength[k]-TotalRootLength[k-1])/(age[[i]][t]-age[[i]][t-1])}
      
      # First-order root length
      
      FirstOrderRootLength[k]<-sum(data.split$'1')*cunit
      
      # Growth rate of the first-order root
      
      if (t==1) {GrowthRateFirstOrderRoot[k]<-FirstOrderRootLength[k]/age[[i]][t]}
      if (t>1) {GrowthRateFirstOrderRoot[k]<-(FirstOrderRootLength[k]-FirstOrderRootLength[k-1])/(age[[i]][t]-age[[i]][t-1])}
      
      # Total number of first-order roots
      
      FirstOrderRootNumber[k]<-sum(data.split$'1'!=0)
      
      # Total number of lateral roots
      
      if (FirstOrderRootLength[k]==0) {TotalLateralRootNumber[k]<-0}  else {TotalLateralRootNumber[k]<-sum(DATA[[i]][[paste("Lengths",t,sep="")]]!=0)-FirstOrderRootNumber[k]}
      
      # Total length of lateral roots
      
      TotalLateralRootLength[k]<-TotalRootLength[k]-FirstOrderRootLength[k]
      
      # Number, length and growth rate of lateral roots (by branching order)		
      
      if (maxord>1){
      
      for (l in 1:(maxord-1)){
        if (l<=max(DATA[[i]]$Ord)-1) {latroot[k,l]<-sum(data.split[[l+1]]!=0)} else {latroot[k,l]<-0}
        if (l<=max(DATA[[i]]$Ord)-1) {latroot[k,(l+(maxord-1))]<-sum(data.split[[l+1]])*cunit} else {latroot[k,(l+(maxord-1))]<-0}
        if (latroot[k,l]==0){latroot[k,(l+2*(maxord-1))]<-0} else {latroot[k,(l+2*(maxord-1))]<-latroot[k,(l+(maxord-1))]/latroot[k,l]}
        if (t==1) {latroot[k,(l+3*(maxord-1))]<-latroot[k,(l+(maxord-1))]/age[[i]][t]} else {latroot[k,(l+3*(maxord-1))]<-(latroot[k,(l+(maxord-1))]-latroot[k-1,(l+(maxord-1))])/(age[[i]][t]-age[[i]][t-1])}}}

      # Density of secondary roots on the first-order root
      
      if (FirstOrderRootLength[k]==0|maxord==1) {LateralRootDensity[k]<-0} else {LateralRootDensity[k]<-latroot[k,1]/FirstOrderRootLength[k]}}}	
  
  # Summary results in a data frame
  
  if (maxord>1) {outputresults<-data.frame(FileNames, Time, TotalRootLength, GrowthRateTotalRoot, FirstOrderRootLength, GrowthRateFirstOrderRoot, FirstOrderRootNumber, TotalLateralRootNumber, TotalLateralRootLength, latroot, LateralRootDensity)}
  else {outputresults<-data.frame(FileNames, Time, TotalRootLength, GrowthRateTotalRoot, FirstOrderRootLength, GrowthRateFirstOrderRoot, FirstOrderRootNumber, TotalLateralRootNumber, TotalLateralRootLength, LateralRootDensity)}
  
  if (maxord>1){
  LRnumberheading<-c()
  LRlengthheading<-c()
  LRmeanlengthheading<-c()
  LRgrowthrateheading<-c()
  
  for (h in 2:maxord){
    LRnumberheading[h-1]<-paste("N", h, "LR", sep="")
    LRlengthheading[h-1]<-paste("L", h, "LR", sep="")
    LRmeanlengthheading[h-1]<-paste("ML", h, "LR", sep="")
    LRgrowthrateheading[h-1]<-paste("GR", h, "L", sep="")}}
  
  if (maxord>1) {colnames(outputresults)<-c("FileName", "Time", "TRL", "GRTR", "L1R", "GR1R", "TN1R", "TNLR", "TLRL", t(LRnumberheading), t(LRlengthheading), t(LRmeanlengthheading), t(LRgrowthrateheading), "D2LR")}
  else {colnames(outputresults)<-c("FileName", "Time", "TRL", "GRTR", "L1R", "GR1R", "TN1R", "TNLR", "TLRL", "D2LR")}
  
  return(outputresults)}

##################
#Second algorithm
##################
  
  
if (algo=="tables"){

  if (vertical3d=="x"|vertical3d=="y"|vertical3d=="z") {} else {stop("vertical3d must be x, y, or z")}
  
  if (mode(fitter)!="logical"){stop("fitter must be logical")}
  
  if (is.null(inputrac)==TRUE & is.null(inputrsml)==FALSE) {maxord<-max(inputrsml$order)}
  if (is.null(inputrac)==FALSE & is.null(inputrsml)==TRUE) {maxord<-max(inputrac$order)}
  if (is.null(inputrac)==FALSE & is.null(inputrsml)==FALSE) {maxord<-max(c(inputrac$order, inputrsml$order))}
  
  #Processing DART files
  
  if (is.null(inputrac)==FALSE){
    
    n<-length(unique(paste(inputrac$file, inputrac$time, sep="")))
    
    if (fitter==TRUE){datadart<-data.frame(FileName=rep(NA, n), Time=rep(NA, n), TRL=rep(NA, n), GRTR=rep(NA, n), L1R=rep(NA, n), GR1R=rep(NA, n), TN1R=rep(NA, n), TNLR=rep(NA, n), TLRL=rep(NA, n), D2LR=rep(NA, n), Height=rep(NA, n), Width=rep(NA, n), Magnitude=rep(NA, n), Altitude=rep(NA, n), ExtPathLength=rep(NA, n))}
    
    else {datadart<-data.frame(FileName=rep(NA, n), Time=rep(NA, n), TRL=rep(NA, n), GRTR=rep(NA, n), L1R=rep(NA, n), GR1R=rep(NA, n), TN1R=rep(NA, n), TNLR=rep(NA, n), TLRL=rep(NA, n), D2LR=rep(NA, n), Height=rep(NA, n), Width=rep(NA, n))}
    
    if (maxord>1){latroot<-matrix(ncol=4*(maxord-1), nrow=n)}
    
    diameter<-matrix(ncol=maxord+1, nrow=n)
    
    files<-unique(inputrac$file)
    
    k<-0
    
    for (i in 1:length(files)){
      
      dates<-sort(unique(inputrac$time[inputrac$file==files[i]]))
      
      for (t in 1:length(dates)){
        
        xt<-inputrac[inputrac$file==files[i],]
        
        k<-k+1
        
        #File name
        datadart$FileName[k]<-xt$file[1]
        
        #Time
        datadart$Time[k]<-dates[t]
        
        #TRL (total root length)
        datadart$TRL[k]<-sum(xt$length[xt$time<=dates[t]])
        
        #GRTR (growth rate of the root system)
        GR<-aggregate(xt$growth, by=list(root=xt$root, order=xt$order, time=xt$time), max)
        datadart$GRTR[k]<-sum(GR$x[GR$time==dates[t]])

        #L1R (total first-order root length)
        datadart$L1R[k]<-sum(xt$length[xt$time<=dates[t] & xt$order==1])
        
        #GR1R (growth rate of first-order roots)
        datadart$GR1R[k]<-sum(GR$x[GR$time==dates[t] & GR$order==1])
        
        #TN1R (number of first-order roots)
        datadart$TN1R[k]<-length(unique(xt$root[xt$time<=dates[t] & xt$order==1]))
        
        #TNLR (number of lateral roots)
        datadart$TNLR[k]<-length(unique(xt$root[xt$time<=dates[t] & xt$order>1]))
        
        #TLRL (total lateral root length)
        datadart$TLRL[k]<-sum(xt$length[xt$time<=dates[t] & xt$order>1])
        
        #D2LR (density of second-order roots on the first-order root)
        datadart$D2LR[k]<-length(unique(xt$root[xt$time<=dates[t] & xt$order==2]))/datadart$L1R[k]
          
        #Latroot (number, length, mean length, growth rate)
        
        if (maxord>1){
          
          for (l in 1:(maxord-1)){
            
            latroot[k, l]<-length(unique(xt$root[xt$time<=dates[t] & xt$order==(l+1)])) #Number of lateral roots
            latroot[k, l+(maxord-1)]<-sum(xt$length[xt$time<=dates[t] & xt$order==(l+1)]) #Length of lateral roots
            if (latroot[k,l]==0) {latroot[k, l+2*(maxord-1)]<-0} else {latroot[k, l+2*(maxord-1)]<-latroot[k, l+(maxord-1)]/latroot[k, l]} #Mean lateral root length
            latroot[k, l+3*(maxord-1)]<-sum(GR$x[GR$time==dates[t] & GR$order==(l+1)])}}
    
        #Height
        if (vertical3d=="y") {datadart$Height[k]<-abs(max(xt$y2[xt$time<=dates[t]])-min(xt$y1[xt$time<=dates[t]]))}
        if (vertical3d=="x") {datadart$Height[k]<-abs(max(xt$x2[xt$time<=dates[t]])-min(xt$x1[xt$time<=dates[t]]))}
        if (vertical3d=="z") {datadart$Height[k]<-abs(max(xt$z2[xt$time<=dates[t]])-min(xt$z1[xt$time<=dates[t]]))}
    
        #Width
        if (vertical3d=="y") {
          
          widthx<-abs(max(xt$x2[xt$time<=dates[t]])-min(xt$x2[xt$time<=dates[t]]))
          widthz<-abs(max(xt$z2[xt$time<=dates[t]])-min(xt$z2[xt$time<=dates[t]]))
          datadart$Width[k]<-max(c(widthx, widthz))}
    
        if (vertical3d=="x") {
          
          widthy<-abs(max(xt$y2[xt$time<=dates[t]])-min(xt$y2[xt$time<=dates[t]]))
          widthz<-abs(max(xt$z2[xt$time<=dates[t]])-min(xt$z2[xt$time<=dates[t]]))
          datadart$Width[k]<-max(c(widthy, widthz))}
    
        if (vertical3d=="z") {
          
          widthy<-abs(max(xt$y2[xt$time<=dates[t]])-min(xt$y2[xt$time<=dates[t]]))
          widthx<-abs(max(xt$x2[xt$time<=dates[t]])-min(xt$x2[xt$time<=dates[t]]))
          datadart$Width[k]<-max(c(widthy, widthx))}
        
        #Fitter topological indices
        
        if (fitter==TRUE){
          
        xt<-xt[xt$time<=dates[t],]
        apex<-aggregate(xt$blength, by=list(root=xt$root), max)
        apicindex<-as.vector(apply(apex, 1, function(x){which(xt$root==x["root"] & xt$blength==x["x"])}))
        xt$apic[apicindex]<-"true"
        branindex<-which(xt$bran=="true")
        
        #Magnitude
        datadart$Magnitude[k]<-sum(xt$bran=="true")
        
        #Altitude and external path length
        xt$pathlength<-rep(1, nrow(xt))
        
        if (nrow(xt)>1){
        
            for (l in 1:length(branindex)){
              
              testbran<-which(xt$x1==xt$x2[branindex[l]] & xt$y1==xt$y2[branindex[l]] & xt$z1==xt$z2[branindex[l]])

              if (length(testbran)==0) {} else {
                
                if (length(testbran)>=2) {
                  xt$pathlength[testbran]<-xt$pathlength[branindex[l]]+1
                  index<-which(xt$bran[testbran]=="false")
                  suiv<-testbran[index]}
                else {
                  xt$pathlength[testbran]<-xt$pathlength[branindex[l]]
                  suiv<-testbran}}
              
              while(xt$apic[suiv]=="false"){
                
                testbran<-which(xt$x1==xt$x2[suiv] & xt$y1==xt$y2[suiv] & xt$z1==xt$z2[suiv])
                
                if (length(testbran)>=2) {
                  xt$pathlength[testbran]<-xt$pathlength[suiv]+1
                  index<-which(xt$bran[testbran]=="false")
                  suiv<-testbran[index]}
                else {
                  xt$pathlength[testbran]<-xt$pathlength[suiv]
                  suiv<-testbran}}}
        
        datadart$Altitude[k]<-max(xt$pathlength)
        datadart$ExtPathLength[k]<-sum(xt$pathlength[xt$apic=="true"])}
        
        else{
          
          datadart$Altitude[k]<-1
          datadart$ExtPathLength[k]<-1}}}}
    
        # Results in a dataframe
        
        diameter<-as.data.frame(diameter)
        colnames(diameter)<-c(paste(rep("MD", maxord), c(1:maxord), sep=""), "MDLR")
        
        if (maxord>1){
          
          latroot<-as.data.frame(latroot)
        
        for (l in 1:(maxord-1)){
          colnames(latroot)[l]<-paste("N", l+1, "LR", sep="")
          colnames(latroot)[l+(maxord-1)]<-paste("L", l+1, "LR", sep="")
          colnames(latroot)[l+2*(maxord-1)]<-paste("ML", l+1, "LR", sep="")
          colnames(latroot)[l+3*(maxord-1)]<-paste("GR", l+1, "L", sep="")}
        
        if (fitter==TRUE){datadart<-data.frame(datadart[,1:9], as.data.frame(latroot), as.data.frame(diameter), datadart[,10:15])}  
            
        else {datadart<-data.frame(datadart[,1:9], as.data.frame(latroot), as.data.frame(diameter), datadart[,10:12])}}}
  
  #Processing RSML files
  
  if (is.null(inputrsml)==FALSE){
    
    inputrsml$diameter3<-inputrsml$diameter1
    inputrsml$diameter3[inputrsml$apic=="true"]<-inputrsml$diameter1[inputrsml$apic=="true"]+inputrsml$diameter2[inputrsml$apic=="true"]
    
    n<-length(unique(paste(inputrsml$file, inputrsml$plant, inputrsml$time, sep="")))
    
    if (fitter==TRUE){datarsml<-data.frame(FileName=rep(NA, n), Time=rep(NA, n), TRL=rep(NA, n), GRTR=rep(NA, n), L1R=rep(NA, n), GR1R=rep(NA, n), TN1R=rep(NA, n), TNLR=rep(NA, n), TLRL=rep(NA, n), D2LR=rep(NA, n), Height=rep(NA, n), Width=rep(NA, n), Magnitude=rep(NA, n), Altitude=rep(NA, n), ExtPathLength=rep(NA, n))}
    
    else {datarsml<-data.frame(FileName=rep(NA, n), Time=rep(NA, n), TRL=rep(NA, n), GRTR=rep(NA, n), L1R=rep(NA, n), GR1R=rep(NA, n), TN1R=rep(NA, n), TNLR=rep(NA, n), TLRL=rep(NA, n), D2LR=rep(NA, n), Height=rep(NA, n), Width=rep(NA, n))}
    
    if (maxord>1){latroot<-matrix(ncol=4*(maxord-1), nrow=n)}
    
    diameter<-matrix(ncol=maxord+1, nrow=n)
    
    files<-unique(inputrsml$file)
    
    k<-0
    
    for (i in 1:length(files)){
      
      plants<-unique(inputrsml$plant[inputrsml$file==files[i]])
      
      for (j in 1:length(plants)){
      
      dates<-sort(unique(inputrsml$time[inputrsml$file==files[i] & inputrsml$plant==plants[j]]))
      
      for (t in 1:length(dates)){
        
        xt<-inputrsml[inputrsml$file==files[i] & inputrsml$plant==plants[j],] #x is a subset of the rsmlToTable object
        
        k<-k+1
        
        #File name
        datarsml$FileName[k]<-paste(xt$file[1], xt$plant[1], sep="_")
        
        #Time
        datarsml$Time[k]<-dates[t]
        
        #TRL (total root length)
        datarsml$TRL[k]<-sum(xt$length[xt$time<=dates[t]])
        
        #GRTR (growth rate of the root system)
        GR<-aggregate(xt$growth, by=list(root=xt$root, order=xt$order, time=xt$time), max)
        datarsml$GRTR[k]<-sum(GR$x[GR$time==dates[t]])
        
        #L1R (total first-order root length)
        datarsml$L1R[k]<-sum(xt$length[xt$time<=dates[t] & xt$order==1])
        
        #GR1R (growth rate of first-order roots)
        datarsml$GR1R[k]<-sum(GR$x[GR$time==dates[t] & GR$order==1])

        #TN1R (number of first-order roots)
        datarsml$TN1R[k]<-length(unique(xt$root[xt$time<=dates[t] & xt$order==1]))
        
        #TNLR (number of lateral roots)
        datarsml$TNLR[k]<-length(unique(xt$root[xt$time<=dates[t] & xt$order>1]))
        
        #TLRL (total lateral root length)
        datarsml$TLRL[k]<-sum(xt$length[xt$time<=dates[t] & xt$order>1])
        
        #D2LR (density of second-order roots on the first-order root)
        datarsml$D2LR[k]<-length(unique(xt$root[xt$time<=dates[t] & xt$order==2]))/datarsml$L1R[k]
        
        #Latroot (number, length, mean length, growth rate)
        
        if (maxord>1){
          
          for (l in 1:(maxord-1)){
            
            latroot[k, l]<-length(unique(xt$root[xt$time<=dates[t] & xt$order==(l+1)])) #Number of lateral roots
            latroot[k, l+(maxord-1)]<-sum(xt$length[xt$time<=dates[t] & xt$order==(l+1)]) #Length of lateral roots
            if (latroot[k,l]==0) {latroot[k, l+2*(maxord-1)]<-0} else {latroot[k, l+2*(maxord-1)]<-latroot[k, l+(maxord-1)]/latroot[k, l]} #Mean lateral root length
            latroot[k, l+3*(maxord-1)]<-sum(GR$x[GR$time==dates[t] & GR$order==(l+1)])}}
    
        #Diameter
        
        if (t==length(dates)){#Dianeter in rsml file is for the last observation date
        
        maxordfile<-max(inputrsml$order[inputrsml$file==files[i]])
        
        for (l in 1:maxordfile){diameter[k,l]<-sum(inputrsml$diameter3[inputrsml$file==files[i] & inputrsml$order==l])/(nrow(inputrsml[inputrsml$file==files[i] & inputrsml$order==l,])+sum(inputrsml$apic[inputrsml$file==files[i] & inputrsml$order==l]=="true"))}
        
        diameter[k, ncol(diameter)]<-sum(inputrsml$diameter3[inputrsml$file==files[i] & inputrsml$order>1])/(nrow(inputrsml[inputrsml$file==files[i] & inputrsml$order>1,])+sum(inputrsml$apic[inputrsml$file==files[i] & inputrsml$order>1]=="true"))} #Mean diameter for all lateral roots
        
        #Height
        if (vertical3d=="y") {datarsml$Height[k]<-abs(max(xt$y2[xt$time<=dates[t]])-min(xt$y1[xt$time<=dates[t]]))}
        if (vertical3d=="x") {datarsml$Height[k]<-abs(max(xt$x2[xt$time<=dates[t]])-min(xt$x1[xt$time<=dates[t]]))}
        if (vertical3d=="z") {datarsml$Height[k]<-abs(max(xt$z2[xt$time<=dates[t]])-min(xt$z1[xt$time<=dates[t]]))}
        
        #Width
        if (vertical3d=="y") {
          
          widthx<-abs(max(xt$x2[xt$time<=dates[t]])-min(xt$x2[xt$time<=dates[t]]))
          widthz<-abs(max(xt$z2[xt$time<=dates[t]])-min(xt$z2[xt$time<=dates[t]]))
          datarsml$Width[k]<-max(c(widthx, widthz))}
        
        if (vertical3d=="x") {
          
          widthy<-abs(max(xt$y2[xt$time<=dates[t]])-min(xt$y2[xt$time<=dates[t]]))
          widthz<-abs(max(xt$z2[xt$time<=dates[t]])-min(xt$z2[xt$time<=dates[t]]))
          datarsml$Width[k]<-max(c(widthy, widthz))}
        
        if (vertical3d=="z") {
          
          widthy<-abs(max(xt$y2[xt$time<=dates[t]])-min(xt$y2[xt$time<=dates[t]]))
          widthx<-abs(max(xt$x2[xt$time<=dates[t]])-min(xt$x2[xt$time<=dates[t]]))
          datarsml$Width[k]<-max(c(widthy, widthx))}
        
        #Fitter topological indices
        
        if (fitter==TRUE){
          
          xt<-xt[xt$time<=dates[t],]
          apex<-aggregate(xt$blength, by=list(root=xt$root), max)
          apicindex<-as.vector(apply(apex, 1, function(x){which(xt$root==x["root"] & xt$blength==x["x"])}))
          xt$apic[apicindex]<-"true"
          branindex<-which(xt$bran=="true") #Equal to the number of roots
          
          #Magnitude
          datarsml$Magnitude[k]<-sum(xt$bran=="true")
          
          #Altitude and external path length
          xt$pathlength<-rep(1, nrow(xt))
          
          if (nrow(xt)>1){
            
            for (l in 1:length(branindex)){ #For each root
              
              root<-xt$root[branindex[l]]
              
              testbran<-which(xt$x1==xt$x2[branindex[l]] & xt$y1==xt$y2[branindex[l]] & xt$z1==xt$z2[branindex[l]]) #Is it a crossing?
              testbran<-testbran[which(xt$root[testbran]==root | xt$parentroot[testbran]==root)] #Select segments based on root ID and parentrootID

              if (length(testbran)==0) {} else {
                
                if (length(testbran)>=2) {
                  xt$pathlength[testbran]<-xt$pathlength[branindex[l]]+1
                  index<-which(xt$bran[testbran]=="false")
                  suiv<-testbran[index]}
                else {
                  xt$pathlength[testbran]<-xt$pathlength[branindex[l]]
                  suiv<-testbran}}
              
              while(xt$apic[suiv]=="false"){
                
                testbran<-which(xt$x1==xt$x2[suiv] & xt$y1==xt$y2[suiv] & xt$z1==xt$z2[suiv]) #Is it a crossing?
                testbran<-testbran[which(xt$root[testbran]==root | xt$parentroot[testbran]==root)] #Select segments based on root ID and parentrootID

                if (length(testbran)>=2) {
                  xt$pathlength[testbran]<-xt$pathlength[suiv]+1
                  index<-which(xt$bran[testbran]=="false")
                  suiv<-testbran[index]}
                else {
                  xt$pathlength[testbran]<-xt$pathlength[suiv]
                  suiv<-testbran}}}
            
            datarsml$Altitude[k]<-max(xt$pathlength)
            datarsml$ExtPathLength[k]<-sum(xt$pathlength[xt$apic=="true"])}
          
          else{
            
            datarsml$Altitude[k]<-1
            datarsml$ExtPathLength[k]<-1}}}}}
        
        # Results in a dataframe
    
        diameter<-as.data.frame(diameter)
        colnames(diameter)<-c(paste(rep("MD", maxord), c(1:maxord), sep=""), "MDLR")
        
        if (maxord>1){
          
          latroot<-as.data.frame(latroot)
          
          for (l in 1:(maxord-1)){
            colnames(latroot)[l]<-paste("N", l+1, "LR", sep="")
            colnames(latroot)[l+(maxord-1)]<-paste("L", l+1, "LR", sep="")
            colnames(latroot)[l+2*(maxord-1)]<-paste("ML", l+1, "LR", sep="")
            colnames(latroot)[l+3*(maxord-1)]<-paste("GR", l+1, "L", sep="")}
          
          if (fitter==TRUE) {datarsml<-data.frame(datarsml[,1:9], as.data.frame(latroot), as.data.frame(diameter), datarsml[,10:15])}
          
          else {datarsml<-data.frame(datarsml[,1:9], as.data.frame(latroot), as.data.frame(diameter), datarsml[,10:12])}}}
  
  if (is.null(inputrac)==TRUE & is.null(inputrsml)==FALSE) {return(datarsml)}
  if (is.null(inputrsml)==TRUE & is.null(inputrac)==FALSE) {return(datadart)}
  
  if (is.null(inputrsml)==FALSE & is.null(inputrac)==FALSE){
    
    if (ncol(datadart)==ncol(datarsml)){return(rbind(datadart, datarsml))}
    
    else{
      
      if (ncol(datadart)>ncol(datarsml)){datarsml<-data.frame(datarsml, Magnitude=rep(NA, nrow(datarsml)), Altitude=rep(NA, nrow(datarsml)), ExtPathLength=rep(NA, nrow(datarsml)))}
      
      if (ncol(datadart)<ncol(datarsml)){datadart<-data.frame(datadart, Magnitude=rep(NA, nrow(datadart)), Altitude=rep(NA, nrow(datadart)), ExtPathLength=rep(NA, nrow(datadart)))}
      
      return(rbind(datadart, datarsml))}}}}