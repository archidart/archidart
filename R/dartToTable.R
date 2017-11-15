dartToTable<-function(inputrac, inputlie, inputtps, res=NULL, unitlength="px", unitangle="d", fitter=FALSE){

  if (mode(inputrac)!="character"){stop("mode(inputrac) must be character")}
  
  if (mode(inputlie)!="character"){stop("mode(inputlie) must be character")}
  
  if (mode(inputtps)!="character"){stop("mode(inputtps) must be character")}
  
  if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
  if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
  
  if (is.null(res)==FALSE){
    if (mode(res)!="numeric"){stop("mode(res) must be numeric")}
    if (res<=0){stop("res must be a positive value")}}
  
  if (mode(unitangle)!="character"){stop("mode(unitangle) must be character")}
  if(unitangle=="d"|unitangle=="r") {} else {stop("unitangle must be either d (degrees) or r (radians)")}
  
  if (mode(fitter)!="logical"){stop("fitter must be logical")}
  
  #Load DART files
  
  filenames.rac<-list.files(path=inputrac, pattern="\\.rac$")
  path.rac<-rep(inputrac, length.out=length(filenames.rac))
  filenamesrac<-sub(x=filenames.rac, pattern="\\.rac$", replacement="")
  message(paste("Number of DART rac files in inputrac:", length(filenames.rac), sep=" "))
  
  filenames.lie<-list.files(path=inputlie, pattern="\\.lie$")
  path.lie<-rep(inputlie, length.out=length(filenames.lie))
  filenameslie<-sub(x=filenames.lie, pattern="\\.lie$", replacement="")
  message(paste("Number of DART lie files in inputlie:", length(filenames.lie), sep=" "))
  
  filenames.tps<-list.files(path=inputtps, pattern="\\.tps$")
  path.tps<-rep(inputtps, length.out=length(filenames.tps))
  filenamestps<-sub(x=filenames.tps, pattern="\\.tps$", replacement="")
  message(paste("Number of DART tps files in inputtps:", length(filenames.tps), sep=" "))
  
  if (length(filenames.rac)==0){stop("There is no rac file in inputrac")}
  if (length(filenames.tps)==0){stop("There is no tps file in inputtps")}
  if (length(filenames.lie)==0){stop("There is no lie file in inputlie")}
  
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
    for (i in 1:length(DATA)) {if(length(TIME[[1]]$Date)!=(ncol(DATA[[i]])-5)){stop("The number of observation dates between corresponding rac et tps files must be equal")}}
    TIME<-rep(list(TIME[[1]]), length(DATA))}
  
  else {
    if (length(TIME)!=length(DATA)) {stop("If there is more than one tps file in inputtps, the number of rac/lie files in inputrac/inputlie and tps files in inputtps must be equal")}
    else {
      for (i in 1:length(DATA)) {if (filenamesrac[i]!=filenamestps[i]) {stop("Input rac/lie files and their corresponding tps files must have the same name")}}
      for (i in 1:length(DATA)) {if (length(TIME[[i]]$Date)!=(ncol(DATA[[i]])-5)) {stop("The number of observation dates between corresponding rac et tps files must be equal")}}}}
  
 
  nodes<-0 #nodes is the number of rows (sum rows of each lie)

  for (i in 1:length(LIE)){nodes<-nodes+nrow(LIE[[i]])}
  
  # Vertical direction vector
  
  dirvert<-c(0,1)
  
  if (unitangle=="r") {cunitangle<-1}
  if (unitangle=="d") {cunitangle<-180/pi}
  
  #Construct dart table (1 line per segment)
  
  table<-data.frame(file=rep(NA, nodes), root=rep(NA, nodes), dbase=rep(NA, nodes), dbasecum=rep(NA, nodes), order=rep(NA, nodes), time=rep(NA, nodes), deltaage=rep(NA, nodes), bran=rep(NA, nodes), apic=rep(NA, nodes), x1=rep(NA, nodes), y1=rep(NA, nodes), z1=rep(NA, nodes), x2=rep(NA, nodes), y2=rep(NA, nodes), z2=rep(NA, nodes), length=rep(NA, nodes), blength=rep(NA, nodes), orientation=rep(NA, nodes), growth=rep(NA, nodes), geodesic=rep(NA, nodes))

  rowsintable<-0
  n<-0 #n is the number of lie files
  
  for (i in 1:length(DATA)){ #For each root system
    
      n<-n+1
      lie<-LIE[[i]]
      rac<-DATA[[i]]
      tps<-TIME[[i]]
      
      #Add dbasecum column in rac file
      
        rac$CumDBase<-rep(NA, nrow(rac))
      
        for (l in 1:nrow(rac)){
        
          if (rac$Ord[l]==1) {rac$CumDBase[l]<-0}
          else {rac$CumDBase[l]<-rac$CumDBase[rac$Mother[l]+1]+rac$DBase[l]}}
      
      s<-0 #Count number of segments added to table

      brantrue<-which(lie$Bran=="true")
      
      for (l in 1:length(brantrue)){ #For each root in dart file
          
          cumsum<-0
          
          r<-brantrue[l]
          
          if (rac$Ord[l]==1 & lie$Bran[r]=="true") {r<-lie$Suiv[r]}
          
          while (lie$Apic[r]=="false"){
          
            s<-s+1
            prec<-lie$Prec[r]
            table$file[rowsintable+s]<-filenamesrac[i]
            table$root[rowsintable+s]<-l
            table$dbase[rowsintable+s]<-rac$DBase[l]*cunit
            table$dbasecum[rowsintable+s]<-rac$CumDBase[l]*cunit
            table$order[rowsintable+s]<-rac$Ord[l]
            table$time[rowsintable+s]<-tps$Date[lie$Date[r]]
            
            if (lie$Bran[prec]=="true" & rac$Ord[l]==1) {table$bran[rowsintable+s]<-as.character(lie$Bran[prec])} else {table$bran[rowsintable+s]<-as.character(lie$Bran[r])}
            
            table$apic[rowsintable+s]<-as.character(lie$Apic[r])

            if (tps$Date[lie$Date[r]]!=min(tps$Date)) {table$deltaage[rowsintable+s]<-tps$Date[lie$Date[r]]-tps$Date[lie$Date[r]-1]} else {table$deltaage[rowsintable+s]<-tps$Date[lie$Date[r]]}
            
            table$x1[rowsintable+s]<-lie$X[prec]*cunit
            table$y1[rowsintable+s]<-lie$Y[prec]*cunit
            table$z1[rowsintable+s]<-0
            
            table$x2[rowsintable+s]<-lie$X[r]*cunit
            table$y2[rowsintable+s]<-lie$Y[r]*cunit
            table$z2[rowsintable+s]<-0
  
            table$length[rowsintable+s]<-distance2D(x1=lie$X[prec], y1=lie$Y[prec], x2=lie$X[r], y2=lie$Y[r])*cunit
            cumsum<-cumsum+table$length[rowsintable+s]
            
            table$blength[rowsintable+s]<-cumsum
            
            dirsegment<-c(lie$X[r]-lie$X[prec], lie$Y[r]-lie$Y[prec])*cunit
            table$orientation[rowsintable+s]<-acos(as.numeric(dirvert%*%dirsegment)/table$length[rowsintable+s])*cunitangle
            
            r<-lie$Suiv[r]}
          
          if (lie$Apic[r]=="true"){
            
            s<-s+1
            prec<-lie$Prec[r]
            table$file[rowsintable+s]<-filenamesrac[i]
            table$root[rowsintable+s]<-l
            table$dbase[rowsintable+s]<-rac$DBase[l]*cunit
            table$dbasecum[rowsintable+s]<-rac$CumDBase[l]*cunit
            table$order[rowsintable+s]<-rac$Ord[l]
            table$time[rowsintable+s]<-tps$Date[lie$Date[r]]
            
            if (lie$Bran[prec]=="true" & rac$Ord[l]==1) {table$bran[rowsintable+s]<-as.character(lie$Bran[prec])} else {table$bran[rowsintable+s]<-as.character(lie$Bran[r])}
            
            table$apic[rowsintable+s]<-as.character(lie$Apic[r])
            
            if (tps$Date[lie$Date[r]]!=min(tps$Date)) {table$deltaage[rowsintable+s]<-tps$Date[lie$Date[r]]-tps$Date[lie$Date[r]-1]} else {table$deltaage[rowsintable+s]<-tps$Date[lie$Date[r]]}
            
            table$x1[rowsintable+s]<-lie$X[prec]*cunit
            table$y1[rowsintable+s]<-lie$Y[prec]*cunit
            table$z1[rowsintable+s]<-0
            
            table$x2[rowsintable+s]<-lie$X[r]*cunit
            table$y2[rowsintable+s]<-lie$Y[r]*cunit
            table$z2[rowsintable+s]<-0
            
            table$length[rowsintable+s]<-distance2D(x1=lie$X[prec], y1=lie$Y[prec], x2=lie$X[r], y2=lie$Y[r])*cunit
            cumsum<-cumsum+table$length[rowsintable+s]
            
            table$blength[rowsintable+s]<-cumsum
            
            dirsegment<-c(lie$X[r]-lie$X[prec], lie$Y[r]-lie$Y[prec])*cunit
            table$orientation[rowsintable+s]<-acos(as.numeric(dirvert%*%dirsegment)/table$length[rowsintable+s])*cunitangle}}
      
      rowsintable<-rowsintable+s
      
      #Add parentroot column and reorder table
      index<-which(table$file==filenamesrac[i])
      table[index, "parentroot"]<-rac$Mother[table$root[index]]+1
      table<-table[,c("file", "root", "dbase", "dbasecum", "order", "parentroot", "time", "deltaage", "bran", "apic", "x1", "y1", "z1", "x2", "y2", "z2", "length", "blength", "orientation", "growth", "geodesic")]}
  
  index<-which(is.na(table$file)==TRUE)
  table<-table[-index,] #Remove lines with NA values

  rownames(table)<-c(1:nrow(table))
  
  table$geodesic<-table$dbasecum+table$blength
  
  #Calculate growth rate of each segment and fill growth column
  sum<-aggregate(table$length, by=list(table$file, table$time, table$root), sum)
  colnames(sum)<-c("file", "time", "root", "length")
  
  index<-as.vector(apply(table, 1, function(x){which(sum$file==x["file"] & sum$time==as.numeric(x["time"]) & sum$root==as.numeric(x["root"]))}))
  length1<-sum$length[index]
  table$growth<-length1/table$deltaage
  
  table<-table[,-c(3,4,8)]
  
  if (fitter==TRUE){
  
  #Compute Fitter topological indices
  
  table$magnitude<-rep(0, nrow(table))
  table$pathlength<-rep(1, nrow(table))
  
  RSlevels<-unique(table$file)
  n<-length(RSlevels)
  
  for (i in 1:n){#For each root system in the table
    
    subtable<-table[table$file==RSlevels[i],] #Create a subtable subset

    apicindex<-which(subtable$apic=="true") #apicindex and branindex should have the same length
    branindex<-which(subtable$bran=="true")

    for (l in 1:length(apicindex)){#For each apic point of a root system
      
      #Magnitude
      subtable$magnitude[apicindex[l]]<-1
      
      root<-subtable$root[apicindex[l]]
      parentroot<-subtable$parentroot[apicindex[l]]
      
      indexprec<-which(subtable$x2==subtable$x1[apicindex[l]] & subtable$y2==subtable$y1[apicindex[l]] & subtable$z2==subtable$z1[apicindex[l]])
      if (length(indexprec)>1) {indexprec<-indexprec[which(subtable$root[indexprec]==root | subtable$root[indexprec]==parentroot)]}
      
      if (length(indexprec)>0){
        
        root<-subtable$root[indexprec]
        parentroot<-subtable$parentroot[indexprec]
        
        subtable$magnitude[indexprec]<-subtable$magnitude[indexprec]+1
        
        while(subtable$order[indexprec]>=1){
          
          segment1<-which(subtable$x2==subtable$x1[indexprec] & subtable$y2==subtable$y1[indexprec] & subtable$z2==subtable$z1[indexprec])
          if (length(segment1)>1) {indexprec<-segment1[which(subtable$root[segment1]==root | subtable$root[segment1]==parentroot)]} else {indexprec<-segment1}
          if (length(indexprec)==0){break}
          root<-subtable$root[indexprec]
          parentroot<-subtable$parentroot[indexprec]
          
          subtable$magnitude[indexprec]<-subtable$magnitude[indexprec]+1}}
    
      #Path length
      
      root<-subtable$root[branindex[l]]
      
      testbran<-which(subtable$x1==subtable$x2[branindex[l]] & subtable$y1==subtable$y2[branindex[l]] & subtable$z1==subtable$z2[branindex[l]]) #Is it a crossing?
      if (length(testbran)>1) {testbran<-testbran[which(subtable$root[testbran]==root | subtable$parentroot[testbran]==root)]} #Select segments based on root ID and parentroot ID
      
      if (length(testbran)==0) {} else{
        
        if (length(testbran)>=2) {
          subtable$pathlength[testbran]<-subtable$pathlength[branindex[l]]+1
          index<-which(subtable$bran[testbran]=="false")
          suiv<-testbran[index]}
        else {
          subtable$pathlength[testbran]<-subtable$pathlength[branindex[l]]
          suiv<-testbran}}
      
      while(subtable$apic[suiv]=="false"){
        
        testbran<-which(subtable$x1==subtable$x2[suiv] & subtable$y1==subtable$y2[suiv] & subtable$z1==subtable$z2[suiv]) #Is it a crossing?
        if (length(testbran)>1) {testbran<-testbran[which(subtable$root[testbran]==root | subtable$parentroot[testbran]==root)]} #Select segments based on root ID and parentroot ID
        
        if (length(testbran)>=2) {
          subtable$pathlength[testbran]<-subtable$pathlength[suiv]+1
          index<-which(subtable$bran[testbran]=="false")
          suiv<-testbran[index]}
        else {
          subtable$pathlength[testbran]<-subtable$pathlength[suiv]
          suiv<-testbran}}}
    
    table$magnitude[table$file==RSlevels[i]]<-subtable$magnitude
    table$pathlength[table$file==RSlevels[i]]<-subtable$pathlength}}

class(table)<-c("data.frame", "dartToTable")
return(table)}