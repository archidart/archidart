rsmlToTable<-function(inputrsml, unitlength="px", rsml.date=NULL, rsml.connect=TRUE, vertical3d="y", unitangle="d", fitter=FALSE, show.progress=FALSE){

  if (mode(inputrsml)!="character"){stop("mode(inputrsml) must be character")}
  
  if (mode(unitlength)!="character"){stop("mode(unitlength) must be character")}
  if (unitlength=="px"|unitlength=="mm"|unitlength=="cm") {} else {stop("unitlength must be either px (pixels), mm (millimeters) or cm (centimeters)")}
  
  if (is.null(rsml.date)==FALSE){
    if (is.character(rsml.date)==TRUE|is.numeric(rsml.date)==TRUE){} else {stop("If rsml.date is not NULL, rsml.date must be a character string or a positive numeric value")}
    if (is.numeric(rsml.date)==TRUE){if (rsml.date<=0|length(rsml.date)>1){stop("If mode(rsml.date) is numeric, rsml.date must be a single positive value")}}}
  
  if (mode(rsml.connect)!="logical"){stop("mode(rsml.connect) must be logical")}
  
  if (vertical3d=="x"|vertical3d=="y"|vertical3d=="z") {} else {stop("vertical3d must be x, y, or z")}
  
  if (mode(unitangle)!="character"){stop("mode(unitangle) must be character")}
  if(unitangle=="d"|unitangle=="r") {} else {stop("unitangle must be either d (degrees) or r (radians)")}
  
  if (mode(fitter)!="logical"){stop("fitter must be logical")}
  
  if (mode(show.progress)!="logical"){stop("show.progress must be logical")}
  
  if (fitter==TRUE & rsml.connect==FALSE){stop("If fitter is TRUE, rsml.connect must be TRUE too")}
  
  #Load rsml files
  
  filenames.rsml<-list.files(path=inputrsml, pattern="\\.rsml$")
  filenamesrsml<-sub(x=filenames.rsml, pattern="\\.rsml$", replacement="")
  message(paste("Number of rsml files in inputrsml:", length(filenames.rsml), sep=" "))
  
  if (show.progress==TRUE) {pb<-txtProgressBar(min=1, max=length(filenames.rsml), style=3)}
  
  TABLE<-vector("list", length(filenamesrsml))
  
  for (f in 1:length(filenames.rsml)){
    
    if (show.progress==TRUE) {setTxtProgressBar(pb, f)}
    
    RSML<-rsmlToDART(rsml.path=paste(inputrsml, "/", filenames.rsml[f], sep=""), final.date=rsml.date, connect=rsml.connect)
    
    nodes<-0 #nodes is the number of rows (sum rows of each lie)
    
    res1<-rep(as.numeric(RSML$resolution), length(RSML$lie))
    filenamesrac<-rep(filenamesrsml[f], length(RSML$lie))
    unitlength1<-rep(as.character(RSML$length), length(RSML$lie))
    for (j in 1:length(RSML$lie)) {nodes<-nodes+nrow(RSML$lie[[j]])}
  
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
    
    if (unitlength=="px"){cunit1[i]<-1}}
  
  # Vertical direction vector
  
  if (vertical3d=="x") {dirvert<-c(1,0,0)}
  if (vertical3d=="y") {dirvert<-c(0,1,0)}
  if (vertical3d=="z") {dirvert<-c(0,0,1)}
  
  if (unitangle=="r") {cunitangle<-1}
  if (unitangle=="d") {cunitangle<-180/pi}
  
  #Construct rsml table (1 line per segment)
  
  table<-matrix(nrow=nodes, ncol=22)
  
  rowsintable<-0
  n<-0 #n is the number of lie files
  
  for (j in 1:length(RSML$lie)){ #For each plant in rsml
      
      n<-n+1
      lie<-RSML$lie[[j]]
      rac<-RSML$rac[[j]]
      tps<-RSML$tps[[j]]
      
      if ("Z" %in% colnames(lie)) {} else {dirvert<-c(0,1)}
      
      #Add dbasecum column in rac file
      
      rac$CumDBase<-rep(NA, nrow(rac))
    
      for (l in 1:nrow(rac)){
      
        if (rac$Ord[l]==1) {rac$CumDBase[l]<-0}
        else {rac$CumDBase[l]<-rac$CumDBase[rac$Mother[l]+1]+rac$DBase[l]}}
      
      s<-0 #Count number of segments added to table
      
      for (l in 1:nrow(lie)){ #For each line in rsml
        
        if (lie$Prec[l]!=0){
          
          s<-s+1
          prec<-lie$Prec[l]
          table[rowsintable+s, 1]<-f #file
          table[rowsintable+s, 2]<-j #plant
          table[rowsintable+s, 3]<-lie$root[l] #root
          table[rowsintable+s, 4]<-rac$CumDBase[lie$root[l]]*cunit1[n] #dbasecum
          table[rowsintable+s, 5]<-lie$ord[l] #order
          table[rowsintable+s, 6]<-tps$Date[lie$Date[l]] #time
          
          if (lie$Bran[prec]=="true" & lie$ord[l]==1) {table[rowsintable+s, 8]<-1} #bran
          else {
            if (lie$Bran[l]=="true") {table[rowsintable+s, 8]<-1}
            if (lie$Bran[l]=="false") {table[rowsintable+s, 8]<-0}}
          
          if (lie$Apic[l]=="true") {table[rowsintable+s, 9]<-1} #apic
          if (lie$Apic[l]=="false") {table[rowsintable+s, 9]<-0}
          
          if (tps$Date[lie$Date[l]] != min(tps$Date)) {table[rowsintable + s, 7] <- tps$Date[lie$Date[l]] - tps$Date[lie$Date[l] - 1]}
          else {table[rowsintable + s, 7] <- tps$Date[lie$Date[l]]} #deltaage
          
          table[rowsintable+s, 10]<-lie$X[prec]*cunit1[n] #x1
          table[rowsintable+s, 11]<-lie$Y[prec]*cunit1[n] #y1
          if ("Z" %in% colnames(lie)) {table[rowsintable+s, 12]<-lie$Z[prec]*cunit1[n]} else {table[rowsintable+s, 12]<-0} #z1
          
          table[rowsintable+s, 13]<-lie$X[l]*cunit1[n] #x2
          table[rowsintable+s, 14]<-lie$Y[l]*cunit1[n] #y2
          if ("Z" %in% colnames(lie)) {table[rowsintable+s, 15]<-lie$Z[l]*cunit1[n]} else {table[rowsintable+s, 15]<-0} #z2
          
          if (lie$Bran[l]=="true") {table[rowsintable+s, 16]<-lie$diameter[l]*cunit1[n]} else {table[rowsintable+s, 16]<-lie$diameter[prec]*cunit1[n]} #diameter1
          table[rowsintable+s, 17]<-lie$diameter[l]*cunit1[n] #diameter2
          if ("Z" %in% colnames(lie)) {table[rowsintable+s, 18]<-distance3D(x1=lie$X[prec], y1=lie$Y[prec], z1=lie$Z[prec], x2=lie$X[l], y2=lie$Y[l], z2=lie$Z[l])*cunit1[n]} else {table[rowsintable+s, 18]<-distance2D(x1=lie$X[prec], y1=lie$Y[prec], x2=lie$X[l], y2=lie$Y[l])*cunit1[n]} #length
          table[rowsintable+s, 19]<-lie$dist[l]*cunit1[n] #blength
          dirsegment<-c(lie$X[l]-lie$X[prec], lie$Y[l]-lie$Y[prec], lie$Z[l]-lie$Z[prec])*cunit1[n]
          table[rowsintable+s, 20]<-acos(as.numeric(dirvert%*%dirsegment)/table[rowsintable+s, 18])*cunitangle}} #orientation
      
      rowsintable<-rowsintable+s}
  
  index<-which(is.na(table[,1])==TRUE)
  table<-table[-index,] #Remove lines with NA values
  
  #Calculate growth rate of each segment and fill growth column

  sum<-aggregate(table[,18], by=list(table[,2], table[,6], table[,3]), sum)
  colnames(sum)<-c("plant", "time", "root", "length")
  
  index<-as.vector(apply(table, 1, function(x){which(sum$plant==as.numeric(x[2]) & sum$time==as.numeric(x[6]) & sum$root==as.numeric(x[3]))}))
  length1<-sum$length[index]
  table[,21]<-length1/table[,7]
  
  #Geodesic distance
  if (rsml.connect==TRUE) {table[,22]<-table[,4]+table[,19]} else {table<-table[,-22]}
  
  #Check if segments have length=0
  
  if (sum(table[,18]==0)>0){
    
    index<-which(table[,18]==0)
    
    for (i in 1:length(index)){
      
      indexsuiv<-which(table[,3]==table[index[i], 3] & table[,10]==table[index[i], 13] & table[,11]==table[index[i], 14] & table[,12]==table[index[i], 15])
      indexsuiv<-indexsuiv[indexsuiv!=index[i]]
      table[indexsuiv, 8]<-table[index[i], 8]
      table[indexsuiv, 16]<-table[index[i], 16]}
    
    table<-table[-index,]
    
    rownames(table)<-c(1:nrow(table))}
  
  table<-table[,-c(4,7)] #Remove dbasecum and deltaage
  
  #Add parentroot column
  parentroot<-rac$Mother[table[,3]]+1
  table<-cbind(table, parentroot)

  #Fitter
  if (rsml.connect==TRUE & fitter==TRUE) {table<-fitter(table)}
  
  #Reorder columns
  if (rsml.connect==TRUE){
    if (fitter==FALSE){table<-table[,c(1:4, 21, 5:20)]}
    else {table<-table[,c(1:4, 21, 5:20, 22:23)]}}
  else{table<-table[,c(1:4, 20, 5:19)]}
  
  #Store table in a list
  TABLE[[f]]<-table
  }
  
  #Convert list to a data frame
  TABLE<-as.data.frame(do.call(rbind, TABLE))
  
  if (rsml.connect==TRUE) {
    if (fitter==FALSE) {colnames(TABLE)<-c("file", "plant", "root", "order", "parentroot", "time", "bran", "apic", "x1", "y1", "z1", "x2", "y2", "z2", "diameter1", "diameter2", "length", "blength", "orientation", "growth", "geodesic")}
    else {colnames(TABLE)<-c("file", "plant", "root", "order", "parentroot", "time", "bran", "apic", "x1", "y1", "z1", "x2", "y2", "z2", "diameter1", "diameter2", "length", "blength", "orientation", "growth", "geodesic", "magnitude", "pathlength")}}
  
  else {colnames(TABLE)<-c("file", "plant", "root", "order", "parentroot", "time", "bran", "apic", "x1", "y1", "z1", "x2", "y2", "z2", "diameter1", "diameter2", "length", "blength", "orientation", "growth")}
  
  TABLE$file<-filenamesrsml[TABLE$file]
  
  TABLE$bran[which(TABLE$bran==0)]<-"false"
  TABLE$bran[which(TABLE$bran==1)]<-"true"
  
  TABLE$apic[which(TABLE$apic==0)]<-"false"
  TABLE$apic[which(TABLE$apic==1)]<-"true"

rownames(TABLE)<-c(1:nrow(TABLE))
class(TABLE)<-c("data.frame", "rsmlToTable")
return(TABLE)}