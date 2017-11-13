perhomology<-function(x, FUN="geodesic", show.progress=FALSE){
  
  #x must be an object of class "rsmlToTable" or "dartToTable"
  if ("rsmlToTable" %in% class(x) | "dartToTable" %in% class(x)) {} else {stop("x must be a rsmlToTable or dartToTable object")}
  
  if (FUN=="geodesic"|FUN=="depth") {} else {stop("FUN must be either geodesic or depth")}
  
  if (mode(show.progress)!="logical"){stop("show.progress must be logical")}
  
  if ("dartToTable" %in% class(x)) {x$rootsystem<-x$file}
  if ("rsmlToTable" %in% class(x)) {x$rootsystem<-paste(x$file, x$plant, sep="_")}
  RSlevels<-levels(as.factor(x$rootsystem))
  
  n<-length(RSlevels) #Number of root systems in x
  
  #Compute persistent homology (geodesic distance function)
  
  if (show.progress==TRUE) {pb<-txtProgressBar(min=1, max=n, style=3)}
  
  if (FUN=="geodesic"){
  
  results<-list()
  
  for (i in 1:n){#For each root system in the table
    
    if (show.progress==TRUE) {setTxtProgressBar(pb, i)}
    
    table<-x[x$rootsystem==RSlevels[i],] #Create a table subset for computing persistent homology
    results[[i]]<-matrix(nrow=sum(table$apic=="true"), ncol=3) #Create matrix to store the results (birth and death of each axis)
    colnames(results[[i]])<-c("dimension", "birth", "death") #Rename columns
    table<-table[order(table$geodesic, decreasing=TRUE),]#Order table lines by geodesic values
    table$hzero<-rep(NA, nrow(table))#Create a new column in table to store the ID number of each zero order homology
    
    apicindex<-which(table$apic=="true")
    results[[i]][,1]<-rep(0, nrow(results[[i]]))
    
    for (l in 1:length(apicindex)){#For each apic point of a root system
      
      table$hzero[apicindex[l]]<-l
      
      results[[i]][l,2]<-table$geodesic[apicindex[l]]
      results[[i]][l,3]<-table$geodesic[apicindex[l]]-table$length[apicindex[l]]

      if (table$order[apicindex[l]]==1 & table$bran[apicindex[l]]=="true") {break}
      
      testbran<-which(table$x1==table$x1[apicindex[l]] & table$y1==table$y1[apicindex[l]] & table$z1==table$z1[apicindex[l]])
      
      if (length(testbran)==2){#We are at a crossing
        
        segment1<-which(table$x2==table$x1[apicindex[l]] & table$y2==table$y1[apicindex[l]] & table$z2==table$z1[apicindex[l]])
        segment2<-testbran[testbran!=apicindex[l]]
        geo1<-table$geodesic[segment1]
        geo2<-table$geodesic[segment2]
        if (geo1>geo2){indexprec<-segment2}
        if (geo1<geo2){indexprec<-segment1}}
      
      else {indexprec<-which(table$x2==table$x1[apicindex[l]] & table$y2==table$y1[apicindex[l]] & table$z2==table$z1[apicindex[l]])}
      
      while(is.na(table$hzero[indexprec])==TRUE){
        
        table$hzero[indexprec]<-l
        
        results[[i]][l,3]<-table$geodesic[indexprec]-table$length[indexprec]
        
        if (table$order[indexprec]==1 & table$bran[indexprec]=="true") {break}

        testbran<-which(table$x1==table$x1[indexprec] & table$y1==table$y1[indexprec] & table$z1==table$z1[indexprec])
        
        if (length(testbran)==2){#We are at a crossing
          
          segment1<-which(table$x2==table$x1[indexprec] & table$y2==table$y1[indexprec] & table$z2==table$z1[indexprec])
          segment2<-testbran[testbran!=indexprec]
          geo1<-table$geodesic[segment1]
          geo2<-table$geodesic[segment2]
          if (geo1>geo2){indexprec<-segment2}
          if (geo1<geo2){indexprec<-segment1}}
        
        else {indexprec<-which(table$x2==table$x1[indexprec] & table$y2==table$y1[indexprec] & table$z2==table$z1[indexprec])}}}
    
    results[[i]]<-results[[i]][order(results[[i]][,3], decreasing=FALSE),]
    class(results[[i]])<-c("matrix", "barcode")}}
  
  #Compute persistent homology (depth distance function)
  
  if (FUN=="depth"){
    
    results<-list()
    
    for (i in 1:n){#For each root system in the table
      
      if (show.progress==TRUE) {setTxtProgressBar(pb, i)}
      
      table<-x[x$rootsystem==RSlevels[i],] #Create a table subset for computing persistent homology
      results[[i]]<-matrix(nrow=sum(table$apic=="true"), ncol=3) #Create matrix to store the results (birth and death of each axis)
      colnames(results[[i]])<-c("dimension", "birth", "death") #Rename columns
      table<-table[order(table$y2, decreasing=TRUE),]#Order table lines by depth values
      table$hzero<-rep(NA, nrow(table))#Create a new column in table to store the ID number of each zero order homology
      
      apicindex<-which(table$apic=="true")
      results[[i]][,1]<-rep(0, nrow(results[[i]]))
      
      for (l in 1:length(apicindex)){#For each apic point of a root system
        
        table$hzero[apicindex[l]]<-l
        
        results[[i]][l,2]<-table$y2[apicindex[l]] #birth
        results[[i]][l,3]<-table$y1[apicindex[l]] #death

        if (table$order[apicindex[l]]==1 & table$bran[apicindex[l]]=="true") {break}
        
        testbran<-which(table$x1==table$x1[apicindex[l]] & table$y1==table$y1[apicindex[l]] & table$z1==table$z1[apicindex[l]])
        
        if (length(testbran)==2){#We are at a crossing
          
          segment1<-which(table$x2==table$x1[apicindex[l]] & table$y2==table$y1[apicindex[l]] & table$z2==table$z1[apicindex[l]])
          segment2<-testbran[testbran!=apicindex[l]]
          depth1<-table$y2[segment1]
          depth2<-table$y2[segment2]
          if (depth1>depth2){indexprec<-segment2}
          if (depth1<depth2){indexprec<-segment1}}
        
        else {indexprec<-which(table$x2==table$x1[apicindex[l]] & table$y2==table$y1[apicindex[l]] & table$z2==table$z1[apicindex[l]])}
        
        while(is.na(table$hzero[indexprec])==TRUE){
          
          table$hzero[indexprec]<-l
          
          results[[i]][l,3]<-table$y1[indexprec]
          
          if (table$order[indexprec]==1 & table$bran[indexprec]=="true") {break}
          
          testbran<-which(table$x1==table$x1[indexprec] & table$y1==table$y1[indexprec] & table$z1==table$z1[indexprec])
          
          if (length(testbran)==2){#We are at a crossing
            
            segment1<-which(table$x2==table$x1[indexprec] & table$y2==table$y1[indexprec] & table$z2==table$z1[indexprec])
            segment2<-testbran[testbran!=indexprec]
            depth1<-table$y2[segment1]
            depth2<-table$y2[segment2]
            if (depth1>depth2){indexprec<-segment2}
            if (depth1<depth2){indexprec<-segment1}}
          
          else {indexprec<-which(table$x2==table$x1[indexprec] & table$y2==table$y1[indexprec] & table$z2==table$z1[indexprec])}}}
      
      results[[i]]<-results[[i]][order(results[[i]][,3], decreasing=FALSE),]
      class(results[[i]])<-c("matrix", "barcode")}}
  
  names(results)<-as.character(RSlevels)
  class(results)<-c("perhomology", "list")
  return(results)}