perhomology<-function(x, show.progress=FALSE){
  
  #x must be an object of class "rsmlToTable" or "dartToTable"
  if ("rsmlToTable" %in% class(x) | "dartToTable" %in% class(x)) {} else {stop("x must be a rsmlToTable or dartToTable object")}
  
  if (mode(show.progress)!="logical"){stop("show.progress must be logical")}
  
  if ("dartToTable" %in% class(x)) {x$rootsystem<-x$file}
  if ("rsmlToTable" %in% class(x)) {x$rootsystem<-paste(x$file, x$plant, sep="_")}
  RSlevels<-unique(x$rootsystem)
  
  n<-length(RSlevels) #Number of root systems in x
  
  #Compute persistent homology (geodesic distance function)
  
  if (show.progress==TRUE) {pb<-txtProgressBar(min=1, max=n, style=3)}
  
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
      
      root<-table$root[apicindex[l]]
      parentroot<-table$parentroot[apicindex[l]]
      
      table$hzero[apicindex[l]]<-l
      
      results[[i]][l,2]<-table$geodesic[apicindex[l]]
      
      results[[i]][l,3]<-table$geodesic[apicindex[l]]-table$length[apicindex[l]]
      
      indexprec<-which(table$x2==table$x1[apicindex[l]] & table$y2==table$y1[apicindex[l]] & table$z2==table$z1[apicindex[l]])
      
      if (length(indexprec)>1) {
        
        indexprec<-indexprec[which(table$root[indexprec]==root | table$root[indexprec]==parentroot)]
        
        if (length(indexprec)>1){
          indexprec<-indexprec[which(is.na(table[indexprec, "hzero"])==TRUE)]
          
          if (length(indexprec)>1){
            max.geodesic<-max(table[indexprec, "geodesic"])
            indexprec<-indexprec[which(table[indexprec, "geodesic"]==max.geodesic)]}}}
      
      if (length(indexprec)>0){
        
        root<-table$root[indexprec]
        parentroot<-table$parentroot[indexprec]
        
        while(is.na(table$hzero[indexprec])==TRUE){
          
          table$hzero[indexprec]<-l
          
          results[[i]][l,3]<-table$geodesic[indexprec]-table$length[indexprec]
          
          segment1<-which(table$x2==table$x1[indexprec] & table$y2==table$y1[indexprec] & table$z2==table$z1[indexprec])
          
          #print(paste("root: ", root, ", parentroot: ", parentroot, ", segment1-pre: ", segment1, sep=""))
          
          #if (length(segment1)>1){print(table[segment1,])}
          
          if (length(segment1)>1) {
            indexprec<-segment1[which(table$root[segment1]==root | table$root[segment1]==parentroot)]
            
            if (length(indexprec)>1){ 
              segment1<-segment1[which(is.na(table[segment1, "hzero"])==TRUE)]
              
              if (length(segment1)>1){ #Choose segment based on geodesic distance
                max.geodesic<-max(table[segment1, "geodesic"])
                segment1<-segment1[which(table[segment1, "geodesic"]==max.geodesic)]}
              
              if (length(segment1)>1){ #Then, choose the last segment of the list
                segment1<-segment1[length(segment1)]}
              
              indexprec<-segment1}
          } 
          
          else {indexprec<-segment1}
          
          #print(paste("indexprec: ", indexprec, sep=""))
          
          if (length(indexprec)==0){break}
          root<-table$root[indexprec]
          parentroot<-table$parentroot[indexprec]}}}
    
    results[[i]]<-results[[i]][order(results[[i]][,3], decreasing=FALSE),]
    class(results[[i]])<-c("matrix", "barcode")}
  
  names(results)<-as.character(RSlevels)
  class(results)<-c("perhomology", "list")
  return(results)}