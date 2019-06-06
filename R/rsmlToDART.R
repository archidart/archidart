rsmlToDART <- function(rsml.path, final.date, connect){
  
  # Create the plant object
  rsml <- xmlToList(xmlParse(rsml.path))
  plants <- rsml$scene
  resolution<-rsml$metadata$resolution
  unitlength1<-rsml$metadata$unit
  
  # Find unittime in 'property-definition'
  
  unittime1<-NULL
  
  if (is.character(final.date)==TRUE) {
    property<-rsml$metadata$'property-definition'
    
    if (is.null(property)==FALSE){
      for (i in 1:length(property)){if (property[[i]]$label==final.date){unittime1<-as.character(property[[i]]$unit)}}}
    
    if (is.null(unittime1)==TRUE){message(paste("No time unit found in ", sub(basename(rsml.path), pattern=".rsml", replacement=""), " metadata (property-definition)", sep=""))}}
  
  # Create LIE and RAC files for each root system 
  
  n <- 0 # Initialise number of LIE/RAC files (as 1 RSML file can contain more than 1 root system)
  lie.all<-list()
  rac.all<-list()
  tps.all<-list()
  
  for (r0 in plants){# For each plant
    
    #Calculate total number of nodes in scene (! add additional nodes if connect=TRUE)
    
    nodes<-length(grep(x=names(rapply(r0, length, how="unlist")), pattern="point"))
    
    #Calculate total number of roots in scene
    
    roots<-length(grep(x=names(rapply(r0, length, how="unlist")), pattern="root..attrs"))

    if (connect==TRUE){nodes<-nodes+roots}
    
    n<-n+1 #Add one unit for each root system
    r <- 0 # Initialise the number of roots consisting a root system
    timeserie<-FALSE #Does the rsml file contain time series data? (Ex: the root system age for each node)
    
    #For each plant root system, create LIE, RAC and TPS files
    lie<-matrix(nrow=nodes, ncol=13)
    rac<-matrix(nrow=roots, ncol=6)
    lie.lines<-0
    
    for (r1 in r0){ # For each first-order root
      
      if (class(r1)=="list"){
        
        ns<-r1$geometry$polyline
        
        # Extract plant age when nodes appeared for LIE and TPS files
        
        age=NULL
        diameter=NULL
        
        if (is.character(final.date)==TRUE & "functions" %in% names(r1)){
          age<-r1$functions
          for (i in 1:length(age)){
            
            if (final.date %in% age[[i]]$.attrs) {
            timeserie<-TRUE
            time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
            if (time[1]<time[2]){time[1]<-time[2]}}}}
        
        if ("functions" %in% names(r1)){
          age<-r1$functions
          for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
        
        length1<-length(ns)
        r<-r+1 #Add one unit for each root

        if (r==1){#For the first first-order root only
        
        currentMother<-r
        
        #c(Num, Date, Bran, Apic, Prec, Suiv)
        if (timeserie==FALSE) {lie[1:length1,1:6]<-c(c(1:length1),rep(1, length1),rep(0, length1), rep(0, length1), 0:(length1-1), 2:(length1+1))}
        if (timeserie==TRUE) {lie[1:length1,1:6]<-c(c(1:length1),time,rep(0, length1), rep(0, length1), 0:(length1-1), 2:(length1+1))}
        lie[1:length1,11]<-r
        if (is.null(diameter)==TRUE){} else {lie[1:length1, 12]<-diameter}
        lie[1:length1,13]<-1
        
        #c(X,Y,Z)
        lie[1:length1,7]<-sapply(ns, xnodes)
        lie[1:length1,8]<-sapply(ns, ynodes)
        if (length(ns[[1]])==3) {lie[1:length1,9]<-sapply(ns, znodes)} else {lie[1:length1,9]<-0}
 
        #c(dist)
        lie[1:length1, 10]<-c(0, cumsum(sqrt((diff(lie[1:length1, 7]))^2+(diff(lie[1:length1, 8]))^2+(diff(lie[1:length1, 9]))^2)))
        
        #For the first point of the first order root
        if (timeserie==FALSE) {lie[1,c(2:3)]<-c(0,1)} else {lie[1,3]<-1}
        start1<-as.numeric(lie[1, 1])
        
        #For the last point of the first order root
        lie[length1,c(4,6)]<-c(1,0)
        stop1<-as.numeric(lie[length1, 1])
        
        # Fill RAC file for the first order root
        #c(Root, Mother, Ord, DBase, DApp, Length)
        cumulDist<-sum(sqrt((diff(lie[1:length1, 7]))^2+(diff(lie[1:length1, 8]))^2+(diff(lie[1:length1, 9]))^2))
        rac[r,1:6]<-c(r,-1,1,0,0,cumulDist)

        lie.lines<-lie.lines+length1}
        
        else {#For the other first-order roots only
          
          currentMother<-r
          
          #c(Num, Date, Bran, Apic, Prec, Suiv)
          if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length1),1:6]<-c(c((lie.lines+1):(lie.lines+length1)),rep(1, length1),rep(0, length1), rep(0, length1), c(0,(lie.lines+1):(lie.lines+length1-1)), (lie.lines+2):(lie.lines+length1+1))}
          if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length1),1:6]<-c(c((lie.lines+1):(lie.lines+length1)),time,rep(0, length1), rep(0, length1), c(0,(lie.lines+1):(lie.lines+length1-1)), (lie.lines+2):(lie.lines+length1+1))}
          lie[(lie.lines+1):(lie.lines+length1),11]<-r
          if (is.null(diameter)==TRUE){} else {lie[(lie.lines+1):(lie.lines+length1),12]<-diameter}
          lie[(lie.lines+1):(lie.lines+length1),13]<-1
          
          #c(X,Y,Z)
          lie[(lie.lines+1):(lie.lines+length1),7]<-sapply(ns, xnodes)
          lie[(lie.lines+1):(lie.lines+length1),8]<-sapply(ns, ynodes)
          if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length1),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length1),9]<-0}

          #c(dist)
          lie[(lie.lines+1):(lie.lines+length1), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length1), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length1), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length1), 9]))^2)))
          
          #For the first point of the first order root
          if (timeserie==FALSE) {lie[lie.lines+1,c(2:3)]<-c(0,1)} else {lie[lie.lines+1,3]<-1}
          start1<-as.numeric(lie[lie.lines+1, 1])
          
          #For the last point of the first order root
          lie[(lie.lines+length1),c(4,6)]<-c(1,0)
          stop1<-as.numeric(lie[lie.lines+length1, 1])
          
          # Fill RAC file for the first order root
          #c(Root, Mother, Ord, DBase, DApp, Length)
          cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length1), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length1), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length1), 9]))^2))
          rac[r,1:6]<-c(r,-1,1,0,0,cumulDist)
          
          lie.lines<-lie.lines+length1}
        
        #----------------------------------------------------------------------------------------------------
        
        # If there are lateral roots
        
        if ("root" %in% names(r1)){
          
          for (r2 in r1){# For each 2-order root
            
            if ("geometry" %in% names(r2)){
              
              r<-r+1
              ns <- r2$geometry$polyline
              length2<-length(ns)

              age=NULL
              diameter=NULL
              if (is.character(final.date)==TRUE & "functions" %in% names(r2)){
                age<-r2$functions
                for (i in 1:length(age)){
                  if (final.date %in% age[[i]]$.attrs) {
                    time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
                    if (time[1]<time[2]){time[1]<-time[2]}}}}
              
              if ("functions" %in% names(r2)){
                age<-r2$functions
                for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
              
              #c(Num, Date, Bran, Apic, Prec, Suiv)
              if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length2),1:6]<-c((lie.lines+1):(lie.lines+length2), rep(1, length2), rep(0, length2), rep(0, length2), lie.lines:(lie.lines+length2-1), (lie.lines+2):(lie.lines+length2+1))}
              if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length2),1:6]<-c((lie.lines+1):(lie.lines+length2), time, rep(0, length2), rep(0, length2), lie.lines:(lie.lines+length2-1), (lie.lines+2):(lie.lines+length2+1))}
              lie[(lie.lines+1):(lie.lines+length2),11]<-r
              if(is.null(diameter)==TRUE) {} else {lie[(lie.lines+1):(lie.lines+length2),12]<-diameter}
              lie[(lie.lines+1):(lie.lines+length2),13]<-2 #Order 2
              
              #c(X,Y,Z)
              lie[(lie.lines+1):(lie.lines+length2),7]<-sapply(ns, xnodes)
              lie[(lie.lines+1):(lie.lines+length2),8]<-sapply(ns, ynodes)
              if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length2),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length2),9]<-0}
              
              #c(dist)
              lie[(lie.lines+1):(lie.lines+length2), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length2), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length2), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length2), 9]))^2)))
              
              # Search the closest point on the mother root (calculate DBase)
              
              parentnode<-which(lie[start1:stop1, 7]==lie[lie.lines+1, 7] & lie[start1:stop1, 8]==lie[lie.lines+1, 8] & lie[start1:stop1, 9]==lie[lie.lines+1, 9])
              
              if (length(parentnode)>0){
                
                dbase<-lie[start1+parentnode[1]-1, 10]
                if (connect==TRUE){
                  lie[lie.lines+1,5]<-lie[start1+parentnode[1]-1, 1]
                  dist1<-0}}
              
              else { #If no physical connection between parent and daughter root. Interpolate new point or find closest point.
              
              alldist<-sqrt((lie[start1:(stop1-1), 7]-lie[lie.lines+1, 7])^2+(lie[start1:(stop1-1), 8]-lie[lie.lines+1, 8])^2+(lie[start1:(stop1-1), 9]-lie[lie.lines+1, 9])^2)

              scalx<-diff(lie[start1:stop1, 7])*(lie[start1:(stop1-1), 7]-lie[lie.lines+1, 7])
              scaly<-diff(lie[start1:stop1, 8])*(lie[start1:(stop1-1), 8]-lie[lie.lines+1, 8])
              scalz<-diff(lie[start1:stop1, 9])*(lie[start1:(stop1-1), 9]-lie[lie.lines+1, 9])
              d2<-diff(lie[start1:stop1, 7])^2+diff(lie[start1:stop1, 8])^2+diff(lie[start1:stop1, 9])^2
              t<-(-(scalx+scaly+scalz)/d2)

              if (length(which(t>=0 & t<=1))==0){
                
                index<-which(alldist==min(alldist))
                dbase<-lie[start1+index-1, 10]
                if (connect==TRUE){
                  lie[lie.lines+1,5]<-lie[start1+index-1, 1] #Update Prec
                  dist1<-min(alldist)
                  lie[(lie.lines+1):(lie.lines+length2), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length2), 10]}}

              else {
              
              t[t<0]<-NA
              t[t>1]<-NA
              xn<-diff(lie[start1:stop1, 7])*t+lie[start1:(stop1-1), 7]
              yn<-diff(lie[start1:stop1, 8])*t+lie[start1:(stop1-1), 8]
              zn<-diff(lie[start1:stop1, 9])*t+lie[start1:(stop1-1), 9]
              dist1<-sqrt((xn-lie[lie.lines+1, 7])^2+(yn-lie[lie.lines+1, 8])^2+(zn-lie[lie.lines+1, 9])^2)
              
              if (sum(is.na(dist1)==T)>0) {
                index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                dist1<-min(dist1, na.rm=T)} 
              else {
                dist1<-min(dist1)
                index<-as.numeric(match(min(dist1), dist1))}
              
              if (dist1>min(alldist)){
                index<-which(alldist==min(alldist))
                dist1<-min(alldist)
                dbase<-lie[start1+index-1, 10]
                if (connect==TRUE){
                  lie[lie.lines+1,5]<-lie[start1+index-1, 1] #Update Prec
                  dist1<-min(alldist)
                  lie[(lie.lines+1):(lie.lines+length2), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length2), 10]}}
              
              else{
                
                lie[(lie.lines+1):(lie.lines+length2), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length2), 10]
                dbase<-lie[start1+index-1, 10]+distance3D(x1=lie[start1+index-1, 7], x2=xn[index], y1=lie[start1+index-1, 8], y2=yn[index], z1=lie[start1+index-1, 9], z2=zn[index])
    
              if (connect==TRUE){
              
              lie[lie.lines+length2+1,1:13]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie[start1+index-1, 11], NA, 1)
              lie<-lie[order(lie[, 11], lie[, 10]),]
              length1<-length1+1
              lie.lines<-lie.lines+1
              stop1<-stop1+1
              pos<-match(NA, lie[1:(lie.lines),1])
              lie[1:(lie.lines+length2),5]<-match(lie[1:(lie.lines+length2),5], lie[1:(lie.lines+length2),1])
              lie[1:(lie.lines+length2),6]<-match(lie[1:(lie.lines+length2),6], lie[1:(lie.lines+length2),1])
              lie[pos, 5]<-pos-1
              lie[pos, 6]<-pos+1
              lie[pos, 2]<-lie[lie[pos, 6], 2]
              lie[pos+1, 5]<-pos
              lie[pos-1, 6]<-pos
              lie[which(is.na(lie[1:(lie.lines),5])==TRUE), 5]<-0
              lie[which(is.na(lie[1:(lie.lines),6])==TRUE), 6]<-0
              lie[1:(lie.lines+length2),1]<-match(lie[1:(lie.lines+length2),1], lie[1:(lie.lines+length2),1])
              lie[lie.lines+1,5]<-pos
              
              #Calculate diameter for interpolated point (linear function)
              if (is.null(diameter)==FALSE){
                prec<-lie[pos, 5]
                suiv<-lie[pos, 6]
                slope<-(lie[suiv, 12]-lie[prec, 12])/distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[suiv, 7], y2=lie[suiv, 8], z2=lie[suiv, 9])
                intercept<-lie[prec, 12]
                lie[pos, 12]<-intercept+slope*distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[pos, 7], y2=lie[pos, 8], z2=lie[pos, 9])}}}}}
              
              if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
              
              lie[lie.lines+1,3]<-1
              start2<-as.numeric(lie[lie.lines+1, 1])
              
              # Change Suiv and Apic values for the last point of a lateral root
              lie[lie.lines+length2,c(4,6)]<-c(1, 0)
              stop2<-as.numeric(lie[lie.lines+length2, 1])
              
              # Fill RAC file for the 2-order root
              
              if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length2), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length2), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length2), 9]))^2))}
              else {cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length2), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length2), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length2), 9]))^2))}
              
              rac[r, 1:6]<-c(max(rac[,1], na.rm=TRUE)+1, currentMother, 2, dbase, 1, cumulDist)

              currentMother2<-rac[r,1]
              
              lie.lines<-lie.lines+length2
              
              #----------------------------------------------------------------------------------------------------            
              
              if ("root" %in% names(r2)){
                
                for (r3 in r2){# For each 3-order root
                  
                  if ("geometry" %in% names(r3)){
                    
                    r<-r+1
                    ns <- r3$geometry$polyline
                    length3<-length(ns)
                    
                    age=NULL
                    diameter=NULL
                    if (is.character(final.date)==TRUE & "functions" %in% names(r3)){
                      age<-r3$functions
                      for (i in 1:length(age)){
                        if (final.date %in% age[[i]]$.attrs) {
                          time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
                          if (time[1]<time[2]){time[1]<-time[2]}}}}
                    
                    if ("functions" %in% names(r3)){
                      age<-r3$functions
                      for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                    
                    #c(Num, Date, Bran, Apic, Prec, Suiv)
                    if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length3),1:6]<-c((lie.lines+1):(lie.lines+length3), rep(1, length3), rep(0, length3), rep(0, length3), lie.lines:(lie.lines+length3-1), (lie.lines+2):(lie.lines+length3+1))}
                    if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length3),1:6]<-c((lie.lines+1):(lie.lines+length3), time, rep(0, length3), rep(0, length3), lie.lines:(lie.lines+length3-1), (lie.lines+2):(lie.lines+length3+1))}
                    lie[(lie.lines+1):(lie.lines+length3),11]<-r
                    if(is.null(diameter)==TRUE) {} else {lie[(lie.lines+1):(lie.lines+length3),12]<-diameter}
                    lie[(lie.lines+1):(lie.lines+length3),13]<-3
                    
                    #c(X,Y,Z)
                    lie[(lie.lines+1):(lie.lines+length3),7]<-sapply(ns, xnodes)
                    lie[(lie.lines+1):(lie.lines+length3),8]<-sapply(ns, ynodes)
                    if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length3),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length3),9]<-0}
                    
                    #c(dist)
                    lie[(lie.lines+1):(lie.lines+length3), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length3), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length3), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length3), 9]))^2)))
                    
                    # Search the closest point on the mother root (calculate DBase)
                    
                    parentnode<-which(lie[start2:stop2, 7]==lie[lie.lines+1, 7] & lie[start2:stop2, 8]==lie[lie.lines+1, 8] & lie[start2:stop2, 9]==lie[lie.lines+1, 9])
                    
                    if (length(parentnode)>0){
                      
                      dbase<-lie[start2+parentnode[1]-1, 10]
                      if (connect==TRUE){
                        lie[lie.lines+1,5]<-lie[start2+parentnode[1]-1, 1]
                        dist1<-0}}
                    
                    else { #If no physical connection between parent and daughter root. Interpolate new point or find closest point.
                      
                      alldist<-sqrt((lie[start2:(stop2-1), 7]-lie[lie.lines+1, 7])^2+(lie[start2:(stop2-1), 8]-lie[lie.lines+1, 8])^2+(lie[start2:(stop2-1), 9]-lie[lie.lines+1, 9])^2)
                      
                      scalx<-diff(lie[start2:stop2, 7])*(lie[start2:(stop2-1), 7]-lie[lie.lines+1, 7])
                      scaly<-diff(lie[start2:stop2, 8])*(lie[start2:(stop2-1), 8]-lie[lie.lines+1, 8])
                      scalz<-diff(lie[start2:stop2, 9])*(lie[start2:(stop2-1), 9]-lie[lie.lines+1, 9])
                      d2<-diff(lie[start2:stop2, 7])^2+diff(lie[start2:stop2, 8])^2+diff(lie[start2:stop2, 9])^2
                      t<-(-(scalx+scaly+scalz)/d2)
                      
                      if (length(which(t>=0 & t<=1))==0){
                        
                        index<-which(alldist==min(alldist))
                        dbase<-lie[start2+index-1, 10]
                        if (connect==TRUE){
                          lie[lie.lines+1,5]<-lie[start2+index-1, 1] #Update Prec
                          dist1<-min(alldist)
                          lie[(lie.lines+1):(lie.lines+length3), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length3), 10]}}
                      
                      else {
                        
                        t[t<0]<-NA
                        t[t>1]<-NA
                        xn<-diff(lie[start2:stop2, 7])*t+lie[start2:(stop2-1), 7]
                        yn<-diff(lie[start2:stop2, 8])*t+lie[start2:(stop2-1), 8]
                        zn<-diff(lie[start2:stop2, 9])*t+lie[start2:(stop2-1), 9]
                        dist1<-sqrt((xn-lie[lie.lines+1, 7])^2+(yn-lie[lie.lines+1, 8])^2+(zn-lie[lie.lines+1, 9])^2)
                        
                        if (sum(is.na(dist1)==T)>0) {
                          index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                          dist1<-min(dist1, na.rm=T)} 
                        else {
                          index<-as.numeric(match(min(dist1), dist1))
                          dist1<-min(dist1)}
                        
                        if (dist1>min(alldist)){
                          index<-which(alldist==min(alldist))
                          dist1<-min(alldist)
                          dbase<-lie[start2+index-1, 10]
                          if (connect==TRUE){
                            lie[lie.lines+1,5]<-lie[start2+index-1, 1] #Update Prec
                            dist1<-min(alldist)
                            lie[(lie.lines+1):(lie.lines+length3), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length3), 10]}}
                        
                        else{
                          
                          lie[(lie.lines+1):(lie.lines+length3), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length3), 10]
                          dbase<-lie[start2+index-1, 10]+distance3D(x1=lie[start2+index-1, 7], x2=xn[index], y1=lie[start2+index-1, 8], y2=yn[index], z1=lie[start2+index-1, 9], z2=zn[index])
                          
                          if (connect==TRUE){
                            
                            lie[lie.lines+length3+1,1:13]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie[start2+index-1, 11], NA, 2)
                            lie<-lie[order(lie[, 11], lie[, 10]),]
                            length2<-length2+1
                            lie.lines<-lie.lines+1
                            stop2<-stop2+1
                            pos<-match(NA, lie[1:(lie.lines),1])
                            lie[1:(lie.lines+length3),5]<-match(lie[1:(lie.lines+length3),5], lie[1:(lie.lines+length3),1])
                            lie[1:(lie.lines+length3),6]<-match(lie[1:(lie.lines+length3),6], lie[1:(lie.lines+length3),1])
                            lie[pos, 5]<-pos-1
                            lie[pos, 6]<-pos+1
                            lie[pos, 2]<-lie[lie[pos, 6], 2]
                            lie[pos+1, 5]<-pos
                            lie[pos-1, 6]<-pos
                            lie[which(is.na(lie[1:(lie.lines),5])==TRUE), 5]<-0
                            lie[which(is.na(lie[1:(lie.lines),6])==TRUE), 6]<-0
                            lie[1:(lie.lines+length3),1]<-match(lie[1:(lie.lines+length3),1], lie[1:(lie.lines+length3),1])
                            lie[lie.lines+1,5]<-pos
                            
                            #Calculate diameter for interpolated point (linear function)
                            if (is.null(diameter)==FALSE){
                              prec<-lie[pos, 5]
                              suiv<-lie[pos, 6]
                              slope<-(lie[suiv, 12]-lie[prec, 12])/distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[suiv, 7], y2=lie[suiv, 8], z2=lie[suiv, 9])
                              intercept<-lie[prec, 12]
                              lie[pos, 12]<-intercept+slope*distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[pos, 7], y2=lie[pos, 8], z2=lie[pos, 9])}}}}}
                    
                    if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                    
                    lie[lie.lines+1,3]<-1
                    start3<-as.numeric(lie[lie.lines+1, 1])
                    
                    # Change Suiv and Apic values for the last point of a lateral root
                    lie[lie.lines+length3,c(4,6)]<-c(1, 0)
                    stop3<-as.numeric(lie[lie.lines+length3, 1])
                    
                    # Fill RAC file for the 3-order root
                    
                    if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length3), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length3), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length3), 9]))^2))}
                    else {cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length3), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length3), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length3), 9]))^2))}
                    
                    rac[r, 1:6]<-c(max(rac[,1], na.rm=TRUE)+1, currentMother2, 3, dbase, 1, cumulDist)

                    currentMother3<-rac[r,1]
                    
                    lie.lines<-lie.lines+length3
                    
                    #----------------------------------------------------------------------------------------------------
                    
                    if ("root" %in% names(r3)){
                      
                      for (r4 in r3){# For each 4-order root
                        
                        if ("geometry" %in% names(r4)){
                          
                          r<-r+1
                          ns <- r4$geometry$polyline
                          length4<-length(ns)
                          
                          age=NULL
                          diameter=NULL
                          if (is.character(final.date)==TRUE & "functions" %in% names(r4)){
                            age<-r4$functions
                            for (i in 1:length(age)){
                              if (final.date %in% age[[i]]$.attrs) {
                                time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
                                if (time[1]<time[2]){time[1]<-time[2]}}}}
                          
                          if ("functions" %in% names(r4)){
                            age<-r4$functions
                            for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                          
                          #c(Num, Date, Bran, Apic, Prec, Suiv)
                          if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length4),1:6]<-c((lie.lines+1):(lie.lines+length4), rep(1, length4), rep(0, length4), rep(0, length4), lie.lines:(lie.lines+length4-1), (lie.lines+2):(lie.lines+length4+1))}
                          if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length4),1:6]<-c((lie.lines+1):(lie.lines+length4), time, rep(0, length4), rep(0, length4), lie.lines:(lie.lines+length4-1), (lie.lines+2):(lie.lines+length4+1))}
                          lie[(lie.lines+1):(lie.lines+length4),11]<-r
                          if(is.null(diameter)==TRUE) {} else {lie[(lie.lines+1):(lie.lines+length4),12]<-diameter}
                          lie[(lie.lines+1):(lie.lines+length4),13]<-4
                          
                          #c(X,Y,Z)
                          lie[(lie.lines+1):(lie.lines+length4),7]<-sapply(ns, xnodes)
                          lie[(lie.lines+1):(lie.lines+length4),8]<-sapply(ns, ynodes)
                          if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length4),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length4),9]<-0}
                          
                          #c(dist)
                          lie[(lie.lines+1):(lie.lines+length4), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length4), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length4), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length4), 9]))^2)))
                          
                          # Search the closest point on the mother root (calculate DBase)
                          
                          parentnode<-which(lie[start3:stop3, 7]==lie[lie.lines+1, 7] & lie[start3:stop3, 8]==lie[lie.lines+1, 8] & lie[start3:stop3, 9]==lie[lie.lines+1, 9])
                          
                          if (length(parentnode)>0){
                            
                            dbase<-lie[start3+parentnode[1]-1, 10]
                            if (connect==TRUE){
                              lie[lie.lines+1,5]<-lie[start3+parentnode[1]-1, 1]
                              dist1<-0}}
                          
                          else { #If no physical connection between parent and daughter root. Interpolate new point or find closest point.
                            
                            alldist<-sqrt((lie[start3:(stop3-1), 7]-lie[lie.lines+1, 7])^2+(lie[start3:(stop3-1), 8]-lie[lie.lines+1, 8])^2+(lie[start3:(stop3-1), 9]-lie[lie.lines+1, 9])^2)
                            
                            scalx<-diff(lie[start3:stop3, 7])*(lie[start3:(stop3-1), 7]-lie[lie.lines+1, 7])
                            scaly<-diff(lie[start3:stop3, 8])*(lie[start3:(stop3-1), 8]-lie[lie.lines+1, 8])
                            scalz<-diff(lie[start3:stop3, 9])*(lie[start3:(stop3-1), 9]-lie[lie.lines+1, 9])
                            d2<-diff(lie[start3:stop3, 7])^2+diff(lie[start3:stop3, 8])^2+diff(lie[start3:stop3, 9])^2
                            t<-(-(scalx+scaly+scalz)/d2)
                            
                            if (length(which(t>=0 & t<=1))==0){
                              
                              index<-which(alldist==min(alldist))
                              dbase<-lie[start3+index-1, 10]
                              if (connect==TRUE){
                                lie[lie.lines+1,5]<-lie[start3+index-1, 1] #Update Prec
                                dist1<-min(alldist)
                                lie[(lie.lines+1):(lie.lines+length4), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length4), 10]}}
                            
                            else {
                              
                              t[t<0]<-NA
                              t[t>1]<-NA
                              xn<-diff(lie[start3:stop3, 7])*t+lie[start3:(stop3-1), 7]
                              yn<-diff(lie[start3:stop3, 8])*t+lie[start3:(stop3-1), 8]
                              zn<-diff(lie[start3:stop3, 9])*t+lie[start3:(stop3-1), 9]
                              dist1<-sqrt((xn-lie[lie.lines+1, 7])^2+(yn-lie[lie.lines+1, 8])^2+(zn-lie[lie.lines+1, 9])^2)
                              
                              if (sum(is.na(dist1)==T)>0) {
                                index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                                dist1<-min(dist1, na.rm=T)} 
                              else {
                                index<-as.numeric(match(min(dist1), dist1))
                                dist1<-min(dist1)}
                              
                              if (dist1>min(alldist)){
                                index<-which(alldist==min(alldist))
                                dist1<-min(alldist)
                                dbase<-lie[start3+index-1, 10]
                                if (connect==TRUE){
                                  lie[lie.lines+1,5]<-lie[start3+index-1, 1] #Update Prec
                                  dist1<-min(alldist)
                                  lie[(lie.lines+1):(lie.lines+length4), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length4), 10]}}
                              
                              else{
                                
                                lie[(lie.lines+1):(lie.lines+length4), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length4), 10]
                                dbase<-lie[start3+index-1, 10]+distance3D(x1=lie[start3+index-1, 7], x2=xn[index], y1=lie[start3+index-1, 8], y2=yn[index], z1=lie[start3+index-1, 9], z2=zn[index])
                                
                                if (connect==TRUE){
                                  
                                  lie[lie.lines+length4+1,1:13]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie[start3+index-1, 11], NA, 3)
                                  lie<-lie[order(lie[, 11], lie[, 10]),]
                                  length3<-length3+1
                                  lie.lines<-lie.lines+1
                                  stop3<-stop3+1
                                  pos<-match(NA, lie[1:(lie.lines),1])
                                  lie[1:(lie.lines+length4),5]<-match(lie[1:(lie.lines+length4),5], lie[1:(lie.lines+length4),1])
                                  lie[1:(lie.lines+length4),6]<-match(lie[1:(lie.lines+length4),6], lie[1:(lie.lines+length4),1])
                                  lie[pos, 5]<-pos-1
                                  lie[pos, 6]<-pos+1
                                  lie[pos, 2]<-lie[lie[pos, 6], 2]
                                  lie[pos+1, 5]<-pos
                                  lie[pos-1, 6]<-pos
                                  lie[which(is.na(lie[1:(lie.lines),5])==TRUE), 5]<-0
                                  lie[which(is.na(lie[1:(lie.lines),6])==TRUE), 6]<-0
                                  lie[1:(lie.lines+length4),1]<-match(lie[1:(lie.lines+length4),1], lie[1:(lie.lines+length4),1])
                                  lie[lie.lines+1,5]<-pos
                                  
                                  #Calculate diameter for interpolated point (linear function)
                                  if (is.null(diameter)==FALSE){
                                    prec<-lie[pos, 5]
                                    suiv<-lie[pos, 6]
                                    slope<-(lie[suiv, 12]-lie[prec, 12])/distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[suiv, 7], y2=lie[suiv, 8], z2=lie[suiv, 9])
                                    intercept<-lie[prec, 12]
                                    lie[pos, 12]<-intercept+slope*distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[pos, 7], y2=lie[pos, 8], z2=lie[pos, 9])}}}}}
                          
                          if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                          
                          lie[lie.lines+1,3]<-1
                          start4<-as.numeric(lie[lie.lines+1, 1])
                          
                          # Change Suiv and Apic values for the last point of a lateral root
                          lie[lie.lines+length4,c(4,6)]<-c(1, 0)
                          stop4<-as.numeric(lie[lie.lines+length4, 1])
                          
                          # Fill RAC file for the 4-order root
                          
                          if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length4), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length4), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length4), 9]))^2))}
                          else {cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length4), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length4), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length4), 9]))^2))}
                          
                          rac[r, 1:6]<-c(max(rac[,1], na.rm=TRUE)+1, currentMother3, 4, dbase, 1, cumulDist)

                          currentMother4<-rac[r,1]
                          
                          lie.lines<-lie.lines+length4
                          
                          #----------------------------------------------------------------------------------------------------
                          
                          if ("root" %in% names(r4)){
                            
                            for (r5 in r4){# For each 5-order root
                              
                              if ("geometry" %in% names(r5)){
                                
                                r<-r+1
                                ns <- r5$geometry$polyline
                                length5<-length(ns)
                                
                                age=NULL
                                diameter=NULL
                                if (is.character(final.date)==TRUE & "functions" %in% names(r5)){
                                  age<-r5$functions
                                  for (i in 1:length(age)){
                                    if (final.date %in% age[[i]]$.attrs) {
                                      time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
                                      if (time[1]<time[2]){time[1]<-time[2]}}}}
                                
                                if ("functions" %in% names(r5)){
                                  age<-r5$functions
                                  for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                                
                                #c(Num, Date, Bran, Apic, Prec, Suiv)
                                if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length5),1:6]<-c((lie.lines+1):(lie.lines+length5), rep(1, length5), rep(0, length5), rep(0, length5), lie.lines:(lie.lines+length5-1), (lie.lines+2):(lie.lines+length5+1))}
                                if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length5),1:6]<-c((lie.lines+1):(lie.lines+length5), time, rep(0, length5), rep(0, length5), lie.lines:(lie.lines+length5-1), (lie.lines+2):(lie.lines+length5+1))}
                                lie[(lie.lines+1):(lie.lines+length5),11]<-r
                                if(is.null(diameter)==TRUE) {} else {lie[(lie.lines+1):(lie.lines+length5),12]<-diameter}
                                lie[(lie.lines+1):(lie.lines+length5),13]<-5
                                
                                #c(X,Y,Z)
                                lie[(lie.lines+1):(lie.lines+length5),7]<-sapply(ns, xnodes)
                                lie[(lie.lines+1):(lie.lines+length5),8]<-sapply(ns, ynodes)
                                if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length5),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length5),9]<-0}
                                
                                #c(dist)
                                lie[(lie.lines+1):(lie.lines+length5), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length5), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length5), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length5), 9]))^2)))
                                
                                # Search the closest point on the mother root (calculate DBase)
                                
                                parentnode<-which(lie[start4:stop4, 7]==lie[lie.lines+1, 7] & lie[start4:stop4, 8]==lie[lie.lines+1, 8] & lie[start4:stop4, 9]==lie[lie.lines+1, 9])
                                
                                if (length(parentnode)>0){
                                  
                                  dbase<-lie[start4+parentnode[1]-1, 10]
                                  if (connect==TRUE){
                                    lie[lie.lines+1,5]<-lie[start4+parentnode[1]-1, 1]
                                    dist1<-0}}
                                
                                else { #If no physical connection between parent and daughter root. Interpolate new point or find closest point.
                                  
                                  alldist<-sqrt((lie[start4:(stop4-1), 7]-lie[lie.lines+1, 7])^2+(lie[start4:(stop4-1), 8]-lie[lie.lines+1, 8])^2+(lie[start4:(stop4-1), 9]-lie[lie.lines+1, 9])^2)
                                  
                                  scalx<-diff(lie[start4:stop4, 7])*(lie[start4:(stop4-1), 7]-lie[lie.lines+1, 7])
                                  scaly<-diff(lie[start4:stop4, 8])*(lie[start4:(stop4-1), 8]-lie[lie.lines+1, 8])
                                  scalz<-diff(lie[start4:stop4, 9])*(lie[start4:(stop4-1), 9]-lie[lie.lines+1, 9])
                                  d2<-diff(lie[start4:stop4, 7])^2+diff(lie[start4:stop4, 8])^2+diff(lie[start4:stop4, 9])^2
                                  t<-(-(scalx+scaly+scalz)/d2)
                                  
                                  if (length(which(t>=0 & t<=1))==0){
                                    
                                    index<-which(alldist==min(alldist))
                                    dbase<-lie[start4+index-1, 10]
                                    if (connect==TRUE){
                                      lie[lie.lines+1,5]<-lie[start4+index-1, 1] #Update Prec
                                      dist1<-min(alldist)
                                      lie[(lie.lines+1):(lie.lines+length5), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length5), 10]}}
                                  
                                  else {
                                    
                                    t[t<0]<-NA
                                    t[t>1]<-NA
                                    xn<-diff(lie[start4:stop4, 7])*t+lie[start4:(stop4-1), 7]
                                    yn<-diff(lie[start4:stop4, 8])*t+lie[start4:(stop4-1), 8]
                                    zn<-diff(lie[start4:stop4, 9])*t+lie[start4:(stop4-1), 9]
                                    dist1<-sqrt((xn-lie[lie.lines+1, 7])^2+(yn-lie[lie.lines+1, 8])^2+(zn-lie[lie.lines+1, 9])^2)
                                    
                                    if (sum(is.na(dist1)==T)>0) {
                                      index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                                      dist1<-min(dist1, na.rm=T)} 
                                    else {
                                      index<-as.numeric(match(min(dist1), dist1))
                                      dist1<-min(dist1)}
                                    
                                    if (dist1>min(alldist)){
                                      index<-which(alldist==min(alldist))
                                      dist1<-min(alldist)
                                      dbase<-lie[start4+index-1, 10]
                                      if (connect==TRUE){
                                        lie[lie.lines+1,5]<-lie[start4+index-1, 1] #Update Prec
                                        dist1<-min(alldist)
                                        lie[(lie.lines+1):(lie.lines+length5), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length5), 10]}}
                                    
                                    else{
                                      
                                      lie[(lie.lines+1):(lie.lines+length5), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length5), 10]
                                      dbase<-lie[start4+index-1, 10]+distance3D(x1=lie[start4+index-1, 7], x2=xn[index], y1=lie[start4+index-1, 8], y2=yn[index], z1=lie[start4+index-1, 9], z2=zn[index])
                                      
                                      if (connect==TRUE){
                                        
                                        lie[lie.lines+length5+1,1:13]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie[start4+index-1, 11], NA, 4)
                                        lie<-lie[order(lie[, 11], lie[, 10]),]
                                        length4<-length4+1
                                        lie.lines<-lie.lines+1
                                        stop4<-stop4+1
                                        pos<-match(NA, lie[1:(lie.lines),1])
                                        lie[1:(lie.lines+length5),5]<-match(lie[1:(lie.lines+length5),5], lie[1:(lie.lines+length5),1])
                                        lie[1:(lie.lines+length5),6]<-match(lie[1:(lie.lines+length5),6], lie[1:(lie.lines+length5),1])
                                        lie[pos, 5]<-pos-1
                                        lie[pos, 6]<-pos+1
                                        lie[pos, 2]<-lie[lie[pos, 6], 2]
                                        lie[pos+1, 5]<-pos
                                        lie[pos-1, 6]<-pos
                                        lie[which(is.na(lie[1:(lie.lines),5])==TRUE), 5]<-0
                                        lie[which(is.na(lie[1:(lie.lines),6])==TRUE), 6]<-0
                                        lie[1:(lie.lines+length5),1]<-match(lie[1:(lie.lines+length5),1], lie[1:(lie.lines+length5),1])
                                        lie[lie.lines+1,5]<-pos
                                        
                                        #Calculate diameter for interpolated point (linear function)
                                        if (is.null(diameter)==FALSE){
                                          prec<-lie[pos, 5]
                                          suiv<-lie[pos, 6]
                                          slope<-(lie[suiv, 12]-lie[prec, 12])/distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[suiv, 7], y2=lie[suiv, 8], z2=lie[suiv, 9])
                                          intercept<-lie[prec, 12]
                                          lie[pos, 12]<-intercept+slope*distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[pos, 7], y2=lie[pos, 8], z2=lie[pos, 9])}}}}}
                                
                                if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                                
                                lie[lie.lines+1,3]<-1
                                start5<-as.numeric(lie[lie.lines+1, 1])

                                # Change Suiv and Apic values for the last point of a lateral root
                                lie[lie.lines+length5,c(4,6)]<-c(1, 0)
                                stop5<-as.numeric(lie[lie.lines+length5, 1])

                                # Fill RAC file for the 5-order root
                                
                                if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length5), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length5), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length5), 9]))^2))}
                                else {cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length5), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length5), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length5), 9]))^2))}
                                
                                rac[r, 1:6]<-c(max(rac[,1], na.rm=TRUE)+1, currentMother4, 5, dbase, 1, cumulDist)

                                currentMother5<-rac[r,1]
                                
                                lie.lines<-lie.lines+length5
                                
                                #------------------------------------------------------------------------------------
                                
                                if ("root" %in% names(r5)){
                                  
                                  for (r6 in r5){# For each 6-order root
                                    
                                    if ("geometry" %in% names(r6)){
                                      
                                      r<-r+1
                                      ns <- r6$geometry$polyline
                                      length6<-length(ns)
                                      
                                      age=NULL
                                      diameter=NULL
                                      if (is.character(final.date)==TRUE & "functions" %in% names(r6)){
                                        age<-r6$functions
                                        for (i in 1:length(age)){
                                          if (final.date %in% age[[i]]$.attrs) {
                                            time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
                                            if (time[1]<time[2]){time[1]<-time[2]}}}}
                                      
                                      if ("functions" %in% names(r6)){
                                        age<-r6$functions
                                        for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                                      
                                      #c(Num, Date, Bran, Apic, Prec, Suiv)
                                      if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length6),1:6]<-c((lie.lines+1):(lie.lines+length6), rep(1, length6), rep(0, length6), rep(0, length6), lie.lines:(lie.lines+length6-1), (lie.lines+2):(lie.lines+length6+1))}
                                      if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length6),1:6]<-c((lie.lines+1):(lie.lines+length6), time, rep(0, length6), rep(0, length6), lie.lines:(lie.lines+length6-1), (lie.lines+2):(lie.lines+length6+1))}
                                      lie[(lie.lines+1):(lie.lines+length6),11]<-r
                                      if(is.null(diameter)==TRUE) {} else {lie[(lie.lines+1):(lie.lines+length6),12]<-diameter}
                                      lie[(lie.lines+1):(lie.lines+length6),13]<-6
                                      
                                      #c(X,Y,Z)
                                      lie[(lie.lines+1):(lie.lines+length6),7]<-sapply(ns, xnodes)
                                      lie[(lie.lines+1):(lie.lines+length6),8]<-sapply(ns, ynodes)
                                      if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length6),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length6),9]<-0}
                                      
                                      #c(dist)
                                      lie[(lie.lines+1):(lie.lines+length6), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length6), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length6), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length6), 9]))^2)))
                                      
                                      # Search the closest point on the mother root (calculate DBase)
                                      
                                      parentnode<-which(lie[start5:stop5, 7]==lie[lie.lines+1, 7] & lie[start5:stop5, 8]==lie[lie.lines+1, 8] & lie[start5:stop5, 9]==lie[lie.lines+1, 9])
                                      
                                      if (length(parentnode)>0){
                                        
                                        dbase<-lie[start5+parentnode[1]-1, 10]
                                        if (connect==TRUE){
                                          lie[lie.lines+1,5]<-lie[start5+parentnode[1]-1, 1]
                                          dist1<-0}}
                                      
                                      else { #If no physical connection between parent and daughter root. Interpolate new point or find closest point.
                                        
                                        alldist<-sqrt((lie[start5:(stop5-1), 7]-lie[lie.lines+1, 7])^2+(lie[start5:(stop5-1), 8]-lie[lie.lines+1, 8])^2+(lie[start5:(stop5-1), 9]-lie[lie.lines+1, 9])^2)
                                        
                                        scalx<-diff(lie[start5:stop5, 7])*(lie[start5:(stop5-1), 7]-lie[lie.lines+1, 7])
                                        scaly<-diff(lie[start5:stop5, 8])*(lie[start5:(stop5-1), 8]-lie[lie.lines+1, 8])
                                        scalz<-diff(lie[start5:stop5, 9])*(lie[start5:(stop5-1), 9]-lie[lie.lines+1, 9])
                                        d2<-diff(lie[start5:stop5, 7])^2+diff(lie[start5:stop5, 8])^2+diff(lie[start5:stop5, 9])^2
                                        t<-(-(scalx+scaly+scalz)/d2)
                                        
                                        if (length(which(t>=0 & t<=1))==0){
                                          
                                          index<-which(alldist==min(alldist))
                                          dbase<-lie[start5+index-1, 10]
                                          if (connect==TRUE){
                                            lie[lie.lines+1,5]<-lie[start5+index-1, 1] #Update Prec
                                            dist1<-min(alldist)
                                            lie[(lie.lines+1):(lie.lines+length6), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length6), 10]}}
                                        
                                        else {
                                          
                                          t[t<0]<-NA
                                          t[t>1]<-NA
                                          xn<-diff(lie[start5:stop5, 7])*t+lie[start5:(stop5-1), 7]
                                          yn<-diff(lie[start5:stop5, 8])*t+lie[start5:(stop5-1), 8]
                                          zn<-diff(lie[start5:stop5, 9])*t+lie[start5:(stop5-1), 9]
                                          dist1<-sqrt((xn-lie[lie.lines+1, 7])^2+(yn-lie[lie.lines+1, 8])^2+(zn-lie[lie.lines+1, 9])^2)
                                          
                                          if (sum(is.na(dist1)==T)>0) {
                                            index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                                            dist1<-min(dist1, na.rm=T)} 
                                          else {
                                            index<-as.numeric(match(min(dist1), dist1))
                                            dist1<-min(dist1)}
                                          
                                          if (dist1>min(alldist)){
                                            index<-which(alldist==min(alldist))
                                            dist1<-min(alldist)
                                            dbase<-lie[start5+index-1, 10]
                                            if (connect==TRUE){
                                              lie[lie.lines+1,5]<-lie[start5+index-1, 1] #Update Prec
                                              dist1<-min(alldist)
                                              lie[(lie.lines+1):(lie.lines+length6), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length6), 10]}}
                                          
                                          else{
                                            
                                            lie[(lie.lines+1):(lie.lines+length6), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length6), 10]
                                            dbase<-lie[start5+index-1, 10]+distance3D(x1=lie[start5+index-1, 7], x2=xn[index], y1=lie[start5+index-1, 8], y2=yn[index], z1=lie[start5+index-1, 9], z2=zn[index])
                                            
                                            if (connect==TRUE){
                                              
                                              lie[lie.lines+length6+1,1:13]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie[start5+index-1, 11], NA, 5)
                                              lie<-lie[order(lie[, 11], lie[, 10]),]
                                              length5<-length5+1
                                              lie.lines<-lie.lines+1
                                              stop5<-stop5+1
                                              pos<-match(NA, lie[1:(lie.lines),1])
                                              lie[1:(lie.lines+length6),5]<-match(lie[1:(lie.lines+length6),5], lie[1:(lie.lines+length6),1])
                                              lie[1:(lie.lines+length6),6]<-match(lie[1:(lie.lines+length6),6], lie[1:(lie.lines+length6),1])
                                              lie[pos, 5]<-pos-1
                                              lie[pos, 6]<-pos+1
                                              lie[pos, 2]<-lie[lie[pos, 6], 2]
                                              lie[pos+1, 5]<-pos
                                              lie[pos-1, 6]<-pos
                                              lie[which(is.na(lie[1:(lie.lines),5])==TRUE), 5]<-0
                                              lie[which(is.na(lie[1:(lie.lines),6])==TRUE), 6]<-0
                                              lie[1:(lie.lines+length6),1]<-match(lie[1:(lie.lines+length6),1], lie[1:(lie.lines+length6),1])
                                              lie[lie.lines+1,5]<-pos
                                              
                                              #Calculate diameter for interpolated point (linear function)
                                              if (is.null(diameter)==FALSE){
                                                prec<-lie[pos, 5]
                                                suiv<-lie[pos, 6]
                                                slope<-(lie[suiv, 12]-lie[prec, 12])/distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[suiv, 7], y2=lie[suiv, 8], z2=lie[suiv, 9])
                                                intercept<-lie[prec, 12]
                                                lie[pos, 12]<-intercept+slope*distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[pos, 7], y2=lie[pos, 8], z2=lie[pos, 9])}}}}}
                                      
                                      if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                                      
                                      lie[lie.lines+1,3]<-1
                                      start6<-as.numeric(lie[lie.lines+1, 1])
                                      
                                      # Change Suiv and Apic values for the last point of a lateral root
                                      lie[lie.lines+length6,c(4,6)]<-c(1, 0)
                                      stop6<-as.numeric(lie[lie.lines+length6, 1])
                                      
                                      # Fill RAC file for the 6-order root
                                      
                                      if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length6), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length6), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length6), 9]))^2))}
                                      else {cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length6), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length6), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length6), 9]))^2))}
                                      
                                      rac[r, 1:6]<-c(max(rac[,1], na.rm=TRUE)+1, currentMother5, 6, dbase, 1, cumulDist)

                                      lie.lines<-lie.lines+length6
                                      
                                      #----------------------------------------------
                                      
                                      if ("root" %in% names(r6)){
                                        
                                        for (r7 in r6){# For each 7-order root
                                          
                                          if ("geometry" %in% names(r7)){
                                            
                                            r<-r+1
                                            ns <- r7$geometry$polyline
                                            length7<-length(ns)
                                            
                                            age=NULL
                                            diameter=NULL
                                            if (is.character(final.date)==TRUE & "functions" %in% names(r7)){
                                              age<-r7$functions
                                              for (i in 1:length(age)){
                                                if (final.date %in% age[[i]]$.attrs) {
                                                  time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
                                                  if (time[1]<time[2]){time[1]<-time[2]}}}}
                                            
                                            if ("functions" %in% names(r7)){
                                              age<-r7$functions
                                              for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                                            
                                            #c(Num, Date, Bran, Apic, Prec, Suiv)
                                            if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length7),1:6]<-c((lie.lines+1):(lie.lines+length7), rep(1, length7), rep(0, length7), rep(0, length7), lie.lines:(lie.lines+length7-1), (lie.lines+2):(lie.lines+length7+1))}
                                            if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length7),1:6]<-c((lie.lines+1):(lie.lines+length7), time, rep(0, length7), rep(0, length7), lie.lines:(lie.lines+length7-1), (lie.lines+2):(lie.lines+length7+1))}
                                            lie[(lie.lines+1):(lie.lines+length7),11]<-r
                                            if(is.null(diameter)==TRUE) {} else {lie[(lie.lines+1):(lie.lines+length7),12]<-diameter}
                                            lie[(lie.lines+1):(lie.lines+length7),13]<-7 #order
                                            
                                            #c(X,Y,Z)
                                            lie[(lie.lines+1):(lie.lines+length7),7]<-sapply(ns, xnodes)
                                            lie[(lie.lines+1):(lie.lines+length7),8]<-sapply(ns, ynodes)
                                            if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length7),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length7),9]<-0}
                                            
                                            #c(dist)
                                            lie[(lie.lines+1):(lie.lines+length7), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length7), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length7), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length7), 9]))^2)))
                                            
                                            # Search the closest point on the mother root (calculate DBase)
                                            
                                            parentnode<-which(lie[start6:stop6, 7]==lie[lie.lines+1, 7] & lie[start6:stop6, 8]==lie[lie.lines+1, 8] & lie[start6:stop6, 9]==lie[lie.lines+1, 9])
                                            
                                            if (length(parentnode)>0){
                                              
                                              dbase<-lie[start6+parentnode[1]-1, 10]
                                              if (connect==TRUE){
                                                lie[lie.lines+1,5]<-lie[start6+parentnode[1]-1, 1]
                                                dist1<-0}}
                                            
                                            else { #If no physical connection between parent and daughter root. Interpolate new point or find closest point.
                                              
                                              alldist<-sqrt((lie[start6:(stop6-1), 7]-lie[lie.lines+1, 7])^2+(lie[start6:(stop6-1), 8]-lie[lie.lines+1, 8])^2+(lie[start6:(stop6-1), 9]-lie[lie.lines+1, 9])^2)
                                              
                                              scalx<-diff(lie[start6:stop6, 7])*(lie[start6:(stop6-1), 7]-lie[lie.lines+1, 7])
                                              scaly<-diff(lie[start6:stop6, 8])*(lie[start6:(stop6-1), 8]-lie[lie.lines+1, 8])
                                              scalz<-diff(lie[start6:stop6, 9])*(lie[start6:(stop6-1), 9]-lie[lie.lines+1, 9])
                                              d2<-diff(lie[start6:stop6, 7])^2+diff(lie[start6:stop6, 8])^2+diff(lie[start6:stop6, 9])^2
                                              t<-(-(scalx+scaly+scalz)/d2)
                                              
                                              if (length(which(t>=0 & t<=1))==0){
                                                
                                                index<-which(alldist==min(alldist))
                                                dbase<-lie[start6+index-1, 10]
                                                if (connect==TRUE){
                                                  lie[lie.lines+1,5]<-lie[start6+index-1, 1] #Update Prec
                                                  dist1<-min(alldist)
                                                  lie[(lie.lines+1):(lie.lines+length7), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length7), 10]}}
                                              
                                              else {
                                                
                                                t[t<0]<-NA
                                                t[t>1]<-NA
                                                xn<-diff(lie[start6:stop6, 7])*t+lie[start6:(stop6-1), 7]
                                                yn<-diff(lie[start6:stop6, 8])*t+lie[start6:(stop6-1), 8]
                                                zn<-diff(lie[start6:stop6, 9])*t+lie[start6:(stop6-1), 9]
                                                dist1<-sqrt((xn-lie[lie.lines+1, 7])^2+(yn-lie[lie.lines+1, 8])^2+(zn-lie[lie.lines+1, 9])^2)
                                                
                                                if (sum(is.na(dist1)==T)>0) {
                                                  index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                                                  dist1<-min(dist1, na.rm=T)} 
                                                else {
                                                  index<-as.numeric(match(min(dist1), dist1))
                                                  dist1<-min(dist1)}
                                                
                                                if (dist1>min(alldist)){
                                                  index<-which(alldist==min(alldist))
                                                  dist1<-min(alldist)
                                                  dbase<-lie[start6+index-1, 10]
                                                  if (connect==TRUE){
                                                    lie[lie.lines+1,5]<-lie[start6+index-1, 1] #Update Prec
                                                    dist1<-min(alldist)
                                                    lie[(lie.lines+1):(lie.lines+length7), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length7), 10]}}
                                                
                                                else{
                                                  
                                                  lie[(lie.lines+1):(lie.lines+length7), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length7), 10]
                                                  dbase<-lie[start6+index-1, 10]+distance3D(x1=lie[start6+index-1, 7], x2=xn[index], y1=lie[start6+index-1, 8], y2=yn[index], z1=lie[start6+index-1, 9], z2=zn[index])
                                                  
                                                  if (connect==TRUE){
                                                    
                                                    lie[lie.lines+length7+1,1:13]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie[start6+index-1, 11], NA, 6)
                                                    lie<-lie[order(lie[, 11], lie[, 10]),]
                                                    length6<-length6+1
                                                    lie.lines<-lie.lines+1
                                                    stop6<-stop6+1
                                                    pos<-match(NA, lie[1:(lie.lines),1])
                                                    lie[1:(lie.lines+length7),5]<-match(lie[1:(lie.lines+length7),5], lie[1:(lie.lines+length7),1])
                                                    lie[1:(lie.lines+length7),6]<-match(lie[1:(lie.lines+length7),6], lie[1:(lie.lines+length7),1])
                                                    lie[pos, 5]<-pos-1
                                                    lie[pos, 6]<-pos+1
                                                    lie[pos, 2]<-lie[lie[pos, 6], 2]
                                                    lie[pos+1, 5]<-pos
                                                    lie[pos-1, 6]<-pos
                                                    lie[which(is.na(lie[1:(lie.lines),5])==TRUE), 5]<-0
                                                    lie[which(is.na(lie[1:(lie.lines),6])==TRUE), 6]<-0
                                                    lie[1:(lie.lines+length7),1]<-match(lie[1:(lie.lines+length7),1], lie[1:(lie.lines+length7),1])
                                                    lie[lie.lines+1,5]<-pos
                                                    
                                                    #Calculate diameter for interpolated point (linear function)
                                                    if (is.null(diameter)==FALSE){
                                                      prec<-lie[pos, 5]
                                                      suiv<-lie[pos, 6]
                                                      slope<-(lie[suiv, 12]-lie[prec, 12])/distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[suiv, 7], y2=lie[suiv, 8], z2=lie[suiv, 9])
                                                      intercept<-lie[prec, 12]
                                                      lie[pos, 12]<-intercept+slope*distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[pos, 7], y2=lie[pos, 8], z2=lie[pos, 9])}}}}}
                                            
                                            if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                                            
                                            lie[lie.lines+1,3]<-1
                                            start7<-as.numeric(lie[lie.lines+1, 1])
                                            
                                            # Change Suiv and Apic values for the last point of a lateral root
                                            lie[lie.lines+length7,c(4,6)]<-c(1, 0)
                                            stop7<-as.numeric(lie[lie.lines+length7, 1])
                                            
                                            # Fill RAC file for the 7-order root
                                            
                                            if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length7), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length7), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length7), 9]))^2))}
                                            else {cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length7), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length7), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length7), 9]))^2))}
                                            
                                            rac[r, 1:6]<-c(max(rac[,1], na.rm=TRUE)+1, currentMother5, 7, dbase, 1, cumulDist)
                                            
                                            lie.lines<-lie.lines+length7
                                            
                                            #-------------------------------------------
                                            
                                            # if ("root" %in% names(r7)){
                                            #   
                                            #   for (r8 in r7){# For each 8-order root
                                            #     
                                            #     if ("geometry" %in% names(r8)){
                                            #       
                                            #       r<-r+1
                                            #       ns <- r8$geometry$polyline
                                            #       length8<-length(ns)
                                            #       
                                            #       age=NULL
                                            #       diameter=NULL
                                            #       if (is.character(final.date)==TRUE & "functions" %in% names(r8)){
                                            #         age<-r8$functions
                                            #         for (i in 1:length(age)){
                                            #           if (final.date %in% age[[i]]$.attrs) {
                                            #             time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
                                            #             if (time[1]<time[2]){time[1]<-time[2]}}}}
                                            #       
                                            #       if ("functions" %in% names(r8)){
                                            #         age<-r8$functions
                                            #         for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                                            #       
                                            #       #c(Num, Date, Bran, Apic, Prec, Suiv)
                                            #       if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length8),1:6]<-c((lie.lines+1):(lie.lines+length8), rep(1, length8), rep(0, length8), rep(0, length8), lie.lines:(lie.lines+length8-1), (lie.lines+2):(lie.lines+length8+1))}
                                            #       if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length8),1:6]<-c((lie.lines+1):(lie.lines+length8), time, rep(0, length8), rep(0, length8), lie.lines:(lie.lines+length8-1), (lie.lines+2):(lie.lines+length8+1))}
                                            #       lie[(lie.lines+1):(lie.lines+length8),11]<-r
                                            #       if(is.null(diameter)==TRUE) {} else {lie[(lie.lines+1):(lie.lines+length8),12]<-diameter}
                                            #       lie[(lie.lines+1):(lie.lines+length8),13]<-8 #order
                                            #       
                                            #       #c(X,Y,Z)
                                            #       lie[(lie.lines+1):(lie.lines+length8),7]<-sapply(ns, xnodes)
                                            #       lie[(lie.lines+1):(lie.lines+length8),8]<-sapply(ns, ynodes)
                                            #       if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length8),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length8),9]<-0}
                                            #       
                                            #       #c(dist)
                                            #       lie[(lie.lines+1):(lie.lines+length8), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length8), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length8), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length8), 9]))^2)))
                                            #       
                                            #       # Search the closest point on the mother root (calculate DBase)
                                            #       
                                            #       parentnode<-which(lie[start7:stop7, 7]==lie[lie.lines+1, 7] & lie[start7:stop7, 8]==lie[lie.lines+1, 8] & lie[start7:stop7, 9]==lie[lie.lines+1, 9])
                                            #       
                                            #       if (length(parentnode)>0){
                                            #         
                                            #         dbase<-lie[start7+parentnode[1]-1, 10]
                                            #         if (connect==TRUE){
                                            #           lie[lie.lines+1,5]<-lie[start7+parentnode[1]-1, 1]
                                            #           dist1<-0}}
                                            #       
                                            #       else { #If no physical connection between parent and daughter root. Interpolate new point or find closest point.
                                            #         
                                            #         alldist<-sqrt((lie[start7:(stop7-1), 7]-lie[lie.lines+1, 7])^2+(lie[start7:(stop7-1), 8]-lie[lie.lines+1, 8])^2+(lie[start7:(stop7-1), 9]-lie[lie.lines+1, 9])^2)
                                            #         
                                            #         scalx<-diff(lie[start7:stop7, 7])*(lie[start7:(stop7-1), 7]-lie[lie.lines+1, 7])
                                            #         scaly<-diff(lie[start7:stop7, 8])*(lie[start7:(stop7-1), 8]-lie[lie.lines+1, 8])
                                            #         scalz<-diff(lie[start7:stop7, 9])*(lie[start7:(stop7-1), 9]-lie[lie.lines+1, 9])
                                            #         d2<-diff(lie[start7:stop7, 7])^2+diff(lie[start7:stop7, 8])^2+diff(lie[start7:stop7, 9])^2
                                            #         t<-(-(scalx+scaly+scalz)/d2)
                                            #         
                                            #         if (length(which(t>=0 & t<=1))==0){
                                            #           
                                            #           index<-which(alldist==min(alldist))
                                            #           dbase<-lie[start7+index-1, 10]
                                            #           if (connect==TRUE){
                                            #             lie[lie.lines+1,5]<-lie[start7+index-1, 1] #Update Prec
                                            #             dist1<-min(alldist)
                                            #             lie[(lie.lines+1):(lie.lines+length8), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length8), 10]}}
                                            #         
                                            #         else {
                                            #           
                                            #           t[t<0]<-NA
                                            #           t[t>1]<-NA
                                            #           xn<-diff(lie[start7:stop7, 7])*t+lie[start7:(stop7-1), 7]
                                            #           yn<-diff(lie[start7:stop7, 8])*t+lie[start7:(stop7-1), 8]
                                            #           zn<-diff(lie[start7:stop7, 9])*t+lie[start7:(stop7-1), 9]
                                            #           dist1<-sqrt((xn-lie[lie.lines+1, 7])^2+(yn-lie[lie.lines+1, 8])^2+(zn-lie[lie.lines+1, 9])^2)
                                            #           
                                            #           if (sum(is.na(dist1)==T)>0) {
                                            #             index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                                            #             dist1<-min(dist1, na.rm=T)} 
                                            #           else {
                                            #             index<-as.numeric(match(min(dist1), dist1))
                                            #             dist1<-min(dist1)}
                                            #           
                                            #           if (dist1>min(alldist)){
                                            #             index<-which(alldist==min(alldist))
                                            #             dist1<-min(alldist)
                                            #             dbase<-lie[start7+index-1, 10]
                                            #             if (connect==TRUE){
                                            #               lie[lie.lines+1,5]<-lie[start7+index-1, 1] #Update Prec
                                            #               dist1<-min(alldist)
                                            #               lie[(lie.lines+1):(lie.lines+length8), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length8), 10]}}
                                            #           
                                            #           else{
                                            #             
                                            #             lie[(lie.lines+1):(lie.lines+length8), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length8), 10]
                                            #             dbase<-lie[start7+index-1, 10]+distance3D(x1=lie[start7+index-1, 7], x2=xn[index], y1=lie[start7+index-1, 8], y2=yn[index], z1=lie[start7+index-1, 9], z2=zn[index])
                                            #             
                                            #             if (connect==TRUE){
                                            #               
                                            #               lie[lie.lines+length8+1,1:13]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie[start7+index-1, 11], NA, 7)
                                            #               lie<-lie[order(lie[, 11], lie[, 10]),]
                                            #               length7<-length7+1
                                            #               lie.lines<-lie.lines+1
                                            #               stop7<-stop7+1
                                            #               pos<-match(NA, lie[1:(lie.lines),1])
                                            #               lie[1:(lie.lines+length8),5]<-match(lie[1:(lie.lines+length8),5], lie[1:(lie.lines+length8),1])
                                            #               lie[1:(lie.lines+length8),6]<-match(lie[1:(lie.lines+length8),6], lie[1:(lie.lines+length8),1])
                                            #               lie[pos, 5]<-pos-1
                                            #               lie[pos, 6]<-pos+1
                                            #               lie[pos, 2]<-lie[lie[pos, 6], 2]
                                            #               lie[pos+1, 5]<-pos
                                            #               lie[pos-1, 6]<-pos
                                            #               lie[which(is.na(lie[1:(lie.lines),5])==TRUE), 5]<-0
                                            #               lie[which(is.na(lie[1:(lie.lines),6])==TRUE), 6]<-0
                                            #               lie[1:(lie.lines+length8),1]<-match(lie[1:(lie.lines+length8),1], lie[1:(lie.lines+length8),1])
                                            #               lie[lie.lines+1,5]<-pos
                                            #               
                                            #               #Calculate diameter for interpolated point (linear function)
                                            #               if (is.null(diameter)==FALSE){
                                            #                 prec<-lie[pos, 5]
                                            #                 suiv<-lie[pos, 6]
                                            #                 slope<-(lie[suiv, 12]-lie[prec, 12])/distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[suiv, 7], y2=lie[suiv, 8], z2=lie[suiv, 9])
                                            #                 intercept<-lie[prec, 12]
                                            #                 lie[pos, 12]<-intercept+slope*distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[pos, 7], y2=lie[pos, 8], z2=lie[pos, 9])}}}}}
                                            #       
                                            #       if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                                            #       
                                            #       lie[lie.lines+1,3]<-1
                                            #       start8<-as.numeric(lie[lie.lines+1, 1])
                                            #       
                                            #       # Change Suiv and Apic values for the last point of a lateral root
                                            #       lie[lie.lines+length8,c(4,6)]<-c(1, 0)
                                            #       stop8<-as.numeric(lie[lie.lines+length8, 1])
                                            #       
                                            #       # Fill RAC file for the 8-order root
                                            #       
                                            #       if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length8), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length8), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length8), 9]))^2))}
                                            #       else {cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length8), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length8), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length8), 9]))^2))}
                                            #       
                                            #       rac[r, 1:6]<-c(max(rac[,1], na.rm=TRUE)+1, currentMother5, 8, dbase, 1, cumulDist)
                                            #       
                                            #       lie.lines<-lie.lines+length8
                                            #       
                                            #       #---------------------------------------------
                                            #       
                                            #       # if ("root" %in% names(r8)){
                                            #       #   
                                            #       #   for (r9 in r8){# For each 9-order root
                                            #       #     
                                            #       #     if ("geometry" %in% names(r9)){
                                            #       #       
                                            #       #       r<-r+1
                                            #       #       ns <- r9$geometry$polyline
                                            #       #       length9<-length(ns)
                                            #       #       
                                            #       #       age=NULL
                                            #       #       diameter=NULL
                                            #       #       if (is.character(final.date)==TRUE & "functions" %in% names(r9)){
                                            #       #         age<-r9$functions
                                            #       #         for (i in 1:length(age)){
                                            #       #           if (final.date %in% age[[i]]$.attrs) {
                                            #       #             time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
                                            #       #             if (time[1]<time[2]){time[1]<-time[2]}}}}
                                            #       #       
                                            #       #       if ("functions" %in% names(r9)){
                                            #       #         age<-r9$functions
                                            #       #         for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                                            #       #       
                                            #       #       #c(Num, Date, Bran, Apic, Prec, Suiv)
                                            #       #       if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length9),1:6]<-c((lie.lines+1):(lie.lines+length9), rep(1, length9), rep(0, length9), rep(0, length9), lie.lines:(lie.lines+length9-1), (lie.lines+2):(lie.lines+length9+1))}
                                            #       #       if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length9),1:6]<-c((lie.lines+1):(lie.lines+length9), time, rep(0, length9), rep(0, length9), lie.lines:(lie.lines+length9-1), (lie.lines+2):(lie.lines+length9+1))}
                                            #       #       lie[(lie.lines+1):(lie.lines+length9),11]<-r
                                            #       #       if(is.null(diameter)==TRUE) {} else {lie[(lie.lines+1):(lie.lines+length9),12]<-diameter}
                                            #       #       lie[(lie.lines+1):(lie.lines+length9),13]<-9 #order
                                            #       #       
                                            #       #       #c(X,Y,Z)
                                            #       #       lie[(lie.lines+1):(lie.lines+length9),7]<-sapply(ns, xnodes)
                                            #       #       lie[(lie.lines+1):(lie.lines+length9),8]<-sapply(ns, ynodes)
                                            #       #       if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length9),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length9),9]<-0}
                                            #       #       
                                            #       #       #c(dist)
                                            #       #       lie[(lie.lines+1):(lie.lines+length9), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length9), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length9), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length9), 9]))^2)))
                                            #       #       
                                            #       #       # Search the closest point on the mother root (calculate DBase)
                                            #       #       
                                            #       #       parentnode<-which(lie[start8:stop8, 7]==lie[lie.lines+1, 7] & lie[start8:stop8, 8]==lie[lie.lines+1, 8] & lie[start8:stop8, 9]==lie[lie.lines+1, 9])
                                            #       #       
                                            #       #       if (length(parentnode)>0){
                                            #       #         
                                            #       #         dbase<-lie[start8+parentnode[1]-1, 10]
                                            #       #         if (connect==TRUE){
                                            #       #           lie[lie.lines+1,5]<-lie[start8+parentnode[1]-1, 1]
                                            #       #           dist1<-0}}
                                            #       #       
                                            #       #       else { #If no physical connection between parent and daughter root. Interpolate new point or find closest point.
                                            #       #         
                                            #       #         alldist<-sqrt((lie[start8:(stop8-1), 7]-lie[lie.lines+1, 7])^2+(lie[start8:(stop8-1), 8]-lie[lie.lines+1, 8])^2+(lie[start8:(stop8-1), 9]-lie[lie.lines+1, 9])^2)
                                            #       #         
                                            #       #         scalx<-diff(lie[start8:stop8, 7])*(lie[start8:(stop8-1), 7]-lie[lie.lines+1, 7])
                                            #       #         scaly<-diff(lie[start8:stop8, 8])*(lie[start8:(stop8-1), 8]-lie[lie.lines+1, 8])
                                            #       #         scalz<-diff(lie[start8:stop8, 9])*(lie[start8:(stop8-1), 9]-lie[lie.lines+1, 9])
                                            #       #         d2<-diff(lie[start8:stop8, 7])^2+diff(lie[start8:stop8, 8])^2+diff(lie[start8:stop8, 9])^2
                                            #       #         t<-(-(scalx+scaly+scalz)/d2)
                                            #       #         
                                            #       #         if (length(which(t>=0 & t<=1))==0){
                                            #       #           
                                            #       #           index<-which(alldist==min(alldist))
                                            #       #           dbase<-lie[start8+index-1, 10]
                                            #       #           if (connect==TRUE){
                                            #       #             lie[lie.lines+1,5]<-lie[start8+index-1, 1] #Update Prec
                                            #       #             dist1<-min(alldist)
                                            #       #             lie[(lie.lines+1):(lie.lines+length9), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length9), 10]}}
                                            #       #         
                                            #       #         else {
                                            #       #           
                                            #       #           t[t<0]<-NA
                                            #       #           t[t>1]<-NA
                                            #       #           xn<-diff(lie[start8:stop8, 7])*t+lie[start8:(stop8-1), 7]
                                            #       #           yn<-diff(lie[start8:stop8, 8])*t+lie[start8:(stop8-1), 8]
                                            #       #           zn<-diff(lie[start8:stop8, 9])*t+lie[start8:(stop8-1), 9]
                                            #       #           dist1<-sqrt((xn-lie[lie.lines+1, 7])^2+(yn-lie[lie.lines+1, 8])^2+(zn-lie[lie.lines+1, 9])^2)
                                            #       #           
                                            #       #           if (sum(is.na(dist1)==T)>0) {
                                            #       #             index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                                            #       #             dist1<-min(dist1, na.rm=T)} 
                                            #       #           else {
                                            #       #             index<-as.numeric(match(min(dist1), dist1))
                                            #       #             dist1<-min(dist1)}
                                            #       #           
                                            #       #           if (dist1>min(alldist)){
                                            #       #             index<-which(alldist==min(alldist))
                                            #       #             dist1<-min(alldist)
                                            #       #             dbase<-lie[start8+index-1, 10]
                                            #       #             if (connect==TRUE){
                                            #       #               lie[lie.lines+1,5]<-lie[start8+index-1, 1] #Update Prec
                                            #       #               dist1<-min(alldist)
                                            #       #               lie[(lie.lines+1):(lie.lines+length9), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length9), 10]}}
                                            #       #           
                                            #       #           else{
                                            #       #             
                                            #       #             lie[(lie.lines+1):(lie.lines+length9), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length9), 10]
                                            #       #             dbase<-lie[start8+index-1, 10]+distance3D(x1=lie[start8+index-1, 7], x2=xn[index], y1=lie[start8+index-1, 8], y2=yn[index], z1=lie[start8+index-1, 9], z2=zn[index])
                                            #       #             
                                            #       #             if (connect==TRUE){
                                            #       #               
                                            #       #               lie[lie.lines+length9+1,1:13]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie[start8+index-1, 11], NA, 8)
                                            #       #               lie<-lie[order(lie[, 11], lie[, 10]),]
                                            #       #               length8<-length8+1
                                            #       #               lie.lines<-lie.lines+1
                                            #       #               stop8<-stop8+1
                                            #       #               pos<-match(NA, lie[1:(lie.lines),1])
                                            #       #               lie[1:(lie.lines+length9),5]<-match(lie[1:(lie.lines+length9),5], lie[1:(lie.lines+length9),1])
                                            #       #               lie[1:(lie.lines+length9),6]<-match(lie[1:(lie.lines+length9),6], lie[1:(lie.lines+length9),1])
                                            #       #               lie[pos, 5]<-pos-1
                                            #       #               lie[pos, 6]<-pos+1
                                            #       #               lie[pos, 2]<-lie[lie[pos, 6], 2]
                                            #       #               lie[pos+1, 5]<-pos
                                            #       #               lie[pos-1, 6]<-pos
                                            #       #               lie[which(is.na(lie[1:(lie.lines),5])==TRUE), 5]<-0
                                            #       #               lie[which(is.na(lie[1:(lie.lines),6])==TRUE), 6]<-0
                                            #       #               lie[1:(lie.lines+length9),1]<-match(lie[1:(lie.lines+length9),1], lie[1:(lie.lines+length9),1])
                                            #       #               lie[lie.lines+1,5]<-pos
                                            #       #               
                                            #       #               #Calculate diameter for interpolated point (linear function)
                                            #       #               if (is.null(diameter)==FALSE){
                                            #       #                 prec<-lie[pos, 5]
                                            #       #                 suiv<-lie[pos, 6]
                                            #       #                 slope<-(lie[suiv, 12]-lie[prec, 12])/distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[suiv, 7], y2=lie[suiv, 8], z2=lie[suiv, 9])
                                            #       #                 intercept<-lie[prec, 12]
                                            #       #                 lie[pos, 12]<-intercept+slope*distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[pos, 7], y2=lie[pos, 8], z2=lie[pos, 9])}}}}}
                                            #       #       
                                            #       #       if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                                            #       #       
                                            #       #       lie[lie.lines+1,3]<-1
                                            #       #       start9<-as.numeric(lie[lie.lines+1, 1])
                                            #       #       
                                            #       #       # Change Suiv and Apic values for the last point of a lateral root
                                            #       #       lie[lie.lines+length9,c(4,6)]<-c(1, 0)
                                            #       #       stop9<-as.numeric(lie[lie.lines+length9, 1])
                                            #       #       
                                            #       #       # Fill RAC file for the 9-order root
                                            #       #       
                                            #       #       if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length9), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length9), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length9), 9]))^2))}
                                            #       #       else {cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length9), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length9), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length9), 9]))^2))}
                                            #       #       
                                            #       #       rac[r, 1:6]<-c(max(rac[,1], na.rm=TRUE)+1, currentMother5, 9, dbase, 1, cumulDist)
                                            #       #       
                                            #       #       lie.lines<-lie.lines+length9
                                            #       #       
                                            #       #       #---------------------------------------------
                                            #       #       
                                            #       #       if ("root" %in% names(r9)){
                                            #       #         
                                            #       #         for (r10 in r9){# For each 10-order root
                                            #       #           
                                            #       #           if ("geometry" %in% names(r10)){
                                            #       #             
                                            #       #             r<-r+1
                                            #       #             ns <- r10$geometry$polyline
                                            #       #             length10<-length(ns)
                                            #       #             
                                            #       #             age=NULL
                                            #       #             diameter=NULL
                                            #       #             if (is.character(final.date)==TRUE & "functions" %in% names(r10)){
                                            #       #               age<-r10$functions
                                            #       #               for (i in 1:length(age)){
                                            #       #                 if (final.date %in% age[[i]]$.attrs) {
                                            #       #                   time<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)
                                            #       #                   if (time[1]<time[2]){time[1]<-time[2]}}}}
                                            #       #             
                                            #       #             if ("functions" %in% names(r10)){
                                            #       #               age<-r10$functions
                                            #       #               for (i in 1:length(age)){if ("diameter" %in% age[[i]]$.attrs) {diameter<-sapply(age[[i]][1:(length(age[[i]])-1)], xnodes)}}}
                                            #       #             
                                            #       #             #c(Num, Date, Bran, Apic, Prec, Suiv)
                                            #       #             if (timeserie==FALSE) {lie[(lie.lines+1):(lie.lines+length10),1:6]<-c((lie.lines+1):(lie.lines+length10), rep(1, length10), rep(0, length10), rep(0, length10), lie.lines:(lie.lines+length10-1), (lie.lines+2):(lie.lines+length10+1))}
                                            #       #             if (timeserie==TRUE) {lie[(lie.lines+1):(lie.lines+length10),1:6]<-c((lie.lines+1):(lie.lines+length10), time, rep(0, length10), rep(0, length10), lie.lines:(lie.lines+length10-1), (lie.lines+2):(lie.lines+length10+1))}
                                            #       #             lie[(lie.lines+1):(lie.lines+length10),11]<-r
                                            #       #             if(is.null(diameter)==TRUE) {} else {lie[(lie.lines+1):(lie.lines+length10),12]<-diameter}
                                            #       #             lie[(lie.lines+1):(lie.lines+length10),13]<-10 #order
                                            #       #             
                                            #       #             #c(X,Y,Z)
                                            #       #             lie[(lie.lines+1):(lie.lines+length10),7]<-sapply(ns, xnodes)
                                            #       #             lie[(lie.lines+1):(lie.lines+length10),8]<-sapply(ns, ynodes)
                                            #       #             if (length(ns[[1]])==3) {lie[(lie.lines+1):(lie.lines+length10),9]<-sapply(ns, znodes)} else {lie[(lie.lines+1):(lie.lines+length10),9]<-0}
                                            #       #             
                                            #       #             #c(dist)
                                            #       #             lie[(lie.lines+1):(lie.lines+length10), 10]<-c(0, cumsum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length10), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length10), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length10), 9]))^2)))
                                            #       #             
                                            #       #             # Search the closest point on the mother root (calculate DBase)
                                            #       #             
                                            #       #             parentnode<-which(lie[start9:stop9, 7]==lie[lie.lines+1, 7] & lie[start9:stop9, 8]==lie[lie.lines+1, 8] & lie[start9:stop9, 9]==lie[lie.lines+1, 9])
                                            #       #             
                                            #       #             if (length(parentnode)>0){
                                            #       #               
                                            #       #               dbase<-lie[start9+parentnode[1]-1, 10]
                                            #       #               if (connect==TRUE){
                                            #       #                 lie[lie.lines+1,5]<-lie[start9+parentnode[1]-1, 1]
                                            #       #                 dist1<-0}}
                                            #       #             
                                            #       #             else { #If no physical connection between parent and daughter root. Interpolate new point or find closest point.
                                            #       #               
                                            #       #               alldist<-sqrt((lie[start9:(stop9-1), 7]-lie[lie.lines+1, 7])^2+(lie[start9:(stop9-1), 8]-lie[lie.lines+1, 8])^2+(lie[start9:(stop9-1), 9]-lie[lie.lines+1, 9])^2)
                                            #       #               
                                            #       #               scalx<-diff(lie[start9:stop9, 7])*(lie[start9:(stop9-1), 7]-lie[lie.lines+1, 7])
                                            #       #               scaly<-diff(lie[start9:stop9, 8])*(lie[start9:(stop9-1), 8]-lie[lie.lines+1, 8])
                                            #       #               scalz<-diff(lie[start9:stop9, 9])*(lie[start9:(stop9-1), 9]-lie[lie.lines+1, 9])
                                            #       #               d2<-diff(lie[start9:stop9, 7])^2+diff(lie[start9:stop9, 8])^2+diff(lie[start9:stop9, 9])^2
                                            #       #               t<-(-(scalx+scaly+scalz)/d2)
                                            #       #               
                                            #       #               if (length(which(t>=0 & t<=1))==0){
                                            #       #                 
                                            #       #                 index<-which(alldist==min(alldist))
                                            #       #                 dbase<-lie[start9+index-1, 10]
                                            #       #                 if (connect==TRUE){
                                            #       #                   lie[lie.lines+1,5]<-lie[start9+index-1, 1] #Update Prec
                                            #       #                   dist1<-min(alldist)
                                            #       #                   lie[(lie.lines+1):(lie.lines+length10), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length10), 10]}}
                                            #       #               
                                            #       #               else {
                                            #       #                 
                                            #       #                 t[t<0]<-NA
                                            #       #                 t[t>1]<-NA
                                            #       #                 xn<-diff(lie[start9:stop9, 7])*t+lie[start9:(stop9-1), 7]
                                            #       #                 yn<-diff(lie[start9:stop9, 8])*t+lie[start9:(stop9-1), 8]
                                            #       #                 zn<-diff(lie[start9:stop9, 9])*t+lie[start9:(stop9-1), 9]
                                            #       #                 dist1<-sqrt((xn-lie[lie.lines+1, 7])^2+(yn-lie[lie.lines+1, 8])^2+(zn-lie[lie.lines+1, 9])^2)
                                            #       #                 
                                            #       #                 if (sum(is.na(dist1)==T)>0) {
                                            #       #                   index<-as.numeric(match(min(dist1, na.rm=T), dist1))
                                            #       #                   dist1<-min(dist1, na.rm=T)} 
                                            #       #                 else {
                                            #       #                   index<-as.numeric(match(min(dist1), dist1))
                                            #       #                   dist1<-min(dist1)}
                                            #       #                 
                                            #       #                 if (dist1>min(alldist)){
                                            #       #                   index<-which(alldist==min(alldist))
                                            #       #                   dist1<-min(alldist)
                                            #       #                   dbase<-lie[start9+index-1, 10]
                                            #       #                   if (connect==TRUE){
                                            #       #                     lie[lie.lines+1,5]<-lie[start9+index-1, 1] #Update Prec
                                            #       #                     dist1<-min(alldist)
                                            #       #                     lie[(lie.lines+1):(lie.lines+length10), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length10), 10]}}
                                            #       #                 
                                            #       #                 else{
                                            #       #                   
                                            #       #                   lie[(lie.lines+1):(lie.lines+length10), 10]<-dist1+lie[(lie.lines+1):(lie.lines+length10), 10]
                                            #       #                   dbase<-lie[start9+index-1, 10]+distance3D(x1=lie[start9+index-1, 7], x2=xn[index], y1=lie[start9+index-1, 8], y2=yn[index], z1=lie[start9+index-1, 9], z2=zn[index])
                                            #       #                   
                                            #       #                   if (connect==TRUE){
                                            #       #                     
                                            #       #                     lie[lie.lines+length10+1,1:13]<-c(NA, NA, 0, 0, NA, NA, xn[index], yn[index], zn[index], dbase, lie[start9+index-1, 11], NA, 9)
                                            #       #                     lie<-lie[order(lie[, 11], lie[, 10]),]
                                            #       #                     length9<-length9+1
                                            #       #                     lie.lines<-lie.lines+1
                                            #       #                     stop9<-stop9+1
                                            #       #                     pos<-match(NA, lie[1:(lie.lines),1])
                                            #       #                     lie[1:(lie.lines+length10),5]<-match(lie[1:(lie.lines+length10),5], lie[1:(lie.lines+length10),1])
                                            #       #                     lie[1:(lie.lines+length10),6]<-match(lie[1:(lie.lines+length10),6], lie[1:(lie.lines+length10),1])
                                            #       #                     lie[pos, 5]<-pos-1
                                            #       #                     lie[pos, 6]<-pos+1
                                            #       #                     lie[pos, 2]<-lie[lie[pos, 6], 2]
                                            #       #                     lie[pos+1, 5]<-pos
                                            #       #                     lie[pos-1, 6]<-pos
                                            #       #                     lie[which(is.na(lie[1:(lie.lines),5])==TRUE), 5]<-0
                                            #       #                     lie[which(is.na(lie[1:(lie.lines),6])==TRUE), 6]<-0
                                            #       #                     lie[1:(lie.lines+length10),1]<-match(lie[1:(lie.lines+length10),1], lie[1:(lie.lines+length10),1])
                                            #       #                     lie[lie.lines+1,5]<-pos
                                            #       #                     
                                            #       #                     #Calculate diameter for interpolated point (linear function)
                                            #       #                     if (is.null(diameter)==FALSE){
                                            #       #                       prec<-lie[pos, 5]
                                            #       #                       suiv<-lie[pos, 6]
                                            #       #                       slope<-(lie[suiv, 12]-lie[prec, 12])/distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[suiv, 7], y2=lie[suiv, 8], z2=lie[suiv, 9])
                                            #       #                       intercept<-lie[prec, 12]
                                            #       #                       lie[pos, 12]<-intercept+slope*distance3D(x1=lie[prec, 7], y1=lie[prec, 8], z1=lie[prec, 9], x2=lie[pos, 7], y2=lie[pos, 8], z2=lie[pos, 9])}}}}}
                                            #       #             
                                            #       #             if (connect==FALSE) {lie[lie.lines+1,5]<-lie.lines+1}
                                            #       #             
                                            #       #             lie[lie.lines+1,3]<-1
                                            #       #             start10<-as.numeric(lie[lie.lines+1, 1])
                                            #       #             
                                            #       #             # Change Suiv and Apic values for the last point of a lateral root
                                            #       #             lie[lie.lines+length10,c(4,6)]<-c(1, 0)
                                            #       #             stop10<-as.numeric(lie[lie.lines+length10, 1])
                                            #       #             
                                            #       #             # Fill RAC file for the 9-order root
                                            #       #             
                                            #       #             if (connect==TRUE){cumulDist<-dist1 + sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length10), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length10), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length10), 9]))^2))}
                                            #       #             else {cumulDist<-sum(sqrt((diff(lie[(lie.lines+1):(lie.lines+length10), 7]))^2+(diff(lie[(lie.lines+1):(lie.lines+length10), 8]))^2+(diff(lie[(lie.lines+1):(lie.lines+length10), 9]))^2))}
                                            #       #             
                                            #       #             rac[r, 1:6]<-c(max(rac[,1], na.rm=TRUE)+1, currentMother5, 10, dbase, 1, cumulDist)
                                            #       #             
                                            #       #             lie.lines<-lie.lines+length10
                                            #       #             
                                            #       # 
                                            #       #           }
                                            #       #         }
                                            #       #       }
                                            #       #       
                                            #       #     }
                                            #       #   }
                                            #       # }
                                            #       
                                            #     }
                                            #   }
                                            # }
                                            
                                          }
                                        }
                                      }
                                      
                                    }
                                  }
                                }
                              } 
                            }
                          } 
                        } 
                      }
                    } 
                  } 
                }
              } 
            } 
          }
        }
      }
    }
    
    #Remove lines with NA values in lie
    index<-which(is.na(lie[,1])==TRUE)
    if (length(index)>0) {lie<-lie[-index,]}
    
    #Remove lines with NA values in rac
    index<-which(is.na(rac[,1])==TRUE)
    if (length(index)>0) {
      rac<-rac[-index,]
      message(paste("Roots with a branching order greater than 7 have been skipped in ", sub(basename(rsml.path), pattern=".rsml", replacement=""), sep=""))}
    
    #Convert matrices to data frames
    lie<-as.data.frame(lie)
    rac<-as.data.frame(rac)
    colnames(lie)<-c("Num","Date", "Bran", "Apic", "Prec", "Suiv", "X", "Y", "Z", "dist", "root", "diameter", "ord")
    colnames(rac)<-c("Root", "Mother", "Ord", "DBase", "DApp", "Lengths1")
    
    lie$Bran[lie$Bran==0]<-"false"
    lie$Bran[lie$Bran==1]<-"true"
    lie$Apic[lie$Apic==0]<-"false"
    lie$Apic[lie$Apic==1]<-"true"
    
    if (sum(lie$Z)==0){lie<-lie[,-9]}
    
    #Round dates to integers
    lie[,2]<-ceiling(lie[,2])
    
    #Create TPS file
    
    dates<-as.numeric(unique(lie[,2]))
    dates<-sort(dates)
    
    if (is.null(unittime1)==TRUE){unittime1<-"unittime"}
    
    if (timeserie==FALSE) {
      if (is.null(final.date)==TRUE) {tps<-data.frame(Num=1, Date=1)}
      if (is.null(final.date)==FALSE) {tps<-data.frame(Num=1, Date=final.date)}}
    
    if (timeserie==TRUE) {tps<-data.frame(Num=c(1:length(dates)), Date=dates)}
    
    #Replace Date by Num in LIE file
    date.lie<-lie$Date
    
    lie$Date<-match(lie$Date, tps$Date)
    
    lie$Date[which(is.na(lie$Date)==TRUE)]<-0

    #Make RAC file for time series
    
    rac$Root<-rac$Root-1
    rac$Mother[rac$Mother!=-1]<-rac$Mother[rac$Mother!=-1]-1
    
    if (timeserie==TRUE){
      
      rac<-rac[,-6] #Delete Lengths column because must be replaced by a matrix
      cols<-nrow(tps)
      rows<-nrow(rac)
      length1<-matrix(0, nrow=rows, ncol=cols)
      
      maxlength<-aggregate(lie$dist, by=list(lie$root, lie$Date), max)
      length1[as.matrix(maxlength[,1:2])]<-maxlength$x
      
      length1<-t(apply(length1, 1, function(x){for (i in 2:length(x)){if (x[i]<x[i-1]){x[i]<-x[i-1]}};return(x)}))
      
      colnames(length1)<-paste(rep("Lengths", cols), c(1:cols), sep="")
      
      rac<-cbind(rac, as.data.frame(length1))}
    
    #Export RAC, LIE, TPS
    lie.all[[n]]<-lie
    rac.all[[n]]<-rac
    tps.all[[n]]<-tps
    
  }
  result<-list(resolution=resolution, length=unitlength1, time=unittime1, lie=lie.all, rac=rac.all, tps=tps.all)
  return(result)
}
