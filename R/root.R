root<-function(x, unitangle="d", vertical3d="y", last=TRUE, show.progress=FALSE){
  
  #x must be an object of class "rsmlToTable" or "dartToTable"
  if ("rsmlToTable" %in% class(x) | "dartToTable" %in% class(x)) {} else {stop("x must be a rsmlToTable or dartToTable object")}
  
  if (mode(unitangle)!="character"){stop("mode(unitangle) must be character")}
  if(unitangle=="d"|unitangle=="r") {} else {stop("unitangle must be either d (degrees) or r (radians)")}
  
  if (mode(show.progress)!="logical"){stop("show.progress must be logical")}
  
  if (mode(last)!="logical"){stop("last must be logical")}
  
  if (vertical3d=="x"|vertical3d=="y"|vertical3d=="z") {} else {stop("vertical3d must be x, y, or z")}
  
  if ("dartToTable" %in% class(x)) {rootsystem<-x$file}
  if ("rsmlToTable" %in% class(x)) {rootsystem<-paste(x$file, x$plant, sep="_")}
  RSlevels<-unique(rootsystem)
  
  n<-length(RSlevels) #Number of root systems in x
  
  # Unit conversion angles
  
  if (unitangle=="r") {cunitangle<-1}
  if (unitangle=="d") {cunitangle<-180/pi}
  
  # Vertical direction vector
  
  if (vertical3d=="x") {dirvert<-c(1,0,0)}
  if (vertical3d=="y") {dirvert<-c(0,1,0)}
  if (vertical3d=="z") {dirvert<-c(0,0,1)}
  
  if (last==FALSE){
  
      #Compute root parameters
      
      if (show.progress==TRUE) {pb<-txtProgressBar(min=1, max=n, style=3)}
      
      results<-list()
      filenames<-c()
      
      allfiles<-x$file
      x$bran[x$bran=="true"]<-1
      x$bran[x$bran=="false"]<-0
      x$apic[x$apic=="true"]<-1
      x$apic[x$apic=="false"]<-0
      
      if ("rsmlToTable" %in% class(x)){
          x<-data.matrix(x[,2:ncol(x)])
          class(x)<-c("matrix", "rsmlToTable")}
      
      if ("dartToTable" %in% class(x)){
        x<-data.matrix(x[,2:ncol(x)])
        class(x)<-c("matrix", "dartToTable")}
        
      for (i in 1:n){#For each root system in x
        
        table<-x[which(rootsystem==RSlevels[i]),] #Create a matrix subset for computing root stats
        dates<-unique(table[,"time"]) #Search for all observation dates
        allroots<-unique(table[,"root"]) #Search for all roots in root system
        filenames<-append(filenames, rep(unique(allfiles[which(rootsystem==RSlevels[i])]), length(allroots)*length(dates)))
        
        if ("rsmlToTable" %in% class(x)){
          results[[i]]<-matrix(nrow=length(allroots)*length(dates), ncol=16)} #Create matrix to store the results

        if ("dartToTable" %in% class(x)){
          results[[i]]<-matrix(nrow=length(allroots)*length(dates), ncol=15)} #Create matrix to store the results

        line<-0 #Count line number in matrix storing the results
        
        if (show.progress==TRUE) {setTxtProgressBar(pb, i)}
        
        for (root in 1:length(allroots)){#For each root
          
          index<-which(table[,"root"]==allroots[root])
          order<-unique(table[index, "order"])
          parent<-unique(table[index, "parentroot"])
          if ("rsmlToTable" %in% class(x)){plant<-unique(table[index, "plant"])}
          
          #DBase
          if (order==1) {dbase<-0}
          else {
            start<-table[which(table[,"root"]==allroots[root] & table[, "bran"]==1), c("x1", "y1", "z1")]
            
            if (length(which(table[,"root"]==parent & table[,"x2"]==start[1] & table[,"y2"]==start[2] & table[,"z2"]==start[3]))>0){
              
              dbase<-table[which(table[,"root"]==parent & table[,"x2"]==start[1] & table[,"y2"]==start[2] & table[,"z2"]==start[3]),"blength"]}
            
            else {#Lateral roots not connected to mother
              
              posparent<-which(table[,"root"]==parent)
              scalx<-(table[posparent,"x2"]-table[posparent,"x1"])*(table[posparent, "x1"]-start[1])
              scaly<-(table[posparent,"y2"]-table[posparent,"y1"])*(table[posparent, "y1"]-start[2])
              scalz<-(table[posparent,"z2"]-table[posparent,"z1"])*(table[posparent, "z1"]-start[3])
              d2<-(table[posparent, "x2"]-table[posparent, "x1"])^2+(table[posparent, "y2"]-table[posparent, "y1"])^2+(table[posparent, "z2"]-table[posparent, "z1"])^2
              t<-(-(scalx+scaly+scalz)/d2)
              
              if (length(which(t>=0 & t<=1))==0){
                
                index<-min(which(t<0))
                dbase<-table[posparent[1]+index-1, "blength"]}

              else {
                
                t[t<0]<-NA
                t[t>1]<-NA
                xn<-(table[posparent, "x2"]-table[posparent, "x1"])*t+table[posparent, "x1"]
                yn<-(table[posparent, "y2"]-table[posparent, "y1"])*t+table[posparent, "y1"]
                zn<-(table[posparent, "z2"]-table[posparent, "z1"])*t+table[posparent, "z1"]
                dist1<-sqrt((xn-start[1])^2+(yn-start[2])^2+(zn-start[3])^2)
                
                if (sum(is.na(dist1)==T)>0) {index<-as.numeric(match(min(dist1, na.rm=T), dist1))}
                else {index<-as.numeric(match(min(dist1), dist1))}
              
                dbase<-table[posparent[1]+index-1, "blength"]-table[posparent[1]+index-1, "length"]+distance3D(x1=table[posparent[1]+index-1, "x1"], x2=xn[index], y1=table[posparent[1]+index-1, "y1"], y2=yn[index], z1=table[posparent[1]+index-1, "z1"], z2=zn[index])}}}
          
          #Branching angle
          if (order==1){ #If root branching order = 1
            nodes<-sum(table[,"root"]==allroots[root])+1
            
            d1<-table[which(table[, "root"]==allroots[root] & table[, "bran"]==1), c("x1", "y1", "z1")]
            if (nodes>=4) {d2<-table[which(table[, "root"]==allroots[root] & table[, "bran"]==1)+2, c("x2", "y2", "z2")]}
            if (nodes<4) {d2<-table[which(table[, "root"]==allroots[root] & table[, "apic"]==1), c("x2", "y2", "z2")]}
            dird<-as.vector(d2-d1)
            angle<-acos((dird%*%dirvert)/sqrt(sum(dird^2)))*cunitangle}
          
          else{ #If root branching order > 1 
            
            nodes<-sum(table[,"root"]==allroots[root])+1
            d1<-table[which(table[, "root"]==allroots[root] & table[, "bran"]==1), c("x1", "y1", "z1")]
            if (nodes>=4) {d2<-table[which(table[, "root"]==allroots[root] & table[, "bran"]==1)+2, c("x2", "y2", "z2")]}
            if (nodes<4) {d2<-table[which(table[, "root"]==allroots[root] & table[, "apic"]==1), c("x2", "y2", "z2")]}
            xyzparent<-table[which(table[, "root"]==parent) , c("x1", "y1", "z1")]
            distance<-sqrt((xyzparent[,1]-d1[1])^2+(xyzparent[,2]-d1[2])^2+(xyzparent[,3]-d1[3])^2)
            pos<-which(table[, "root"]==parent & table[, "bran"]==1)+which(distance==min(distance))-1
            
            if (table[pos, "bran"]!=1){
              m2<-table[pos, c("x2", "y2", "z2")]
              m1<-table[pos-1, c("x1", "y1", "z1")]}
            else{
              m1<-table[pos, c("x1", "y1", "z1")]
              m2<-table[pos+1, c("x2", "y2", "z2")]}
            
            dird<-as.vector(d2-d1)
            dirm<-as.vector(m2-m1)
            angle<-acos((dird%*%dirm)/(sqrt(sum(dird^2))*sqrt(sum(dirm^2))))*cunitangle}
          
            for (t in 1:length(dates)){#For each observation dates
              
                line<-line+1
                length<-sum(table[which(table[,"root"]==allroots[root] & table[,"time"]<=dates[t]),"length"])
                nlat<-length(unique(table[which(table[,"parentroot"]==allroots[root] & table[,"time"]<=dates[t]),"root"]))
              
              if ("rsmlToTable" %in% class(x)){
                
                surface<-sum(table[which(table[,"root"]==allroots[root] & table[,"time"]<=dates[t]), "surface"])
                volume<-sum(table[which(table[,"root"]==allroots[root] & table[,"time"]<=dates[t]), "volume"])

                #Diameter
                if (length>0){
                  meandiam<-mean(c(table[which(table[,"root"]==allroots[root] & table[,"time"]<=dates[t]),"diameter2"], table[which(table[,"root"]==allroots[root] & table[,"bran"]==1),"diameter1"]))
                  sddiam<-sd(c(table[which(table[,"root"]==allroots[root] & table[,"time"]<=dates[t]),"diameter2"], table[which(table[,"root"]==allroots[root] & table[,"bran"]==1),"diameter1"]))}
                else {
                  meandiam<-NA
                  sddiam<-NA}}
                
              if ("dartToTable" %in% class(x)){
                
                surface<-NA
                volume<-NA
                meandiam<-NA
                sddiam<-NA}
                
                #Tortuosity
                if (length==0){tortuosity<-NA}
                else {
                  xyzstart<-as.vector(table[which(table[,"root"]==allroots[root] & table[,"bran"]==1), c("x1", "y1", "z1")])
                  xyzend<-as.vector(table[max(which(table[,"root"]==allroots[root] & table[,"time"]<=dates[t])), c("x2", "y2", "z2")])
                  diff<-xyzend-xyzstart
                  tortuosity<-length/sqrt(sum(diff^2))}
                
                #Growth rate
                if (t==1){growth<-length/dates[t]} else {growth<-(length-results[[i]][line-1, 6])/(dates[t]-dates[t-1])}
                
                #Length apical unbranched zone (luaz). Cfr Pagès et al (2010) Plant and Soil, 328, 35-44.
                #Can only be calculated if rsml.connect=TRUE
                if (length==0) {luaz<-0}
                else {
                  if (nlat==0){luaz<-length}
                  else{
                    
                    laterals<-table[which(table[,"parentroot"]==allroots[root] & table[,"time"]<=dates[t] & table[,"bran"]==1), c("x1", "y1", "z1")]
                    if (is.matrix(laterals)==FALSE){laterals<-matrix(laterals, ncol=3)}
                    index<-apply(laterals, 1, function(x){which(table[,"root"]==allroots[root] & table[,"x1"]==x[1] & table[,"y1"]==x[2] & table[,"z1"]==x[3])})
                    if (length(index)>0) {luaz<-length-max(table[index, "blength"]-table[index, "length"])} else {luaz<-NA}}}
                
                #Fill matrix
                if  ("rsmlToTable" %in% class(x)) {results[[i]][line, c(1:16)]<-c(allroots[root], dates[t], order, parent, dbase, length, meandiam, sddiam, nlat, angle, tortuosity, growth, surface, volume, luaz, plant)}
                if  ("dartToTable" %in% class(x)) {results[[i]][line, c(1:15)]<-c(allroots[root], dates[t], order, parent, dbase, length, meandiam, sddiam, nlat, angle, tortuosity, growth, surface, volume, luaz)}}}}}
  
  if (last==TRUE){
    
    #Compute root parameters for last observation date only
    
    if (show.progress==TRUE) {pb<-txtProgressBar(min=1, max=n, style=3)}
    
    results<-list()
    filenames<-c()
    
    allfiles<-x$file
    x$bran[x$bran=="true"]<-1
    x$bran[x$bran=="false"]<-0
    x$apic[x$apic=="true"]<-1
    x$apic[x$apic=="false"]<-0
    
    if ("rsmlToTable" %in% class(x)){
      x<-data.matrix(x[,2:ncol(x)])
      class(x)<-c("matrix", "rsmlToTable")}
    
    if ("dartToTable" %in% class(x)){
      x<-data.matrix(x[,2:ncol(x)])
      class(x)<-c("matrix", "dartToTable")}
    
    for (i in 1:n){#For each root system in x
      
      table<-x[which(rootsystem==RSlevels[i]),] #Create a matrix subset for computing root stats
      allroots<-unique(table[,"root"]) #Search for all roots in root system
      filenames<-append(filenames, rep(unique(allfiles[which(rootsystem==RSlevels[i])]), length(allroots)))

    if ("rsmlToTable" %in% class(x)){
      results[[i]]<-matrix(nrow=length(allroots), ncol=16)} #Create matrix to store the results
      
    if ("dartToTable" %in% class(x)){
      results[[i]]<-matrix(nrow=length(allroots), ncol=15)} #Create matrix to store the results
      
      line<-0 #Count line number in matrix storing the results
      
      if (show.progress==TRUE) {setTxtProgressBar(pb, i)}
      
      for (root in 1:length(allroots)){#For each root
        
        #Branching angle (Just once, does not depend on time!)
        
        index<-which(table[,"root"]==allroots[root])
        order<-unique(table[index, "order"])
        parent<-unique(table[index, "parentroot"])
        if ("rsmlToTable" %in% class(x)){plant<-unique(table[index, "plant"])}
        
        #DBase
        if (order==1) {dbase<-0}
        else {
          start<-table[which(table[,"root"]==allroots[root] & table[, "bran"]==1), c("x1", "y1", "z1")]
          
          if (length(which(table[,"root"]==parent & table[,"x2"]==start[1] & table[,"y2"]==start[2] & table[,"z2"]==start[3]))>0){
            
            dbase<-table[which(table[,"root"]==parent & table[,"x2"]==start[1] & table[,"y2"]==start[2] & table[,"z2"]==start[3]),"blength"]}
          
          else {#Lateral roots not connected to mother
            
            posparent<-which(table[,"root"]==parent)
            scalx<-(table[posparent,"x2"]-table[posparent,"x1"])*(table[posparent, "x1"]-start[1])
            scaly<-(table[posparent,"y2"]-table[posparent,"y1"])*(table[posparent, "y1"]-start[2])
            scalz<-(table[posparent,"z2"]-table[posparent,"z1"])*(table[posparent, "z1"]-start[3])
            d2<-(table[posparent, "x2"]-table[posparent, "x1"])^2+(table[posparent, "y2"]-table[posparent, "y1"])^2+(table[posparent, "z2"]-table[posparent, "z1"])^2
            t<-(-(scalx+scaly+scalz)/d2)
            
            if (length(which(t>=0 & t<=1))==0){
              
              index<-min(which(t<0))
              dbase<-table[posparent[1]+index-1, "blength"]}
            
            else {
              
              t[t<0]<-NA
              t[t>1]<-NA
              xn<-(table[posparent, "x2"]-table[posparent, "x1"])*t+table[posparent, "x1"]
              yn<-(table[posparent, "y2"]-table[posparent, "y1"])*t+table[posparent, "y1"]
              zn<-(table[posparent, "z2"]-table[posparent, "z1"])*t+table[posparent, "z1"]
              dist1<-sqrt((xn-start[1])^2+(yn-start[2])^2+(zn-start[3])^2)
              
              if (sum(is.na(dist1)==T)>0) {index<-as.numeric(match(min(dist1, na.rm=T), dist1))}
              else {index<-as.numeric(match(min(dist1), dist1))}
              
              dbase<-table[posparent[1]+index-1, "blength"]-table[posparent[1]+index-1, "length"]+distance3D(x1=table[posparent[1]+index-1, "x1"], x2=xn[index], y1=table[posparent[1]+index-1, "y1"], y2=yn[index], z1=table[posparent[1]+index-1, "z1"], z2=zn[index])}}}
        
        #Branching angle
        
        if (order==1){
          nodes<-sum(table[,"root"]==allroots[root])+1
          d1<-table[which(table[, "root"]==allroots[root] & table[, "bran"]==1), c("x1", "y1", "z1")]
          if (nodes>=4) {d2<-table[which(table[, "root"]==allroots[root] & table[, "bran"]==1)+2, c("x2", "y2", "z2")]}
          if (nodes<4) {d2<-table[which(table[, "root"]==allroots[root] & table[, "apic"]==1), c("x2", "y2", "z2")]}
          dird<-as.vector(d2-d1)
          angle<-acos((dird%*%dirvert)/sqrt(sum(dird^2)))*cunitangle}
        
        else{
          nodes<-sum(table[,"root"]==allroots[root])+1
          d1<-table[which(table[, "root"]==allroots[root] & table[, "bran"]==1), c("x1", "y1", "z1")]
          if (nodes>=4) {d2<-table[which(table[, "root"]==allroots[root] & table[, "bran"]==1)+2, c("x2", "y2", "z2")]}
          if (nodes<4) {d2<-table[which(table[, "root"]==allroots[root] & table[, "apic"]==1), c("x2", "y2", "z2")]}
          xyzparent<-table[which(table[, "root"]==parent) , c("x1", "y1", "z1")]
          distance<-sqrt((xyzparent[,1]-d1[1])^2+(xyzparent[,2]-d1[2])^2+(xyzparent[,3]-d1[3])^2)
          pos<-which(table[, "root"]==parent & table[, "bran"]==1)+which(distance==min(distance))-1
          
          if (table[pos, "bran"]!=1){
            m2<-table[pos, c("x2", "y2", "z2")]
            m1<-table[pos-1, c("x1", "y1", "z1")]}
          else{
            m1<-table[pos, c("x1", "y1", "z1")]
            m2<-table[pos+1, c("x2", "y2", "z2")]}
          
          dird<-as.vector(d2-d1)
          dirm<-as.vector(m2-m1)
          angle<-acos((dird%*%dirm)/(sqrt(sum(dird^2))*sqrt(sum(dirm^2))))*cunitangle}
        
          line<-line+1
          length<-sum(table[which(table[,"root"]==allroots[root]),"length"])
          nlat<-length(unique(table[which(table[,"parentroot"]==allroots[root]),"root"]))
          
        if ("rsmlToTable" %in% class(x)){
          
          surface<-sum(table[which(table[,"root"]==allroots[root]), "surface"])
          volume<-sum(table[which(table[,"root"]==allroots[root]), "volume"])
          
          #Mean diameter
          if (length>0){
            meandiam<-mean(c(table[which(table[,"root"]==allroots[root]),"diameter2"], table[which(table[,"root"]==allroots[root] & table[,"bran"]==1),"diameter1"]))
            sddiam<-sd(c(table[which(table[,"root"]==allroots[root]),"diameter2"], table[which(table[,"root"]==allroots[root] & table[,"bran"]==1),"diameter1"]))}
          else {
            meandiam<-NA
            sddiam<-NA}}
          
          if ("dartToTable" %in% class(x)){
            surface<-NA
            volume<-NA
            meandiam<-NA
            sddiam<-NA}
          
          #Tortuosity
          if (length==0){tortuosity<-NA}
          else {
            xyzstart<-as.vector(table[which(table[,"root"]==allroots[root] & table[,"bran"]==1), c("x1", "y1", "z1")])
            xyzend<-as.vector(table[max(which(table[,"root"]==allroots[root])), c("x2", "y2", "z2")])
            diff<-xyzend-xyzstart
            tortuosity<-length/sqrt(sum(diff^2))}
          
          #Growth rate
          growth<-length/max(table[,"time"])
          
          #Length apical unbranched zone (lauz). Cfr Pagès et al (2010) Plant and Soil, 328, 35-44.
          if (length==0) {luaz<-0}
          else {
            if (nlat==0){luaz<-length}
            else{
              laterals<-table[which(table[,"parentroot"]==allroots[root] & table[,"bran"]==1), c("x1", "y1", "z1")]
              if (is.matrix(laterals)==FALSE){laterals<-matrix(laterals, ncol=3)}
              index<-apply(laterals, 1, function(x){which(table[,"root"]==allroots[root] & table[,"x1"]==x[1] & table[,"y1"]==x[2] & table[,"z1"]==x[3])})
              if (length(index)>0) {luaz<-length-max(table[index, "blength"]-table[index, "length"])} else {luaz<-NA}}}
          
          #Fill matrix
          if ("rsmlToTable" %in% class(x)) {results[[i]][line, c(1:16)]<-c(allroots[root], max(table[,"time"]), order, parent, dbase, length, meandiam, sddiam, nlat, angle, tortuosity, growth, surface, volume, luaz, plant)}
          if ("dartToTable" %in% class(x)) {results[[i]][line, c(1:15)]<-c(allroots[root], max(table[,"time"]), order, parent, dbase, length, meandiam, sddiam, nlat, angle, tortuosity, growth, surface, volume, luaz)}}}}
  
  results<-as.data.frame(do.call(rbind, results))

  if ("rsmlToTable" %in% class(x)){
    colnames(results)<-c("root", "time", "order", "parentroot", "DBase", "length", "mean.diameter", "sd.diameter", "nlat", "branching.angle", "tortuosity", "growth", "surface", "volume", "lauz", "plant")
    results$file<-filenames
    results<-results[,c(17,16,1,3,4,2,5:15)]}
  
  if ("dartToTable" %in% class(x)){
    colnames(results)<-c("root", "time", "order", "parentroot", "DBase", "length", "mean.diameter", "sd.diameter", "nlat", "branching.angle", "tortuosity", "growth", "surface", "volume", "lauz")
    results$file<-filenames
    results<-results[,c(16,1,3,4,2,5:15)]
    results<-results[,-which(colnames(results)=="mean.diameter"|colnames(results)=="sd.diameter"|colnames(results)=="surface"|colnames(results)=="volume")]}
  
  return(results)}