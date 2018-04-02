bottleneckdist<-function(x, show.progress=FALSE){
  
  #x is a list created with the perhomology function
  
  message("bottleneckdist requires the TDA package to be installed")
  
  if ("perhomology" %in% class(x)) {} else {stop("x must be a perhomology object")}
  
  if (mode(show.progress)!="logical"){stop("show.progress must be logical")}
  
  n<-length(x)
  
  distance<-matrix(nrow=n, ncol=n)
  colnames(distance)<-names(x)
  rownames(distance)<-names(x)
  
  if (requireNamespace("TDA", quietly=TRUE)){
  
      if (show.progress==TRUE) {pb<-txtProgressBar(min=1, max=n, style=3)}
      
      for (r in 1:n){
        
        if (show.progress==TRUE) {setTxtProgressBar(pb, r)}
        
        for (c in r:n){
          
          distance[r,c]<-TDA::bottleneck(Diag1=x[[r]], Diag2=x[[c]], dimension=0)} #From TDA package
        
        distance[,r]<-t(distance[r,])}} #Because the matrix is symmetric
  
  else{TDA::bottleneck}
  
  return(distance)}