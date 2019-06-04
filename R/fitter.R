fitter<-function(x){

x1<-x  #Make a copy

#x is a matrix (for 1 RSML, can contain several plants)
#Add 2 columns to x (magnitude and pathlength)
x1<-cbind(x1, matrix(c(0,1), nrow=nrow(x1), ncol=2, byrow=TRUE))  
plants<-unique(x1[,2]) #Find all plants in RSML

for (n in 1:length(plants)){#For each plant in RSML
  
  x<-x1[which(x1[,2]==plants[n]),] #Subset for each plant
  if (is.null(dim(x))==TRUE){x<-matrix(x, ncol=23)}
  apicindex<-which(x[,7]==1) #apicindex and branindex should have the same length
  branindex<-which(x[,6]==1)

  for (l in 1:length(apicindex)){#For each apic point of a root system
    
    #Magnitude
    
    x[apicindex[l], 22]<-1 #External link has a magnitude of 1
    
    root<-x[apicindex[l], 3] #Take root ID
    parentroot<-x[apicindex[l], 21] #Take parent root ID
    
    indexprec<-which(x[,11]==x[apicindex[l], 8] & x[,12]==x[apicindex[l], 9] & x[,13]==x[apicindex[l], 10]) #Search prec segment
    if (length(indexprec)>1) {indexprec<-indexprec[which(x[indexprec,3]==root | x[indexprec, 3]==parentroot)]} #If more than 1 segment found, it should belong to root or parentroot
    
    if (length(indexprec)>0){
    
      root<-x[indexprec, 3] #Take root ID of prec segment
      parentroot<-x[indexprec, 21] #Take parent root ID of the prec segment
      
      x[indexprec, 22]<-x[indexprec, 22]+1
      
      while(x[indexprec, 4]>=1){
        
        segment1<-which(x[,11]==x[indexprec, 8] & x[,12]==x[indexprec, 9] & x[,13]==x[indexprec, 10])
        if (length(segment1)>1) {indexprec<-segment1[which(x[segment1,3]==root | x[segment1, 3]==parentroot)]} 
        else {indexprec<-segment1}
        if (length(indexprec)==0){break}
        root<-x[indexprec, 3]
        parentroot<-x[indexprec, 21]
          
        x[indexprec, 22]<-x[indexprec, 22]+1}}
    
    #Path length
    
    root<-x[branindex[l], 3]
    
    testbran<-which(x[,8]==x[branindex[l], 11] & x[,9]==x[branindex[l], 12] & x[,10]==x[branindex[l], 13]) #Is it a crossing?
    if (length(testbran)>1) {testbran<-testbran[which(x[testbran, 3]==root | x[testbran, 21]==root)]} #Select segments based on root ID and parentroot ID

    if (length(testbran)==0) {} else{
    
    if (length(testbran)>=2) {
      x[testbran, 23]<-x[branindex[l], 23]+1
      index<-which(x[testbran, 6]==0)
      suiv<-testbran[index]}
    else {
      x[testbran, 23]<-x[branindex[l], 23]
      suiv<-testbran}
    
    while(x[suiv, 7]==0){
      
      testbran<-which(x[,8]==x[suiv, 11] & x[,9]==x[suiv, 12] & x[,10]==x[suiv, 13]) #Is it a crossing?
      if (length(testbran)>1) {testbran<-testbran[which(x[testbran, 3]==root | x[testbran, 21]==root)]} #Select segments based on root ID and parentroot ID

      if (length(testbran)>=2) {
        x[testbran, 23]<-x[suiv, 23]+1
        index<-which(x[testbran, 6]==0)
        suiv<-testbran[index]}
      else {
        x[testbran, 23]<-x[suiv, 23]
        suiv<-testbran}}}}
  
  x1[x1[,2]==plants[n],]<-x}

  return(x1)}