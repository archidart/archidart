fitter<-function(x){

#x is a matrix
#Add 2 columns to x (magnitude and pathlength)
x<-cbind(x, matrix(c(0,1), nrow=nrow(x), ncol=2, byrow=TRUE))  

  apicindex<-which(x[,7]==1) #apicindex and branindex should have the same length
  branindex<-which(x[,6]==1)
  
  for (l in 1:length(apicindex)){#For each apic point of a root system
    
    #Magnitude
    
    x[apicindex[l], 22]<-1
    
    root<-x[apicindex[l], 3]
    
    if (x[apicindex[l], 4]==1 & x[apicindex[l], 6]==1) {break}
    
    testbran<-which(x[,8]==x[apicindex[l], 8] & x[,9]==x[apicindex[l], 9] & x[,10]==x[apicindex[l], 10]) #Is it a crossing?
    testbran<-testbran[which(x[testbran,3]==root | x[testbran,21]==root)] #Select segments based on root and parentroot
    
    if (length(testbran)>=2){#We are at a crossing
      
      segment1<-which(x[,11]==x[apicindex[l], 8] & x[,12]==x[apicindex[l], 9] & x[,13]==x[apicindex[l], 10])
      segment1<-segment1[which(x[segment1,3]==root)]
      segment2<-testbran[testbran!=apicindex[l]]
      geo1<-x[segment1, 20]
      geo2<-x[segment2, 20]
      if (geo1>geo2){
        indexprec<-segment2
        root<-x[indexprec, 3]}
      if (geo1<geo2){
        indexprec<-segment1
        root<-x[indexprec, 3]}}
    
    else {
      indexprec<-which(x[,11]==x[apicindex[l], 8] & x[,12]==x[apicindex[l], 9] & x[,13]==x[apicindex[l], 10])
      root<-x[indexprec, 3]}
    
    x[indexprec, 22]<-x[indexprec, 22]+1
    
    while(x[indexprec, 4]>=1){
      
      if (x[indexprec, 4]==1 & x[indexprec, 6]==1) {break}
      
      testbran<-which(x[,8]==x[indexprec, 8] & x[,9]==x[indexprec, 9] & x[,10]==x[indexprec, 10])
      
      if (length(testbran)==2){#We are at a crossing
        
        segment1<-which(x[,11]==x[indexprec, 8] & x[,12]==x[indexprec, 9] & x[,13]==x[indexprec, 10])
        segment2<-testbran[testbran!=indexprec]
        geo1<-x[segment1, 20]
        geo2<-x[segment2, 20]
        if (geo1>geo2){indexprec<-segment2}
        if (geo1<geo2){indexprec<-segment1}}
      
      else {indexprec<-which(x[,11]==x[indexprec, 8] & x[,12]==x[indexprec, 9] & x[,13]==x[indexprec, 10])}
      
      x[indexprec, 22]<-x[indexprec, 22]+1}
    
    #Path length
    
    testbran<-which(x[,8]==x[branindex[l], 11] & x[,9]==x[branindex[l], 12] & x[,10]==x[branindex[l], 13])
    
    if (length(testbran)==2) {
      x[testbran, 23]<-x[branindex[l], 23]+1
      index<-which(x[testbran, 6]==0)
      suiv<-testbran[index]}
    else {
      x[testbran, 23]<-x[branindex[l], 23]
      suiv<-testbran}
    
    while(x[suiv, 7]==0){
      
      testbran<-which(x[,8]==x[suiv, 11] & x[,9]==x[suiv, 12] & x[,10]==x[suiv, 13])
      
      if (length(testbran)==2) {
        x[testbran, 23]<-x[suiv, 23]+1
        index<-which(x[testbran, 6]==0)
        suiv<-testbran[index]}
      else {
        x[testbran, 23]<-x[suiv, 23]
        suiv<-testbran}}}

  return(x)}