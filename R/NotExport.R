#Compute euclidean distance between 2 points in 2D

distance2D<-function(x1,y1,x2,y2){
  a<-sqrt((x1-x2)^2+(y1-y2)^2)
  return(a)}

#Compute euclidean distance between 2 points in 3D

distance3D<-function(x1,y1,z1,x2,y2,z2){
  a<-sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
  return(a)}

# In 2D
# Calculation of a normal vector to a line (ax+by+c=0)
# u is a direction vector of ax+by+c=0

normal<-function(u){
  n<-c(u[2], -u[1])
  return(n)}

#Export X, Y and Z coordinates of nodes from an RSML file

xnodes<-function(x){
  xnode<-as.numeric(x[[1]])
  return(xnode)}

ynodes<-function(x){
  ynode<-as.numeric(x[[2]])
  return(ynode)}

znodes<-function(x){
  znode<-as.numeric(x[[3]])
  return(znode)}

# In 2D
# Creation of a function to calculate the coordinates of an unknown point (Xn) on a line knowing two points (X1 and X2) 
# of that line and the distance between X1 and Xn

XYcoord2D<-function(x1,x2,y1,y2,d){
  
  distx1x2<-distance2D(x1=x1, y1=y1, x2=x2, y2=y2)
  
  if (d>=(distx1x2/2)) {
    t<-d/distx1x2
    xn<-(x2-x1)*t+x1
    yn<-(y2-y1)*t+y1
    result<-c(xn, yn)
    return(result)}
  
  if (d<(distx1x2/2)) {
    t<-(distx1x2-d)/distx1x2
    xn<-(x1-x2)*t+x2
    yn<-(y1-y2)*t+y2
    result<-c(xn, yn)
    return(result)}}

# In 3D
# Creation of a function to calculate the coordinates of an unknown point (Xn) on a line knowing two points (X1 and X2) 
# of that line and the distance between X1 and Xn

XYcoord3D<-function(x1,x2,y1,y2,z1,z2,d){
  
  distx1x2<-distance3D(x1=x1, y1=y1, z1=z1, x2=x2, y2=y2, z2=z2)
  
  if (d>=(distx1x2/2)) {
    t<-d/distx1x2
    xn<-(x2-x1)*t+x1
    yn<-(y2-y1)*t+y1
    zn<-(z2-z1)*t+z1
    result<-c(xn, yn, zn)
    return(result)}
  
  if (d<(distx1x2/2)) {
    t<-(distx1x2-d)/distx1x2
    xn<-(x1-x2)*t+x2
    yn<-(y1-y2)*t+y2
    zn<-(z1-z2)*t+z2
    result<-c(xn, yn, zn)
    return(result)}}