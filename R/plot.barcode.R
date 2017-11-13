plot.barcode<-function(x, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, ...){
  if (is.null(xlim)==TRUE){xlim<-c(ceiling(max(x[,"birth"])), floor(min(x[,"death"])))}
  if (is.null(ylim)==TRUE){ylim<-c(nrow(x),1)}
  if (is.null(xlab)==TRUE){xlab<-"Distance function"}
  if (is.null(ylab)==TRUE){ylab<-expression(H[0])}
  plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  segments(x0=x[,"birth"], y0=c(1:nrow(x)), x1=x[,"death"], y1=c(1:nrow(x)), ...)}