#' plot.mscp
#'
#' Plot method for class 'mscp'
#'
#' @param x object of class mscp
#' @param cex numeric, global sizes in plot
#' @param plot.legend logical, if TRUE legends are plotted
#' @param ... additional arguments
#' 
#' @return No return value, called for side effects
#' 
#' @seealso \code{\link{mscp}, \link{summary.mscp}}
#'
#' @author Tijana Levajkovic and Michael Messer
#' 
#' @references Multiscale change point detection via gradual bandwidth adjustment in moving sum processes (2021+), Tijana Levajkovic and Michael Messer
#' 
#' @examples 
#' set.seed(1)
#' Tt <- 1000
#' cp <- c(250,500,600,650,750)
#' mu <- c(2,3,6,9,12,15)
#' sd <- c(1,1,2,1,2,1)
#' m  <- rep(mu,diff(c(0,cp,Tt))) 
#' s  <- rep(sd,diff(c(0,cp,Tt)))    
#' x  <- rnorm(Tt,m,s)
#' result <- mscp(x,kappa=4.77) # kappa set manually
#' # result <- mscp(x) # kappa derived in simulations
#' summary(result)
#' plot(result)
#' 
#' @rdname plot.mscp
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end



plot.mscp <- function(x=x,cex=1,plot.legend=TRUE,...)
{
  if(class(x) != "mscp"){stop("object must be of class 'mscp'")}
  if(!is.logical(plot.legend)){stop("plot.legend must be 'logical'")}
  
  # assign statistics
  object     <-  x
  x          <-  object$x
  mindist    <-  object$mindist
  delta      <-  object$delta
  S          <-  object$S
  path       <-  object$path
  cp         <-  object$cp
  mean_sd    <-  object$mean_sd
  
  # set layout parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  layout(matrix(c(1,1,2,3)))
  co <- c("magenta","red","blue","darkgreen")
  
  # plot 1: triangle 
  par(cex=cex,mar=c(1.5,4.5,0.5,0.5))
  plot(1,type="n",xlim=c(1,length(x)),ylim=c(1,length(x)/2),cex.axis=1,xlab="",ylab="",las=1,axes=FALSE)
  axis(2,las=1)
  polygon(x=c(delta,length(x)/2,length(x)-delta),y=c(delta,length(x)/2,delta))
  points(S[,1],S[,2],col=co[1],pch=19,cex=0.05)
  
  if(length(cp)>0){
  for (i in 1:length(cp))
    {
      ti <- path[[i]][,1]
      hi <- path[[i]][,2]
      lines(ti,hi,type="S",col=1,lwd=2)  
      points(ti[1],hi[1],pch=1,col=co[1],lwd=2)
      points(ti[length(ti)],hi[length(hi)],pch=4,col=co[4],cex=1,lwd=2)
    }
  }
    
  if(plot.legend)
  {
    if(length(cp)>0)
      {
      text(cp,rep(-length(x)*cex/30,length(cp)),1:length(cp),cex=seq(1,0.3,length.out = length(cp)),col=co[4],xpd=TRUE)  
      }
    xl <- 0.05; xr <- 0.1; yu <- 0.4; yo <- 0.5 
    arrows(xl*length(x),yu*length(x),xr*length(x),yu*length(x),length=0.05,lwd=2)
    arrows(xl*length(x),yu*length(x),xl*length(x),yo*length(x),length=0.05,lwd=2)
    text(((xl+xr)/2)*length(x),yu*length(x),"t",pos=1)
    text(xl*length(x),((yu+yo)/2)*length(x),"h",pos=2)
    legend(x="topright",legend=c("start","path","end","order"),pch=c(1,124,4,109),col=c(co[1],"black",co[4],co[4]),bty="n")
  }# end if plot.legend
  
  # plot 2: mean
  par(mar=c(0,4.5,1,0.5))
  plot(x,type="l",axes=FALSE,ylab=expression(paste (hat(mu))),col.lab=co[2]); 
  axis(3,padj = 0.5); 
  axis(2,las=1)
  lines(c(1,sort(cp),length(x)),c(mean_sd[1,2],mean_sd[,2]),type="S",col=co[2],lwd=2)
  abline(v=cp,lty=3,lwd=0.8,col=co[4])
  
  # plot 3: variance
  par(mar=c(0.5,4.5,0.5,0.5))
  plot(c(1,sort(cp),length(x)),c(mean_sd[1,3],mean_sd[,3]),type="S",axes=FALSE,ylab=expression(paste (hat(sigma))),col=co[3],col.lab=co[3],lwd=2); 
  axis(2,las=1)
  abline(v=cp,lty=3,lwd=0.8,col=co[4])
  
}#end-plot.mscp
