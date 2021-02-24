#' summary.mscp
#'
#' Summary method for class 'mscp'
#'
#' @param object object of class mscp
#' @param ... additional arguments
#'
#' @return No return value, called for side effects
#'
#' @seealso \code{\link{mscp}, \link{plot.mscp}}
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
#' @rdname summary.mscp
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end


summary.mscp <- function(object,...)
{
  cat("","\n")
  cat("--------------------------------------------------------------")
  cat("","\n")
  cat("Multiscale change point detection")   
  cat("","\n")
  cat("","\n")
  cat("Input: ");
  cat("","\n"); cat("         ")
  cat(paste("Sequence of length = ",length(object$x),sep="")); 
  cat("","\n"); cat("         ")
  cat("","\n"); 
  
  cat("Change point detection: "); 
  cat("","\n"); cat("         ")
  if(length(object$cp)==0){cat("No change points")}
  if(length(object$cp)==1){cat(paste(length(object$cp),"change point: "))} 
  if(length(object$cp)>1){cat(paste(length(object$cp),"change points (time ordered): "))}
  if(length(object$cp)>0){cat(paste(sort(object$cp)),sep=", ")}; 
  if(length(object$cp)>1){cat("","\n"); cat("         "); cat(paste("order of detection: ")); cat(paste(order(object$cp)),sep=", ");}
  cat("","\n");cat("","\n");
  cat("Parameter estimation: "); 
  cat("","\n"); cat("         ")
  if(length(object$mean_sd[,2])==1){cat(paste(length(object$mean_sd[,2]),"section with"))}
  if(length(object$mean_sd[,2])>1){cat(paste(length(object$mean_sd[,2]),"sections (time ordered) with"))}
  cat("","\n"); cat("              ")
  if(length(object$mean_sd[,2])==1){cat(paste("estimated expectation: ",signif(object$mean_sd[,2],2)) )} else{cat("estimated expectations: ")
    cat(paste(signif(object$mean_sd[,2],2)),sep=", ")}; 
  cat("","\n"); cat("              ")
  if(length(object$mean_sd[,3])==1){cat(paste("estimated standard deviation: ",signif(object$mean_sd[,3],2)) )} else{cat("estimated standard deviations: "); cat(paste(signif(object$mean_sd[,3],2)),sep=", ")}
  cat("","\n");#cat("","\n");
  cat("--------------------------------------------------------------")
}# end-summary.mscp
