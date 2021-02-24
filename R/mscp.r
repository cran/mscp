#' mscp
#'
#' Multiscale change point detection via gradual bandwidth adjustment in moving sum processes. A method for the detection of changes in the expectation in univariate sequences.
#'
#' @param x numeric vector. Input sequence of random variables.
#' @param delta integer >=2. Default = 20. Minimal window considered.
#' @param g integer >=1. Default = 20. Spacing between starting points.
#' @param kappa NA or positive real number. Default = NA. Breaking threshold. If NA, then kappa is derived in simulations, using alpha and sim
#' @param alpha numeric in (0,1). Default = 0.01. Significance level, i.e., sets kappa as (1-alpha)-quantile of maximum of Gaussian process limit.  
#' @param sim integer >=1. Default = 500. Number of simulations for kappa. 
#' 
#' @return invisible list
#' \item{cp}{detected change points (ordered according to detection)}
#' \item{mean_sd}{matrix of estimated means and standard deviations}
#' \item{path}{list containing matrices, each matrix describing the path of a detected change point. First column: t-value, second column: h-value, third column: D-value (statistic), first row: starting values, last row: end values}
#' \item{S}{matrix of possible starting values. First column: t-value, second column: h-value, third column: D-value (statistic), fourth column: step when cut out}
#' \item{x}{input sequence}
#' \item{delta}{minimal window size}
#' \item{g}{spacing between starting points}
#' \item{kappa}{threshold}
#' 
#' @seealso \code{\link{plot.mscp}, \link{summary.mscp}}
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
#' @rdname mscp
#' @import stats
#' @import grDevices
#' @import graphics
#' @export
###### end

mscp <-function(x, delta=20, g=20, kappa=NA, alpha=0.01, sim=500){
  
  if(!is.vector(x) | !is.numeric(x)){stop("invalid choice of x.")}
  
  if((delta %% 1)!=0 | delta<2 | length(delta)>1 | !is.numeric(delta)){stop("invalid choice of delta.")}
  
  if((g %% 1)!=0 | g<1 | length(g)>1 | !is.numeric(g)){stop("invalid choice of g.")} 
  
  if(delta + ceiling(g/2) >= floor(length(x)/2)){stop("invalid paramters. Please decrease delta or g.")}
  
  if(!( is.na(kappa) | (kappa>0 & is.numeric(kappa) & length(kappa)==1))){stop("invalid choice of kappa.")}
  
  if(is.na(kappa) & (alpha*(1-alpha) <= 0)){stop("invalid choice alpha.")}  
  
  if(is.na(kappa) & (sim <= 0 | !is.numeric(sim) | length(sim) != 1)){stop("invalid choice of sim")}
  
  
  ###
  ### auxiliary functions
  ###
  
  # Compute D_t,h
  D_marginal <- function(t,x,h)
  {
    xr <- x[(t+1):(t+h)]; xl <- x[(t-h+1):t]
    denominator <- sqrt(sd(xr)^2 + sd(xl)^2)
    ifelse(denominator>0, (h^(1/2))*(mean(xr) - mean(xl))/denominator, 0)
  }#end-D_marginal   
  
  
  # Compute zigzag-path
  zigzagdown <- function(start,x=x,delta=delta)
  {
    xtmp <- start[1]; ytmp <- start[2]
    step <- c(-1,0,1)
    if(ytmp == xtmp)          {step <- c(0,1)} # on left leg avoid jump left
    if(ytmp == length(x)-xtmp){step <- c(-1,0)}# on right leg avoid jump right
    if(ytmp == xtmp & ytmp == length(x)-xtmp){step <- 0}# on upper edge move down
    Dsub   <- apply(as.matrix(xtmp+step),MARGIN=1,FUN=D_marginal,x=x,h=ytmp)
    choice <- which.max(abs(Dsub)) 
    xtmp   <- xtmp + step[choice]
    xpath <- c(xtmp); ypath <- c(ytmp)
    Dpath <- c(D_marginal(x=x,t=xtmp,h=ytmp))
    while(ytmp >= 1+delta)  
    {
      # 1. step-down (minimal stepsize)
      ytmp  <- ytmp - 1
      ypath <- c(ypath,ytmp)
      # 2. step-horizontal (go to max in horizontal neighborhood of distance mininal stepsize)
      tsub   <- xtmp + c(-1,0,1)
      step   <- c(-1,0,1)
      Dsub   <- apply(as.matrix(xtmp+step),MARGIN=1,FUN=D_marginal,x=x,h=ytmp)
      choice <- which.max(abs(Dsub)) # in {1,2,3}
      xtmp   <- xtmp + step[choice]  # move left , stay, or move right
      xpath  <- c(xpath,xtmp)
      Dpath  <- c(Dpath,Dsub[choice])
    }#end-while
    # return path and estimate (end of the path)
    path            <- cbind(xpath,ypath,Dpath)  
    colnames(path)  <- c("t","h","D")
    cpest           <- (path[,1])[length(path[,1])]
    list(cpest=cpest,path=path)
  }#-end-zigzagdown
  
  # Compute starting points
  startingpoints <- function(x=x,delta=delta,g=g){
    S <- c()
    for( h in delta:(length(x)/2)){
      for( t in h:(length(x)-h)){
        #if(t%%floor(g/2) == 0 & h%%floor(g/2)==0)
        if(t%%g == 0 & h%%g==0)
        {
          S <- rbind(S,c(t,h,(h^(-1/2))*D_marginal(t,x=x,h=h)))
        }
      }
    }
    S
  }#end-startingpoints
  
  # Simulate Brownian motion and calculate max(abs(L)) over triangle     
  sim_global_max_of_absL <- function(Tt=Tt,delta=delta)
  {
    W         <- c(0,cumsum(rnorm(Tt)))
    maxabsLh  <- c()
    for(h in delta:(floor(Tt/2)))
    {
      maxabsLht <- max(abs(W[(2*h+1):(Tt+1)] - 2*W[(1+h):(Tt-h+1)] + W[1:(Tt-2*h+1)])/sqrt(2*h)) 
      maxabsLh  <- c(maxabsLh,maxabsLht) 
    }
    max(maxabsLh)
  }#end-sim_global_max_of_absL  
  
  ###
  ### detect change points using auxiliary functions
  ###
  
  Tt <- length(x)
  
  # simulate kappa
  if(is.na(kappa))
  {
    kappa <- quantile(replicate(n=sim,sim_global_max_of_absL(Tt=Tt,delta=delta)),1-alpha)
  }  
  
  # calculate and order starting points
  thD         <- startingpoints(x=x,delta=delta,g=g)
  thD_ordered <- thD[rev(order(abs(thD[,3]))),] 
  n           <- dim(thD)[1]
  
  # placeholder to be updated
  cp              <- c()
  path            <- list()
  considered      <- rep(1,n)
  step            <- rep(Inf,n)
  step            <- rep(0,n)
  rowsconsidered  <- (1:n)
  mcounter        <- 1
  
  # successively detect change points
  while(length(thD_ordered[rowsconsidered,])>0)
  {
    zzd      <- zigzagdown(start=matrix(thD_ordered[rowsconsidered,],ncol=3)[1,1:2],x=x,delta=delta)
    cpcand   <- zzd$cpest # change point candidate
    pathcand <- zzd$path  # path of candidate
    # check three criteria
    # first criterion:
    Ddist <- ifelse(length(cp)>0,min(abs(cpcand-cp)),Inf) 
    go_up <- ifelse(Ddist<=2*(delta-1),TRUE,FALSE)
    if(go_up == TRUE)
    { 
      # cut out cone but do not accept candidate
      a1 <- -cpcand; b1 <- 1
      a2 <- +cpcand; b2 <- -1
      ttmp    <- thD_ordered[,1][considered==1]
      htmp    <- thD_ordered[,2][considered==1]
      keep    <- ifelse(htmp < a1 + b1*ttmp | htmp <=  a2 + b2*ttmp, 1, 0)
      considered[(1:n)[considered==1]] <- keep
      rowsconsidered        <- (1:n)[as.logical(considered)]
    }#end-if-go_up == TRUE
    if(go_up == FALSE)
    { 
      # second criterion
      if(max(abs(pathcand[,3])) < kappa){break} 
      # third criterion
      if(Ddist < g - 2*(delta-1)){break} 
      # accept candidate and cut out cone 
      cp                <- c(cp,cpcand)
      path[[mcounter]]  <- pathcand
      a1 <- -cpcand; b1 <- 1
      a2 <- +cpcand; b2 <- -1
      ttmp    <- thD_ordered[,1][considered==1]
      htmp    <- thD_ordered[,2][considered==1]
      keep    <- ifelse(htmp < a1 + b1*ttmp | htmp <=  a2 + b2*ttmp, 1, 0)
      considered[(1:n)[considered==1]] <- keep
      rowsconsidered        <- (1:n)[as.logical(considered)]
      #step[-rowsconsidered]  <- mcounter
      step[-rowsconsidered]  <- step[-rowsconsidered] +1 
      mcounter <- mcounter+1
    }#end-if-go_up == FALSE   
  }#end-while
  
  mstep           <- max(step)
  step            <- ifelse(step ==  0,-Inf,step)
  S               <- cbind(thD_ordered,1+mstep-step)
  colnames(S)     <- c("t","h","statistic","deleted after step") 
  if(length(cp)>0){names(cp) <- rep("",length(cp))}
  
  # means and standard deviations
  CPs       <- c(1,sort(cp),length(x)) 
  section   <- 1:(length(CPs)-1)
  hat_mean  <- rep(NA,(length(CPs)-1))        
  hat_sd    <- rep(NA,(length(CPs)-1))
  for(i in 1:(length(CPs)-1)){      
    hat_mean[i] <- mean(x[CPs[i]:CPs[i+1]])  
    hat_sd[i]   <-  sd(x[CPs[i]:CPs[i+1]])   
  }# end-for-i
  mean_sd <- cbind(section,hat_mean,hat_sd)
  
  # output
  print(sort(cp))
  lst <- list(cp=cp, mean_sd=mean_sd, path=path, S=S, x=x, delta=delta, g=g, kappa = kappa)
  ###
  ### 
  ###
  class(lst) <- "mscp"
  invisible(lst)
}#end-mscp 

