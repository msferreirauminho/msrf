require(extRemes)
require(boot)
require(evd)
require(exdex)
require(blocklength) # function pwsd() is based on the Automatic Block-Length selection 
# method proposed by Politis and White (2004) and corrected in Patton, Politis, and White (2009).


############################################################################################
############################################################################################
## Runs estimator and block bootstrap CI
############################################################################################

extInd_R <- function(x,threshold,k)
{
  xj <- x
  bj <- pwsd(xj, round = T, correlogram = F)$BlockLength[1]
  ffdir <- function(xj) {
    unj <- threshold
    kj <- k
    resffdir <- exi(xj, u=unj, r=kj)
    return(resffdir)
  }
  bootffdir<-tsboot(xj, ffdir, R = 500, l = bj, sim = "fixed")
  bci<-boot.ci(bootffdir,type = "perc")
  est <- bootffdir$t0
  li <- bci$percent[4]
  ls <- bci$percent[5]
return(c(est, li, ls))
}
############################################################################################

############################################################################################
############################################################################################
## Cycles estimator and block bootstrap CI
############################################################################################

extInd_C <- function(x, threshold, k)
{
  xj <- x
  bj <- pwsd(xj, round = T, correlogram = F)$BlockLength[1]
  ffdir <- function(xj) {
    k <- (k+1)
    unj <- threshold
    z <- numeric(0)
    lim <- floor(n/(k-1))
    for(i in 1:lim){
      z[i]=max(xj[((i-1)*(k-1)+1):(i*(k-1))])
    }
    excxj <- sum(xj>unj)
    fact <- sum(z>unj)/excxj
    lim <- length(z)
    uz <- sum(c(z[1:(lim-1)]<=unj & z[2:lim]>unj))
    return(uz/excxj)
  }
  bootffdir <- tsboot(xj, ffdir, R = 500, l = bj, sim = "fixed")
  bci <- boot.ci(bootffdir,type = "perc")
  est <- bootffdir$t0
  li <- bci$percent[4]
  ls <- bci$percent[5]
return(c(est, li, ls))
}
############################################################################################


############################################################################################
############################################################################################
## Truncated estimator and block bootstrap CI
############################################################################################

extInd_T <- function(x, threshold, k)
{
  xj <- x
  bj <- pwsd(xj, round = T, correlogram = F)$BlockLength[1]
  ffdir <- function(xj) {
    n <- length(xj)
    uj <- threshold
    D <- k
    eid <- xj > uj
    Nu <- sum(eid)
    inst <- (1:n)[eid]
    T <- diff(inst)
    s <- (T-D)
    pos.s <- (s>0)
    ND <- sum(pos.s)
    trf <- ND/sum(s[pos.s])
    Fu <- Nu/n
    trs <- trf/Fu
    trBC <- (((Nu-1)*trs-1)*(1/(Nu-1+D)))
    trj <- (trBC-(Fu/(2*(Nu-1)))*(1+trBC*(Nu-4)-(trBC^2*(Nu-1))))
    return(trj)
  }
  bootffdir <- tsboot(xj, ffdir, R = 500, l = bj, sim = "fixed")
  bci <- boot.ci(bootffdir,type = "perc")
  est <- bootffdir$t0
  li <- bci$percent[4]
  ls <- bci$percent[5]
return(c(est, li, ls))
}
############################################################################################


############################################################################################
############################################################################################
## K-gaps estimator and asymptotic normal CI (package exdex)
############################################################################################

extInd_K <- function(x, threshold, k)
{
  xj <- x
  unj <- threshold
  kj <- k
  reskgs <- kgaps(xj, unj,kj)
  ICi <- (confint(reskgs))
  est <- reskgs$theta
  li <- ICi[2]
  ls <- ICi[4]
return(c(est, li, ls))
}
############################################################################################


############################################################################################
############################################################################################
## Intervals estimator and bootstrap CI (Ferro & Segers 2003)
############################################################################################

extInd_I <- function(x, threshold)
{
  xj <- x
  unj <- threshold
  resfs<-extremalindex(xj, unj, method="intervals")
  ICi<-(ci(resfs)[c(1,7)])
  est <- resfs[1]
  li <- ICi[1]
  ls <- ICi[2]
return(c(est, li, ls))
}
############################################################################################



############################################################################################
############################################################################################
## Northrop estimator and asymptotic normal CI (package exdex)
############################################################################################

extInd_K <- function(x, b)
{
  xj <- x
  nor<-spm(xj,b=b)
  ic<-confint(nor,maxima="sliding",interval_type="norm")
  est <- nor$theta_sl[1]
  li <- ic$cis[1]
  ls <- ic$cis[4]
return(c(est, li, ls))
}
############################################################################################

############################################################################################
############################################################################################
## Ferreira & Ferreira sliding blocks estimator and percentile CI (Ferreira & Ferreira 2022)
############################################################################################

extInd_F <- function(x, b, alpha=0.05, re=100)
{
  xj <- (-1/(log(rank(x)/(n+1))))
  est <- numeric(0)
  for(l in 1:re){
  z <- (-1/log(runif(n))) 
  y1 <- z
  y2 <- (0.5*apply(cbind(z,xj),1,max))
  maxy1 <- 0
  maxy2 <- 0
  fim <- (n-b+1)
  for(i in 1:fim){
    maxy1[i]=max(y1[(i):(i+b-1)])
    maxy2[i]=max(y2[(i):(i+b-1)])
   }
  dimy <- length(maxy1)
  yn <- data.frame(maxy1,maxy2) # Y_n=(Y_1,Y_2)
  w <- cbind(rank(maxy1)/(dimy+1),rank(maxy2)/(dimy+1)) 
  maxy1y2 <- apply(w,1,max)
  ss <- (3-1/(1-mean(maxy1y2)))
  ssc <- max(0.5,ss)
  est[l] <- (1/ssc-1)
  }
  est <- mean(est)
  li <- quantile(est,alpha/2)
  ls <- quantile(est,1-alpha/2)
return(c(est, li, ls))
}
############################################################################################



############################################################################################
### IMT procedure to choose threshold and k in D^{k+1} (Fukutome et al. 2014, 2019)
############################################################################################

IMT <- function(x, quants=seq(0.8, .995, .005), run=1:12)
{
  nrun <- length(run)
  nthreshold <- length(quants)
  n <- length(x)
  trip <- numeric(0)
  for(j in 1:nthreshold)
    {
      u <- quantile(x,quants[j])
      for(i in 1:nrun)
      {
        k <- run[i]
        locExc <- which(x>u)
        N <- length(locExc)
        cl <- clusters(x,u,k,cmax=T)
        Nc <- length(cl)
        if(N>1)
        {
          Tis <- (locExc[2:N]-locExc[1:(N-1)])
          tabcis <- data.frame(Tis-k,rep(0,N-1))
          cis <- (N/n)*apply(tabcis,1,"max")
          
          loglik <- function(theta,cis)
          {
            ll <- ((N-1)-Nc)*log(1-theta)+2*Nc*log(theta)-theta*sum(cis)
            return(ll)
          }
          mloglik <- function(theta,cis)
          {
            ll <- (((N-1)-Nc)*log(1-theta)+2*Nc*log(theta)-theta*sum(cis))
            mll <- (-ll)
            return(mll)
          }
          thetaStart <- 0.5
          maximization <- nlminb(start = thetaStart, objective = mloglik,
                                 cis = cis, lower = c(0), upper = c(.99))
          estimate    <- maximization$par
          results      <- list(thetaEstimate = estimate,
                               logLikelihood = loglik(estimate, cis))
          cisnul<-ifelse(cis==0,1,0)
          cispos<-ifelse(cis>0,1,0)
          
          lilinha <- (-cisnul/(1-estimate)+2*cispos/estimate-cis)
          ii <- (cisnul/(1-estimate)^2+2*cispos/estimate^2)
          di <- (2*cispos/estimate^2+cis^2-4*cis/estimate)
          dilinha <- (-4*cispos/estimate^3+4*cis/estimate^2)
          D <- mean(di)
          I <- mean(ii)
          Dl <- mean(dilinha)
          sub2 <- (Dl*I^(-1)*lilinha)
          V <- mean((di-sub2)^2)
          IMT <- ((N-1)*D^2/V)
          if(is.nan(IMT)) trip <- rbind(c(0,0,0,0,0),trip)
          else
            if(IMT<0.05) trip <- rbind(c(u,k,estimate,Nc,N),trip) else 
              trip <- rbind(c(u,k,estimate,0,0),trip)
        }
        else trip <- rbind(c(0,0,0,0,0),trip)
      }
    }
    posmaxNc <- which.max(trip[,4])
    res <- trip[posmaxNc,]
    return(res)
}
  
############################################################################################