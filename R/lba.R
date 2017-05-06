
#' @rdname rplba
#' @export
rplbaR1 <- function(n=10, pVec=c(A=1.51, b=2.7, v1=3.32, v2=2.24, w1=1.51, w2=3.69, 
                          sv=1, rD=0.31, swt=0.5, t0=0.08)) 
{
  ##  A  b  v1   v2  w1   w2  sv  rD  swt t0
  ##  0  1   2    3   4    5   6   7    8  9 
  T0 <- pVec["swt"] + pVec["rD"]; T0

  ## Stage 1 LBA
  v1 <- rtnorm(n, pVec["v1"], pVec["sv"],0); 
  v2 <- rtnorm(n, pVec["v2"], pVec["sv"],0); 
  sp <- matrix(runif(2*n, 0, pVec["A"]), nrow=2)
  
  ## Race
  dt <- rbind((pVec["b"]-sp[1,])/v1, (pVec["b"]-sp[2,])/v2)
  
  ## dt[dt<0] <- Inf
  choice <- apply(dt, 2, which.min)
  rt <- dt[cbind(choice, 1:n)]  ## convert to vector choose (row, col)
  
  ## Which are finished?
  done <- rt <= T0
  n2 <- sum(!done)
  
  ## Distance left to travel for those not finished
  B1 <- pVec["b"] - (sp[1, !done] + T0*v1[!done]) 
  B2 <- pVec["b"] - (sp[2, !done] + T0*v2[!done]) 
  
  ## Stage 2 LBA
  w1 <- msm::rtnorm(n2, pVec["w1"], pVec["sv"], 0)
  w2 <- msm::rtnorm(n2, pVec["w2"], pVec["sv"], 0)
  
  ## Race  
  dt <- rbind(B1/w1, B2/w2)
  choice[!done] <- apply(dt, 2, which.min)
  rt[!done] <- T0 + dt[cbind(choice[!done], 1:n2)]
  
  ## save results
  cbind(choice=choice,rt=pVec["t0"]+rt)
  
}

#' @rdname rplba
#' @export
rplbaR2 <- function(n=10, pVec=c(A1=1.51, A2=1.51, b1=2.7, b2=2.7, v1=3.32, v2=2.24, 
                          w1=1.51, w2=3.69, sv1=1, sv2=1, sw1=1, 
                          sw2=1,  rD=0.31,  swt=0.5, t0=0.08)) {
  ## A1  A2  b1  b2  v1  v2  w1  w2  sv1  sv2  sw1  sw2  rD swt   t0
  ##  1   2   3   4   5   6   7   8    9   10   11   12  13  14   15

  # Stage 1 LBA
  v1 <- rtnorm(n,pVec["v1"], pVec["sv1"],0)
  v2 <- rtnorm(n,pVec["v2"], pVec["sv2"],0)
  sp <- matrix(runif(2*n, 0, pVec[c("A1","A2")]), nrow=2)
  
  # Race
  dt <- rbind((pVec[c("b1","b2")]-sp[1,])/v1,(pVec[c("b1","b2")]-sp[2,])/v2)
  # dt[dt<0] <- Inf
  choice <- apply(dt,2,which.min)
  rt <- dt[cbind(choice,1:n)]
  
  # Calculate effective switch time
  swt <- pVec["swt"] + pVec["rD"]
  
  # Which are finished?
  done <- rt <= swt
  n2 <- sum(!done)
  
  # Distance left to travel for those not finished
  B1 <- pVec["b1"] - (sp[1,!done] + swt*v1[!done]) 
  B2 <- pVec["b2"] - (sp[2,!done] + swt*v2[!done]) 
  
  # Stage 2 LBA
  w1 <- rtnorm(n2,pVec["w1"],pVec["sw1"],0)
  w2 <- rtnorm(n2,pVec["w2"],pVec["sw2"],0)
  
  # Race  
  dt <- rbind(B1/w1,B2/w2)
  choice[!done] <- apply(dt,2,which.min)
  rt[!done] <- swt+dt[cbind(choice[!done],1:n2)]
  
  # save results
  cbind(choice=choice,rt=pVec["t0"]+rt)
  
}

#' @rdname rplba
#' @export
rplbaR <- function(n, pVec=c(A1=1.5, A2=1.5, B1=1.2, B2=1.2, C1=.3, C2=.3,
                              v1=3.32, v2=2.24, w1=1.51, w2=3.69, sv1=1, sv2=1,
                              sw1=1, sw2=1, rD=0.3, tD=.3, swt=0.5, t0=0.08)) {

  # Stage 1 LBA
  v1 <- rtnorm(n, pVec["v1"], pVec["sv1"], 0)
  v2 <- rtnorm(n, pVec["v2"], pVec["sv2"], 0)
  sp <- matrix(runif(2*n,0,pVec[c("A1","A2")]),nrow=2)
  
  # Calcualte thresholds
  b1 <- sum(pVec[c("A1","B1")])
  b2 <- sum(pVec[c("A2","B2")])
  c1 <- b1 + pVec[c("C1")]
  c2 <- b2 + pVec[c("C2")]
  
  # Race
  dt <- rbind((c(b1,b2)-sp[1,])/v1,(c(b1,b2)-sp[2,])/v2)
  # dt[dt<0] <- Inf
  choice <- apply(dt,2,which.min)
  rt <- dt[cbind(choice,1:n)]
  
  # Calculate effective switch times
  swt_b <- pVec["swt"] + pVec["tD"]
  swt_r <- pVec["swt"] + pVec["rD"]
  
  # Which switch is first
  swt <- pmin(swt_b,swt_r)
  if (swt_b==swt_r) {
    change <- "both"
  } else if (swt_r < swt_b) {
    change <- "rate"
  } else { 
    change <- "threshold"
  }
  
  

  # Which are finished?
  done <- rt <= swt
  n2 <- sum(!done)
  
  # Stage 2 LBA
  
  # Distance left to travel for those not finished
  # threshold - distance already travelled
  if ( change=="rate" ) {
    B1 <- b1 - (sp[1,!done] + swt*v1[!done]) 
    B2 <- b2 - (sp[2,!done] + swt*v2[!done]) 
  } else {
    B1 <- c1 - (sp[1,!done] + swt*v1[!done]) 
    B2 <- c2 - (sp[2,!done] + swt*v2[!done]) 
  }
  
  # Change rates?
  if ( change=="threshold" ) {
    w1 <- v1[!done]; w2 <- v2[!done]
  } else {
    w1 <- rtnorm(n2,pVec["w1"],pVec["sw1"],0)
    w2 <- rtnorm(n2,pVec["w2"],pVec["sw2"],0)
  } 
  
  # Race  
  dt <- rbind(B1/w1,B2/w2)
  # dt[dt<0] <- Inf
  choice[!done] <- apply(dt,2,which.min)
  rt[!done] <- swt+dt[cbind(choice[!done],1:n2)]
  
  if ( change != "both" ) { # Stage 3 LBA
    
    if ( change=="threshold" ) swt1 <- swt_r else swt1 <- swt_b
    t2 <- swt1-swt
    
    # Which are finished?
    done1 <- rt[!done] < swt1
    n2 <- sum(!done1)
    
    if ( !all(done1) ) {
      
      # Distance left to travel for those not finished
      # Distance left at end of stage 1 - further travel
      B1 <- B1[!done1] - t2*w1[!done1]
      B2 <- B2[!done1] - t2*w2[!done1]
      
      if ( change=="threshold" ) {
        w1 <- rtnorm(n2,pVec["w1"],pVec["sw1"],0)
        w2 <- rtnorm(n2,pVec["w2"],pVec["sw2"],0)
      }  else {
        w1 <- w1[!done1]; w2 <- w2[!done1]
        B1 <- B1 + pVec["C1"]
        B2 <- B2 + pVec["C2"]
      }
      
      # Race  
      dt <- rbind(B1/w1,B2/w2)
      # dt[dt<0] <- Inf
      choice[!done][!done1] <- apply(dt,2,which.min)
      rt[!done][!done1] <- swt1+dt[cbind(choice[!done][!done1],1:n2)]
    }
    
  }
  
  # save results
  cbind(choice=choice,rt=pVec["t0"]+rt)
  
  
}



#' @export
rlba_test <- function(n, b=1, A=.5, mean_v=c(2.4, 1.6), sd_v=c(1, 1), t0=.5) {
  out <- data.frame(rlba_internal(n, b, A, mean_v, sd_v, t0))
  names(out) <- c("RT", "R")
  return(out)
}