# Template setup for n-choice LBA, B=b-A parameterization
#   External parameters types: A, B, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

#source("rtdists_extras.R")

# This function transfroms parameters to a form suitbale for the model 
#   being used. Called inside of get.par.mat. 
# "par.df" is a data frame of parameters types , some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used ("norm" etc.) is present.
transform.dmc <- function(par.df) 
{
  # User supplied tranforms go here
  par.df$b <- par.df$B+par.df$A
  
#   # COMMENT OUT this check for speed after debugging
#   if ( !all(type.par.names %in% names(par.df)) )
#     stop("Trasform has not created parameter(s) required by the model.")
  
  par.df[,c("A","b","t0","mean_v","sd_v","st0")]
}


random.dmc <- function(n,p.df,model)
{
  rlba.norm(n,A=p.df$A,b=p.df$b,t0=p.df$t0, 
                               mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
                               posdrift = attr(model,"posdrift"))
}


likelihood.dmc <- function(p.vector,data,ok.types=c("norm"),min.like=1e-10)   
  #We can probably use more ok types to vary the rate parameters and have all LBA rate models in each LBA model 
  #I.e. all lba_B models together or all lba_B1V models together, etc.
# Returns vector of likelihoods for each RT in data (in same order)
# !!! TO DO: types other than norm
{

  common.cell <- function(cell.response) {
    # Gets indexs of model row.names that share cell
    cells <- unlist(lapply(
      strsplit(cell.response,".",fixed=TRUE),function(x){x[[-length(x)]]}))
    list(use=!duplicated(cells),index=as.numeric(factor(cells))) 
  }

  
  likelihood <- numeric(dim(data)[1])
  ui <- common.cell(row.names(attr(data,"model")))
  use.cells <- row.names(attr(data,"model"))[ui$use]
  
  for ( i in 1:length(use.cells) ) { 
    
    cells.to.fill <- row.names(attr(data,"model"))[ui$index==i]
    rts <- sapply(cells.to.fill,function(x){
      data$RT[attr(data,"cell.index")[[x]]]  
    })
    p.df <- p.df.dmc(p.vector,use.cells[i],attributes(data)$model,n1order=TRUE)

    samp <- try(rlba.norm(attr(data,"n.pda"), 
      A=p.df$A, b=p.df$b, t0=p.df$t0, mean_v=p.df$mean_v,
      sd_v=p.df$sd_v, st0=p.df$st0[1],
      posdrift = attr(attr(data,"model"),"posdrift")),silent=TRUE)

    samp <- samp[samp$rt<10,]
    
    for (j in 1:length(cells.to.fill))  if ( length(rts[[j]]) > 0 ) {
      if ((class(samp)!="try-error")) {
        is.samp <- samp$response==j
        if (any(is.samp)) {
       if ( length(rts[[j]]) < 0 ) 
         likelihood[ attr(data,"cell.index")[[ cells.to.fill[j]]] ] <- 
           cpda::logLik_pw(rts[[j]], samp[is.samp,"rt"], n=attr(data,"n.pda")) else
         likelihood[ attr(data,"cell.index")[[ cells.to.fill[j]]] ] <- 
           cpda::logLik_fft2(rts[[j]], samp[is.samp,"rt"], n=attr(data,"n.pda"))[,2]
        } else likelihood[ attr(data,"cell.index")[[ cells.to.fill[j]]] ] <- 0
      } else likelihood[ attr(data,"cell.index")[[ cells.to.fill[j]]] ] <- 0
    }
      
  }
  pmax(likelihood,min.like)
}



likelihood.dmc <- function(p.vector,data,ok.types=c("norm"),min.like=1e-10)   
  #We can probably use more ok types to vary the rate parameters and have all LBA rate models in each LBA model 
  #I.e. all lba_B models together or all lba_B1V models together, etc.
# Returns vector of likelihoods for each RT in data (in same order)
# !!! TO DO: types other than norm
{

  common.cell <- function(cell.response) {
    # Gets indexs of model row.names that share cell
    cells <- unlist(lapply(
      strsplit(cell.response,".",fixed=TRUE),function(x){x[[-length(x)]]}))
    list(use=!duplicated(cells),index=as.numeric(factor(cells))) 
  }

  my.density <- function(rts,srt,n) {
    if (length(srt)<10) return(numeric(length(rts)))
    d <- density(srt)
    d$y[d$y<0] <- 0
    d$y <- d$y*length(srt)/n
    out <- numeric(length(rts))
    ok <- (rts>d$x[1]) & (rts<d$x[length(d$x)])
    out[ok] <- interp1(d$x,d$y,rts[ok])
    out[is.na(out) | !is.finite(out)] <- 0
    out
  }

  
  likelihood <- numeric(dim(data)[1])
  ui <- common.cell(row.names(attr(data,"model")))
  use.cells <- row.names(attr(data,"model"))[ui$use]
  
  for ( i in 1:length(use.cells) ) { 
    
    cells.to.fill <- row.names(attr(data,"model"))[ui$index==i]
    rts <- sapply(cells.to.fill,function(x){
      data$RT[attr(data,"cell.index")[[x]]]  
    })
    p.df <- p.df.dmc(p.vector,use.cells[i],attributes(data)$model,n1order=TRUE)

    samp <- try(rlba.norm(attr(data,"n.pda"), 
      A=p.df$A, b=p.df$b, t0=p.df$t0, mean_v=p.df$mean_v,
      sd_v=p.df$sd_v, st0=p.df$st0[1],
      posdrift = attr(attr(data,"model"),"posdrift")),silent=TRUE)

    if ((class(samp)!="try-error")) samp <- samp[samp$rt<10,]
    
    for (j in 1:length(cells.to.fill))  if ( length(rts[[j]]) > 0 ) {
      if ((class(samp)!="try-error")) {
        is.samp <- samp$response==j
        if (any(is.samp)) {
         likelihood[ attr(data,"cell.index")[[ cells.to.fill[j]]] ] <- 
           my.density(rts[[j]], samp[is.samp,"rt"], n=attr(data,"n.pda")) 
        } else likelihood[ attr(data,"cell.index")[[ cells.to.fill[j]]] ] <- 0
      } else likelihood[ attr(data,"cell.index")[[ cells.to.fill[j]]] ] <- 0
    }
      
  }
  pmax(likelihood,min.like)
}


# # Check density and normalization
# n=1e5; A=.25; B=.35; mean_v=c(1,.25); sd_v=1; t0=0; b=A+B
# sim <- rlba.norm(n=n,A=A,b=b,mean_v=mean_v,sd_v=sd_v,t0=t0,st0=0)
# names(sim) <- c("RT","R")
# par(mfrow=c(1,2))
# dns <- plot.cell.density(sim,xlim=c(0,4),save.density=TRUE)
# dt=dns$"1"$x
# d <- n1PDF(dt,A=A,b=b,mean_v=mean_v,sd_v=sd_v,t0=t0,st0=0)
# # red=black?
# plot(dns$"1"$x,dns$"1"$y,type="l")
# lines(dns$"1"$x,d,col="red")
# 
# rt.c <- sim$RT[sim$R==1]
# rt.c <- rt.c[rt.c<4]
# pc <- mean(sim$R==1)
# system.time({d.c <- density(rt.c)}); d.c$y <- d.c$y*pc
# lines(d.c$x,d.c$y,col="green")
# 
# system.time({d.cpda <- cpda::logLik_pw(d.c$x, rt.c, n=n)})
# lines(d.c$x,d.cpda,col="blue")
# 
# system.time({d.r <- my.density(d.c$x, rt.c, n=n)})
# lines(d.c$x,d.r,col="orange")
