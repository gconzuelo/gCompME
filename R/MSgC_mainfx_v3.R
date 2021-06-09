# g-Computation for Measurement Error
# Gabriel Conzuelo
# Multi-step g-Computation (MSgC) v.1
source("~/Dropbox/MyR_funcs/saving.R")
library(survival)
library(tidyverse)
path88 <- list.files(path="./rawData/8080", pattern="*.txt", full.names=T)

# load data
df88 <- lapply(path88, function(x) read_delim(x, "\t", escape_double=F, trim_ws=T))

# Main function
# step 1 == F means no pre-modeling step 
dat=df88[[6]]; step1=T; seed=1
testMSgC <- function(dat,step1=T,seed=NULL) {
  require(tidyverse)
  require(splines)
  require(parallel)
  require(data.table)
  set.seed(seed)
  expit <- function(x){ exp(x)/(1+exp(x)) }
  pFunc <- function(mod,ndat){
    as.numeric(predict(mod, newdata=ndat, type="response")>runif(1))
  }
  cat("Now Running SEED",seed,'\n')
  d <- dat
  sens <- round(with(dat, sum(tX==1 & mX==1) / sum(tX==1)),2)
  spec <- round(with(dat, sum(tX==0 & mX==0) / sum(tX==0)),2)
  size <- round(with(dat, mean(Vset)),2)
  
  # step 1: model for ME correction in validation set,
  #  subset at time with available GS measure
  
  fitXprime <- glm(
    tX~mX, 
    data = subset(d, Vset==1), 
    family = binomial(link='logit')
  )
  #+bs(time, df=3)
  # fit_Xhat <- glm(
  #   tX~mX+bs(time, df=3), 
  #   data = subset(d, Vset==1), 
  #   family = binomial(link='logit')
  # )
  # fit_mX <- glm(
  #   mX~bs(time, df=3), 
  #   data = subset(d, Vset==1), 
  #   family = binomial(link='logit')
  # )
  
  # step 2: predict X at missing time points  
  #  for those in validation set. After this
  #  step, validation set has complete info 
  
  
  # vDF0 <- d %>%
  #   filter(Vset==1) %>% 
  #   select(ID,time,mX)
  # vDF0$id <- seq(1,nrow(vDF0),1)
  # pgf0 <- function(jj,dat,L) {
  #   set.seed(jj)
  #   d <- dat
  #   xhat <- mXp <- mm <- numeric()
  #   id <- d$id
  #   mm[1] <- time <- 1
  #   if (is.null(exposure)) {
  #     Xp <- d$X
  #   } else {
  #     Xp <- exposure
  #   }
  #   Zp[1] <- d$Z
  #   Yp[1] <- 0
  #   
  # }
  
  
  
  
  
  ####
  if (step1==T) {
    vDF <- d %>% 
      group_by(time) %>%
      filter(Vset==1) %>%
      mutate(
        xhat = ifelse(
          is.na(gX),pFunc(fitXprime,.),tX)
      )
  } else {
    vDF <- d %>% 
      group_by(ID) %>%
      filter(Vset==1) %>%
      mutate(
        xhat = gX
      )
  }

  #  Bootstrap step (optional)
  bootFx <- function(dat,seed) {
    boot <- NULL
    nms <- as.numeric(names(table(dat$ID)))
    idx <- sample(1:length(nms), length(nms), replace=T)
    ref <- table(nms[idx])
    if(is.null(seed)){
      boot <- dat
    } else {
      for(zzz in 1:max(ref)){
        cc <- dat[dat$ID %in% names(ref[ref %in% c(zzz:max(ref))]),] 
        cc$bid <- paste0(cc$ID,zzz)
        boot <- rbind(boot, cc)
      }}
    return(boot)
  }
  boot <- bootFx(vDF,seed)
  boot <- boot %>%
    rename(X = xhat) %>%
    mutate(
      Xm1 = lag(X,n=1L,default=0),
      Zm1 = lag(Z,n=1L,default=0)
    )
  with(subset(boot, gold==0), table(X, tX))
  with(subset(boot, gold==1), table(X, tX))
  with(subset(boot, Vset==1), cor.test(X, tX))
  
  # step 3: make outcome, time-varying exposure
  #  and confounding models.
  #  (1) outcome
  mY <- function(k){
    fitY <- glm(
      Y~X+Z+Xm1+Zm1+bs(time,df=3), 
      data=boot, 
      family=binomial(link="logit"))
    return(fitY)
  }
  #  (2) exposure
  mX <- function(k){
    fitX <- glm(
      X~Z+Xm1+Zm1+bs(time,df=3),
      data=boot,
      family=binomial(link="logit"))
    return(fitX)
  }
  #  (3) Time-varying confounder (Z)
  mZ <- function(k){
    fitZ <- glm(
      Z~Xm1+Zm1+bs(time,df=3), 
      data=boot, 
      family=binomial(link="logit"))
    return(fitZ)
  }
  #  (4) model to impute X
  mXc <- function(k){
    fitXc <- glm(
      X~mX+bs(time,df=3), 
      data=boot, 
      family=binomial(link="logit"))
    return(fitXc)
  }
  mL<-c(mZ,mX,mY,mXc)
  fitR<-lapply(1:4,function(x) mL[[x]](k))

  # step 4: make baseline data for pgf function
  #  Here, we predict X for everyone outside validsation set
  cDF <- dat %>% 
    filter(Vset==0) %>%
    mutate(xhat=NA)
  DF <- rbind(data.frame(cDF), data.frame(vDF))
  DF <- DF %>%
    filter(int==1) %>%
    group_by(ID) %>%
    arrange(ID) %>%
    mutate(
      xhat = ifelse(
        is.na(xhat),pFunc(fitR[[4]],.),xhat)
    ) %>%
    rename(X = xhat)
  DF$id <- seq(1,nrow(DF),1)
  # pgf
  pgf <- function(ii, dat, lngth, exposure=NULL){
    set.seed(ii)
    d <- dat
    d <- d[d$id==ii,]
    # cat("Now running observation",ii,"from Monte Carlo Data",'\n')
    lngth <- lngth
    Zp <- Xp <- Yp <- mm <- numeric()
    id <- d$id
    mm[1] <- time <- 1
    if (is.null(exposure)) {
      Xp <- d$X
    } else {
      Xp <- exposure
    }
    Zp[1] <- d$Z
    Yp[1] <- 0
    
    for (time in 2:lngth) {
      if (Yp[time-1]==0) {
        Xm1=Xp[time - 1]
        Zm1=Zp[time - 1]
        # Generating time-varying confounder 
        dZp <- data.table(Xm1,Zm1,time)
        Zp[time] <- pFunc(fitR[[1]], dZp)
        # Generating exposure
        if (is.null(exposure)) {
          dXp <- data.table(Z=Zp[time],Xm1,Zm1,time)
          Xp[time] <- pFunc(fitR[[2]], dXp)
        } else {
          Xp[time] <- exposure
        }
        # Generating outcome
        dYp <- data.table(X=Xp[time],Z=Zp[time],Xm1,Zm1,time)
        Yp[time] <- pFunc(fitR[[3]], dYp)
        
      } else {
        break
      }
      mm[time] <- time
    }
    boot_num<-seed
    gdat <- data.table(boot_num,id,mm,Zp,Zm1,Xp,Xm1,Yp)
    gdat$last<-as.numeric(!(gdat$Yp==0)|gdat$mm==lngth)
    return(gdat)
  }
  cores <- detectCores()
  res <- mclapply(1:nrow(DF), function(ii) pgf(ii,DF,12,exposure=NULL), mc.cores=cores)
  tmp <- do.call(rbind,res)
  tmp$sens <- sens
  tmp$spec <- spec
  tmp$size <- size

  cat("End of run for SEED",seed,'\n')
  return(tmp)
}

# run
library(pbmcapply)
nsim <- 200
cores <- detectCores()
r <- list()
dL <- df88 # choose your data list
for (i in 1:length(dL)) {
  r[[i]] <- pbmclapply(
    1:nsim, 
    function(x) testMSgC(dL[[i]],step1=T,x),
    mc.cores=cores,
    mc.set.seed=T)
  re <- do.call(rbind,r[[i]])
  nm <- paste0("MSgCres_",Sys.Date(),"_dupCode",sample(1:999,1))
  saveData(re, nm)
}







1# debug
r <- lapply(1:2, function(x) testMSgC(dfL[[1]],x))

