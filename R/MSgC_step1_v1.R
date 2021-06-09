MSgC_s1s2 <- function(dat,seed=NULL) {
  require(splines)
  set.seed(seed)
  expit <- function(x){ exp(x)/(1+exp(x)) }
  pFunc <- function(mod,ndat){
    as.numeric(predict(mod, newdata=ndat, type="response")>runif(1))
  }
  # load data
  d <- dat
  # model to predict tX
  fitX_s1 <- glm(
    gX~mX+bs(time, df=3), 
    data = subset(d, gold==1), 
    family = binomial(link='logit')
  )
  # make validation DF0 (raw)
  vDF_s1 <- d %>%
    group_by(ID) %>%
    filter(Vset==1) %>%
    mutate(
      id_Vset = cur_group_id(),
      mXl = lag(mX,n=1L,default=0)
    ) %>%
    select(
      ID,id_Vset,int,time,tX,gX,mX,mXl
    ) %>%
    ungroup()
  # mini pgf
  mini_pgf <- function(dat, ii) {
    set.seed(ii)
    # bring data from observation ii
    d <- dat[dat$id_Vset==ii,]
    tXp <- mXp <- int <- numeric()
    # predict tX at time 1
    int[1] <- time <- 1
    mXp[1] <- d$mX[time]
    tXp[1] <- pFunc(fitX_s1, d[time,])
    id <- d$id_Vset[time]
    # predict tX at time 2 and onwards
    if (max(d$int) > 1) {
      for (time in 2:max(d$int)) {
        mXl <- mXp[time - 1]
        mXp[time] <- d$mX[time]
        dXp <- data.frame(mX=mXp[time],mXl=mXl,time)
        tXp[time] <- pFunc(fitX_s1, dXp)
        int[time] <- time
      }
      orig_ID<-d$ID
      gdat <- data.frame(ID=orig_ID,int=int,tXp=tXp)
    } else {
      orig_ID<-d$ID
      gdat <- data.frame(ID=orig_ID,int=int,tXp=tXp)
    }
    return(gdat)
  }
  r1 <- pblapply(1:max(vDF_s1$id_Vset), function(ii) mini_pgf(vDF_s1, ii))
  mPGFr1 <- data.frame(do.call(rbind, r1))
  # join results to validation set (vDF0)
  vDF <- left_join(vDF_s1, mPGFr1, by=c("ID","int"))
  vDF <- vDF %>%
    mutate(
      tXp = ifelse(is.na(gX),tXp,gX),
      sens_s0 = sum(tX==1 & mX==1) / sum(tX==1),
      spec_s0 = sum(tX==0 & mX==0) / sum(tX==0),
      sens_s1 = sum(tX==1 & tXp==1) / sum(tX==1),
      spec_s1 = sum(tX==0 & tXp==0) / sum(tX==0)
    ) %>%
    select(ID,int,mXl,tXp,starts_with(c("se","sp")))
  # join with original data
  d1 <- left_join(d,vDF,by=c("ID","int"))
  # end of step 1
  return(d1)
}

# run function
MSgC_s1s2(df[[6]], seed=123)

results <- list()
for (i in 1:length(df)) {
  results[[i]] <- MSgC_s1s2(df[[i]], seed=32145)
}
# save results
re <- list()
for (i in 1:length(results)) {
  saveData(results[[i]], name=paste0("r",sample(1:1e6,1)) )
}



