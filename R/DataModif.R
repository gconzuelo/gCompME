# g-Computation for Measurement Error
# Gabriel Conzuelo
# Data modification v.1

# This function depends on Erica Moodie Data and 
#  will not work on other data
DataModif <- function(dat,size,k,j,seed=seed) {
  require(tidyverse)
  expit <- function(x){ exp(x)/(1+exp(x)) }
  set.seed(seed)
  # Select obs. with gold std. (GS) for some time points
  idx <- sample(
    dat[dat$int==1,"ID"],
    nrow(subset(dat,int==1))*size
  )
  tmp <- dat %>%
    group_by(ID,int) %>%
    mutate(
      Vset = ifelse(
        ID %in% idx,1,0),
      gold = ifelse(
        Vset==0,0,rbinom(1,1,0.5)),
      mXpr = expit(log(k) + log(j)*tvX - time),
      # mXpr = expit(log(k)*tvX - log(j)*tvX - time),
      mX = rbinom(1,1,mXpr),
      gX = ifelse(
        gold==1,tvX,NA)
    ) %>%
    select(
      -c(mXpr)
    ) %>%
    ungroup()
  return(tmp)
}


