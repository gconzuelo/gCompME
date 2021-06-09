# g-Computation for Measurement Error
# Gabriel Conzuelo
# Generation of data with measuremnt error v.1
source("R/DataGen.R")
source("R/DataModif.R")
source("~/Dropbox/MyR_funcs/saving.R")
library(tidyverse)
library(survival)
library(data.table)
library(pbapply)
library(pbmcapply)

# This function uses Erica Moodie's function (DataGen.R) and 
#  modifies the result DF (DataModif.R) to have a mismeasured
#  exposure with a given sens/spec
DataGenME <- function(n,size,k,j, seed) {
  d <- DataGen(n=n,N=12,K=1,seed=seed)
  d <- d %>%
    rename(
      tvX = X,
      tvXl = Xm1,
      tvZ = Z,
      out = Y
    ) %>%
    select(
      -c(Zm1,timem1)
    )
  tmp <- DataModif(d,size,k=k,j=j,seed=seed)
  tmp <- tmp %>%
    rename(tX = tvX,
           Z = tvZ,
           Y = out)
  return(tmp)
}

# This function automatizes process to obtain multiple DFs with 
#  different errors and sizes of Validation Set
specNm <- function(x, n=1000, seed=32145) {
  size <- c(0.05,0.1,0.2,0.3)
  tmp <- crossing(n, size=size, k=x[["k"]], j=x[["j"]], seed) %>%
    group_by(size,j,k) %>%
    mutate(id=cur_group_id()) %>%
    group_split()
  tmp <- lapply(1:length(tmp), function(x) {tmp[[x]][,-6]} )
  return(tmp)
}

# run
# simulation specifications
sp60 <- list(k=72, j=c(15,40,160))
sp70 <- list(k=24, j=c(35,130,600))
sp80 <- list(k=8, j=c(105,400,1500))
spec <- specNm(sp80)   # pick 
df <- list()
for (i in 1:length(spec)) {
  df[[i]] <- do.call(DataGenME, spec[[i]])
}
# compute sens/spec and save results
df <- df %>%
  mutate(Vset_size = round(mean(Vset),2),
         sens_mX = round(sum(tX==1 & mX==1) / sum(tX==1),2),
         spec_mX = round(sum(tX==0 & mX==0) / sum(tX==0),2))
name <- paste0("rawData/dfme_",Sys.Date(),
               "_vs_",df$Vset_size[1],
               "_se_",df$sens_mX[1],
               "_sp_",df$spec_mX[1])
saveData(df, name)


######################## CHECKS ########################  
# test.df <- pika[[1]]
# head(test.df)
# test.df <- test.df %>%
#   mutate(
#     Xm1 = lag(tX,n=1L,default=0),
#     Zm1 = lag(Z,n=1L,default=0)
#     )
# coxph(Surv(time, Y) ~ tX + Z + Xm1 + Zm1, data=test.df)
# # data structure: use Vset==0 or 1
# summary(test.df$ID)
# with(subset(test.df, int==1), table(Vset)/1000)
# head(test.df %>% 
#        filter(Vset==1) %>% 
#        select(ID,int,Vset,gold,tX,mX,gX),
#      12)
# # check summary statistics
# for (i in 0:1) {
#   print(
#     test.df %>% 
#       filter(Vset==i) %>%
#       select(Vset,int,tX,mX,gold,gX) %>%
#       mutate(
#         sens = sum(tX==1 & mX==1) / sum(tX==1),
#         spec = sum(tX==0 & mX==0) / sum(tX==0)
#       ) %>%
#       summarise_all(~mean(.,na.rm=T))
#   )
# }
