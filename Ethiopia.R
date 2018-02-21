library(plyr)
library(dplyr)
library(magrittr)
library(tidyverse)
library(ggrepel)
library(treemap)
library(gdata)
library(knitr)
library(devtools)
library(xtable)
library(gridExtra)
library(gridBase)
library(readxl)
library(gtable)
library(grid)
library(Rmisc)
library(survey)
library(foreign)
library(memisc) #https://cran.r-project.org/web/packages/memisc/vignettes/items.html
library(reshape2)
library (prettyR)
library(survminer)
#Read in data

eth2016_wide <- read.dta("C:/Users/weny/Google Drive/2018/FGM/01 -Survival Analysis/01 -Ethiopia_DHS 2016/STATA/ETIR70FL.DTA", convert.factors=FALSE)

#eth2016 <- read.dta("C:/Users/Kathrin Weny/Google Drive/FGM/Survival Analysis/Ethiopia_DHS 2016/STATA/ETIR70FL.dta", convert.factors=FALSE)

options(scipen = 999) #disable scientif notation in R

#select variables and create dataset for birth reshape and FGM reshape

eth2016_wide_allchildren <- eth2016_wide%>%
  dplyr::select(c(1,3:4,7,25:29,51,4864,4866,4867:4870,4880,66:67,77,87,102,3624,227:246,307:326))

eth2016_wide_fgm <- eth2016_wide %>%
  dplyr::select(c(1,4894:4914,4915:4994))

#Reshape FGM module

eth2016_long_fgm<-reshape(eth2016_wide_fgm, 
                          varying=c(3:102), direction="long", idvar="caseid", sep="_", timevar="order")
#Delete order
eth2016_long_fgm<-eth2016_long_fgm%>%
  dplyr::select(-order)

#Rename GIDX to order
colnames(eth2016_long_fgm)[3] <- "order"

#Create ID for merge
eth2016_long_fgm<-eth2016_long_fgm%>%
  mutate(id_match = paste(caseid,order,sep=""))

#Reshape ALL CHILDREN

eth2016_long_allchildren<-reshape(eth2016_wide_allchildren, 
                                  varying=c(24:63), direction="long", idvar="caseid", sep="_", timevar="order")

#Create ID for merge
eth2016_long_allchildren<-eth2016_long_allchildren%>%
  mutate(id_match = paste(caseid,order,sep=""))

#Merge data sets

df <- merge(eth2016_long_fgm,eth2016_long_allchildren,by="id_match")
colnames(df)[2] <- "caseid"

#filter for age bellow 14
#df <- df%>%
# filter(b8<15)

#add time to event column

time <- ifelse(df$g121 == 0, df$b8,
               ifelse(df$g121==1,df$g122,"na"))

time <- as.numeric(time)
df <- df%>%
  mutate(time = time)

#Recode NA
df$time <- ifelse(df$time==98,NA,df$time)

#Add id column
df <- df%>%
  mutate(id = 1)

#create weight variable from v005
df$wgt <- as.numeric(df$v005/1000000)

#Bring in complex survey format

library(pryr)

memory.limit(550000) 

eth2016_dhs_design <- svydesign(id=~v021,
                                strata=df$v022, #V022 Sample strata for standard errors; same as V023 stratification used for sampling design
                                weights = df$wgt, #weight expressed in 6 decimals
                                data=df)

#create a small sample survey to calculate standard errors
small_survey <- svydesign(id=~v021, strata = df$v022, 
                          variables = ~g121 + time,
                          weight = df$wgt, data=df)

#object.size(small_survey)
s_small <- svykm(Surv(time,g121>0)~1,design=small_survey,se=TRUE)

confint(s_small,parm=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),level=0.95)
confint(s_small,parm=c(1.4),level=0.95)
plot(s_small, pars=NULL, ci=TRUE)

