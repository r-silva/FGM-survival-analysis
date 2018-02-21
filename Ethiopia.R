###########
# copyright statement comment
###########

# @UNFPA

###########
# author comment
###########

# work is unfinished
# standard error calculation takes a long time

###########
# file description comment, including purpose of program, inputs, and outputs
###########

# this code has been written in order to estimate the survival curves for girls age 0-14 in Ethiopia
# including non-sampling errors (confidence intervals)
# the calculations are based on the microdataset of the 2016 DHS in Ethiopia


###########
# source() and library() statements
###########

library(plyr)
library(dplyr)
library(magrittr)
library(tidyverse)
library(gridExtra)
library(gridBase)
library(readxl)
library(survey)
library(foreign)
library(reshape2)
library(survminer)

###########
# function definitions
###########

#time function
#create id

###########
# executed statements
###########


# read in data
eth2016_wide <- read.dta("C:/Users/weny/Google Drive/2018/FGM/01 -Survival Analysis/01 -Ethiopia_DHS 2016/STATA/ETIR70FL.DTA", convert.factors=FALSE)
# eth2016 <- read.dta("C:/Users/Kathrin Weny/Google Drive/FGM/Survival Analysis/Ethiopia_DHS 2016/STATA/ETIR70FL.dta", convert.factors=FALSE)

options(scipen = 999)  # disable scientific notation in R

# select variables and create dataset for birth reshape and FGM reshape

eth2016_wide_allchildren <- eth2016_wide %>%
  dplyr::select(c(1, 3:4, 7, 25:29, 51, 4864, 4866, 4867:4870, 4880, 66:67, 77, 87, 102, 3624, 227:246, 307:326))

eth2016_wide_fgm <- eth2016_wide %>%
  dplyr::select(c(1, 4894:4914, 4915:4994))

# reshape FGM module

eth2016_long_fgm <- reshape(eth2016_wide_fgm, varying = c(3:102), direction = "long", idvar = "caseid", sep = "_", timevar = "order")

# delete order
eth2016_long_fgm <- eth2016_long_fgm %>%
  dplyr::select(-order)

# rename GIDX to order
colnames(eth2016_long_fgm)[3] <- "order"

# create ID for merge
eth2016_long_fgm <- eth2016_long_fgm %>%
  mutate(id.match = paste(caseid, order, sep = ""))

# reshape ALL CHILDREN

eth2016_long_allchildren <- reshape(eth2016_wide_allchildren, varying = c(24:63), direction = "long", idvar = "caseid", sep = "_", timevar="order")

# create ID for merge
eth2016_long_allchildren <- eth2016_long_allchildren %>%
  mutate(id.match = paste(caseid, order, sep = ""))

# merge data sets

df <- merge(eth2016_long_fgm, eth2016_long_allchildren, by = "id.match")

colnames(df)[2] <- "caseid"

# add time to event column

time <- ifelse(df$g121 == 0, df$b8, ifelse(df$g121 == 1, df$g122, "na"))

time <- as.numeric(time)

df <- df %>%
  mutate(time = time)

# add id column
df <- df %>%
  mutate(id = 1)

# create weight variable from v005
df$wgt <- as.numeric(df$v005 / 1000000)

memory.limit(550000) 

# bring in complex survey format
eth2016_dhs_design <- svydesign(id      = ~v021,
                                strata  = df$v022,  # V022 Sample strata for standard errors
                                weights = df$wgt,  # weight expressed in 6 decimals
                                data    = df)

# create a small sample survey to calculate standard errors
SmallSurvey <- svydesign(id             = ~v021, 
                         strata         = df$v022, 
                         variables      = ~g121 + time,
                         weight         = df$wgt,
                         data           = df)

SmallSurvival <- svykm(Surv(time, g121 > 0) ~ 1, design = SmallSurvey, se = TRUE)

confint(SmallSurvival, parm = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), level = 0.95)
confint(SmallSurvival, parm = c(1.4), level = 0.95)
plot(SmallSurvival, pars = NULL, ci = TRUE)

