###########
# copyright statement comment
###########

# @UNFPA

###########
# author comment
###########

# work is unfinished
# standard error calculation takes a long time (more than 1 hour)

###########
# file description comment, including purpose of program, inputs, and outputs
###########

# this code has been written in order to estimate the survival curves for girls age 0-14 in Ethiopia
# including non-sampling errors (confidence intervals)
# the calculations are based on the microdataset of the 2016 DHS in Ethiopia

###########
# source(), library() and options statements
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

options(scipen = 999)  # disable scientific notation in R

###########
# function definitions
###########

CalculateTimetoEvent <- function(x,y,z){
  # calculates the time (years) that has passed until a certain event(FGM) occurs
  # and in the case, there was no event, uses proxy (age)
  #
  #Args:
  #   x = event indicator variable
  #   y = age at wich event occured
  #   z = substitute in case no event
  #
  #Returns:
  # new dataframe with time to event column
  
  time <- ifelse(x == 0, z, ifelse(x == 1, y, NA))
  
  time <- as.numeric(time)
  
  return(time)

}

ReadListofDTA <- function(x){
  #Rreads in all DTA files in the directors and stores them in a named list
  #
  #Args:
  # x = list of character values representing names of all DTA files in the directory
  #Returns:
  #  ldf     = list of all data frames in the directory named by x
  
  ldf <- list() # creates a list
  
  for (k in 1:length(x)){
    ldf[[k]] <- read.dta(x[k], convert.factors=FALSE)
  }
  
  names(ldf) <- x
  
  return(ldf)

}

###########
# executed statements
###########

# set working directory
setwd("C:/Users/weny/Google Drive/2018/FGM/01 -Survival Analysis/03 -Data/DHS")

# read in stata datafiles
listdta <- dir(pattern = "*.DTA") # creates the list of all the csv files in the directory

ldf <- ReadListofDTA(listdta) 

# select relevant data file
  # ethiopia
#wide <- ldf[["ethiopia.DTA"]]

# try a loop

# create empty list to store results
SmallSurvivalList <- list()
PlotList <- list()

for (i in length(ldf)){

wide <- ldf[[i]]

# select variables and create dataset for birth reshape and FGM reshape
wide_allchildren <- wide %>%
  dplyr::select(c(1, 3:4, 7, 25:29, 51, 4864, 4866, 4867:4870, 4880, 66:67, 77, 87, 102, 3624, 227:246, 307:326))

wide_fgm <- wide %>%
  dplyr::select(c(1, 4894:4914, 4915:4994))

# reshape FGM module
long_fgm <- reshape(wide_fgm, varying = c(3:102), direction = "long", idvar = "caseid", sep = "_", timevar = "order")

# delete order
long_fgm <- long_fgm %>%
  dplyr::select(-order)

# rename GIDX to order
colnames(long_fgm)[3] <- "order"

# create ID for merge
long_fgm <- long_fgm %>%
  mutate(id.match = paste(caseid, order, sep = ""))

# reshape ALL CHILDREN

long_allchildren <- reshape(wide_allchildren, varying = c(24:63), direction = "long", idvar = "caseid", sep = "_", timevar="order")

# create ID for merge
long_allchildren <- long_allchildren %>%
  mutate(id.match = paste(caseid, order, sep = ""))

# merge data sets

df <- merge(long_fgm, long_allchildren, by = "id.match")

colnames(df)[2] <- "caseid"

# add time to event column

time <- CalculateTimetoEvent(df$g121,df$g122,df$b8)

df <- df %>%
  mutate(time = time)

# add id column
df <- df %>%
  mutate(id = 1)

# create weight variable from v005
df$wgt <- as.numeric(df$v005 / 1000000)

# increase memory
memory.limit(550000) 

# create a small sample survey to calculate standard errors
SmallSurvey <- svydesign(id             = ~v021, 
                         strata         = df$v022, 
                         variables      = ~g121 + time,
                         weight         = df$wgt,
                         data           = df)

SmallSurvival <- svykm(Surv(time, g121 > 0) ~ 1, design = SmallSurvey, se = TRUE)
SmallSurvivalList[[i]]     <- SmallSurvival

confint(SmallSurvival, parm = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), level = 0.95)

plot <- plot(SmallSurvival, pars = NULL, ci = TRUE)
PlotList[[i]]     <- plot
}
