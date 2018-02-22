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
  # reads in all DTA files in the directory and stores them in a named list
  # uses a lot of memomry, in the case of the 11 DHS files 3 661 612 264 bytes
  # 
  # Args:
  # x = directory
  # Returns:
  #  ldf     = list of all data frames in the directory named by x
  
  listdta <- dir(x,pattern = "*.DTA") # stores names of all DTA files in a list
  
  ldf <- list() # creates a list
  
  for (k in 1 : length(listdta)){
    ldf[[k]] <- read.dta(listdta[k], convert.factors = FALSE)
  }
  
  names(ldf) <- x
  
  return(ldf)
  
}


ReadSingleDTA <- function(x,y){
  # reads in one DTA file from a given directory
  # 
  # Args:
  #  x = A directory
  #  y = list index of file to be read into R
  #
  # Returns:
  #  single dataset
  
  listdta <- dir(x,pattern = "*.DTA") # finds dta files in the directory
  
  data <- read.dta(paste(x,"/",listdta[y], sep=""), convert.factors = FALSE)
  
  return(data)
  
}


###########
# executed statements
###########

# read in stata datafiles

# for pcs with a lot of memory
# ldf <- ReadListofDTA("C:/Users/weny/Google Drive/2018/FGM/01 -Survival Analysis/03 -Data/DHS") 

# create empty list to store results in loop
SmallSurvivalList <- list()
PlotList <- list()

for (i in 1:length(ldf)){

  # for pcs with a lot of memory
  #wide <- ldf[[i]]
  
  # if limited memory
  wide <- ReadSingleDTA("C:/Users/weny/Google Drive/2018/FGM/01 -Survival Analysis/03 -Data/DHS",i) 
  
  # select variables and create dataset for birth reshape and FGM reshape
  wide_allchildren <- wide %>%
    dplyr::select(c("caseid","v001","v002","v005","v021","v022","v023","v024","v025",
                    "bidx_01","bidx_02","bidx_03","bidx_04","bidx_05","bidx_06","bidx_07","bidx_08","bidx_09","bidx_10","bidx_11","bidx_12","bidx_13","bidx_14","bidx_15","bidx_16","bidx_17","bidx_18","bidx_19","bidx_20",
                    "b4_01","b4_02","b4_03","b4_04","b4_05","b4_06","b4_07","b4_08","b4_09","b4_10","b4_11","b4_12","b4_13","b4_14","b4_15","b4_16","b4_17","b4_18","b4_19","b4_20",
                    "b8_01","b8_02","b8_03","b8_04","b8_05","b8_06","b8_07","b8_08","b8_09","b8_10","b8_11","b8_12","b8_13","b8_14","b8_15","b8_16","b8_17","b8_18","b8_19","b8_20"))
  
  wide_fgm <- wide %>%
    dplyr::select(c("caseid","v001",
                    "gidx_01","gidx_02","gidx_03","gidx_04","gidx_05","gidx_06","gidx_07","gidx_08","gidx_09","gidx_10","gidx_11","gidx_12","gidx_13","gidx_14","gidx_15","gidx_16","gidx_17","gidx_18","gidx_19","gidx_20",
                    "g121_01","g121_02","g121_03","g121_04","g121_05","g121_06","g121_07","g121_08","g121_09","g121_10","g121_11","g121_12","g121_13","g121_14","g121_15","g121_16","g121_17","g121_18","g121_19","g121_20",
                    "g122_01","g122_02","g122_03","g122_04","g122_05","g122_06","g122_07","g122_08","g122_09","g122_10","g122_11","g122_12","g122_13","g122_14","g122_15","g122_16","g122_17","g122_18","g122_19","g122_20",
                    "g123_01","g123_02","g123_03","g123_04","g123_05","g123_06","g123_07","g123_08","g123_09","g123_10","g123_11","g123_12","g123_13","g123_14","g123_15","g123_16","g123_17","g123_18","g123_19","g123_20",
                    "g124_01","g124_02","g124_03","g124_04","g124_05","g124_06","g124_07","g124_08","g124_09","g124_10","g124_11","g124_12","g124_13","g124_14","g124_15","g124_16","g124_17","g124_18","g124_19","g124_20"))
  
  
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
  
  long_allchildren <- reshape(wide_allchildren, varying = c(10:69), direction = "long", idvar = "caseid", sep = "_", timevar="order")
  
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
  
  #NAs
  df$time <- ifelse(df$time==98 | df$time==99,NA,df$time)
  
  # create weight variable from v005
  df$wgt <- as.numeric(df$v005 / 1000000)
  
  # increase memory
  memory.limit(10000000) 
  
  # create a small sample survey to calculate standard errors
  SmallSurvey <- svydesign(id             = ~v021, 
                           strata         = df$v022, 
                           variables      = ~g121 + time,
                           weight         = df$wgt,
                           data           = df)
  
  # free memomry for KM estimates
  rm(wide, wide_allchildren,long_allchildren,data, wide_fgm,long_fgm,df,time)
  
  # KM estimate
  SmallSurvival <- svykm(Surv(time, g121 > 0) ~ 1, design = SmallSurvey, se = TRUE)
  SmallSurvivalList[[i]] <- SmallSurvival  # store survival object in list
  
  confint(SmallSurvival, parm = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), level = 0.95)
  
  # plot
  plot <- plot(SmallSurvival, pars = NULL, ci = TRUE)
  PlotList[[i]] <- plot  # store plot in list
}
