
# Copyright statement comment ---------------------------------------------

# @UNFPA


# Author comment ----------------------------------------------------------

# work is progress-
# Standard error calculation takes a long time (more than 1 hour).
# For countries other than Tanzania, Cote d'Ivoire and Togo, complex sample SE have been calculated based on 
# a sample of the actual data available. In addition, Se based on (weighted) random sampling have been added.
# This proivdes an upper an a lower boundary for SEs


# File description, purpose of code, inputs and output --------------------

# this code has been written in order to estimate the survival curves for girls age 0-14 in Ethiopia
# including non-sampling errors (confidence intervals)
# the calculations are based on the microdataset of the 2016 DHS in Ethiopia


# Source codes, libraries, global options and working directory -----------

library(plyr)
library(dplyr)
library(magrittr)
library(tidyverse)
library(survey)
library(foreign)
library(reshape2)
library(survminer)
library(data.table)

options(scipen = 999)  # disable scientific notation in R

#if a primary sapling unit has only a single observation, R will crash
# option adjust calculates conservative standard errors
# reference: http://faculty.washington.edu/tlumley/survey/example-lonely.html
options(survey.lonely.psu = "adjust")

setwd("C:/Users/weny/Google Drive/2018/FGM/01 -Survival Analysis/03 -Data/MICS/fgm")

listdta_fg <- dir(pattern = "*.sav") 

# Function definitions ----------------------------------------------------

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
  # uses a lot of memomry, in the case of the 11 DHS files 3 661 612 264 bytes, so better do not store in RAM
  # 
  # Args:
  # x = directory
  # Returns:
  #  ldf     = list of all data frames in the directory named by x
  
  ldf_fg <- list() # creates a list
  
  for (k in 1 : length(listdta_fg)){
    ldf_fg[[k]] <- read.spss(listdta_fg[k], to.data.frame = TRUE)
  }
  
  names(ldf_fg) <- x
  
  return(ldf_fg)
  
}


ReadSingleDTA <- function(x){
  # reads in one DTA file from a given directory
  # 
  # Args:
  #  x = list index of file to be read into R
  #
  # Returns:
  #  single dataset
  
  data <- read.spss(listdta_fg[x], to.data.frame = TRUE)
  
  return(data)
  
}


randomRows <- function(x,y){
  # randomly selects y records from original filex
  # 
  # Args:
  #  x = dataframe
  #  y = number of lines to be selected
  #
  # Returns:
  #  dataframe with y rows
  #  reports error when more cases than available rows are selected
  
  if(y > nrow(x)) stop("too many cases select")
  
  return(x[sample(nrow(x),y),])
  
}


# Lists -------------------------------------------------------------------

ConInList          <- list()

# Executed statements -----------------------------------------------------

# Benin -------------------------------------------------------------------

df <- ReadSingleDTA(1) %>%
  mutate(id = 1)

# Recode variables of interest 

df$FG13 <- as.numeric(as.character(df$FG13))#Age of daughter

df$FG16 <- as.numeric(as.character(df$FG16))#Age of daughter at circumcisison

df[,'FG15'] <- as.character(df[,'FG15'])#Daughter circumcised or not
df <- df %>%
  mutate(fgm_status = ifelse(as.character(df$FG15)== "Non", 0,
                             ifelse(as.character(df$FG15) =="Oui",1,NA)))


time <- ifelse(df$fgm_status == 0, df$FG13,
               ifelse(df$fgm_status==1,df$FG16,NA))
time <- as.integer(time)

unique(df$fgm_status)
unique(df$time)

df <- df%>%
  mutate(time=time)

# Bring in complex survey format

dhs_design <- svydesign(id = ~HH1 + ~HH2, #primary and secondary sampling units, reference: https://rpubs.com/trjohns/survey-cluster
                        weights = ~wmweight, #weight expressed in 6 decimals
                        data = df)

# Without standard errors (we can use full data set)
s1 <- svykm(Surv(time, fgm_status>0) ~1, design=dhs_design, se = T)

# Plot

a <- c(14,15)
b <- c(0.9854752 , 0.9854752)
c <- c(0.99248095, 0.99248095)
d <- c(0.9994867 , 0.9994867)

ab <- data.frame(a,b)
ac <- data.frame(a,c)
ad <- data.frame(a,d)

plot(s1, main = "KM estimates for Benin",
     xlab = "Study time", ylab = "Probability of not experiencing FGM", xlim = c(0,15),
     ylim = c(0.96,1)) 
axis(1, at = 1:15, labels = 1:15)
lines(ab, lty = 2)
lines(ac)
lines(ad, lty = 2)
legend(0, 0.2, legend=c("KM estimate", "95% confidence interval"),
       col=c("black", "black"), lty = 1:2)
title(sub = "Data: MICS 2014", adj = 1, line = 4, font = 3)

# Standard errors

ses <- confint(s1, parm = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,15), level = 0.95)
ConInList[[1]] <- ses  
ses


# Central African Republic ------------------------------------------------

df <- ReadSingleDTA(2) %>%
  mutate(id = 1)

# Recode variables of interest 

df$FG13 <- as.numeric(as.character(df$FG13))#Age of daughter

df$FG16 <- as.numeric(as.character(df$FG16))#Age of daughter at circumcisison

df[,'FG15'] <- as.character(df[,'FG15'])#Daughter circumcised or not
df <- df %>%
  mutate(fgm_status = ifelse(as.character(df$FG15)== "Non", 0,
                             ifelse(as.character(df$FG15) =="Oui",1,NA)))

# Create time to event or time to censoring variable

time <- ifelse(df$fgm_status == 0, df$FG13,
               ifelse(df$fgm_status == 1, df$FG16, NA))

time <- as.integer(time)

df <- df%>%
  mutate(time=time)

# Bring in complex survey format

dhs_design <- svydesign(id = ~HH1 + ~HH2, #primary and secondary sampling units, reference: https://rpubs.com/trjohns/survey-cluster
                        weights = ~WMWEIGHT, #weight expressed in 6 decimals
                        data = df)

# KM estimates
s2 <- svykm(Surv(time, fgm_status>0) ~1, design=dhs_design, se = T)

# Plot

a <- c(0, 1, 2, 3)
b <- c(1, 1, 1, 1)

ab <- data.frame(a,b)

c <- c(14,15)
d <- c(0.9055838 , 0.9055838)
e <- c(0.92524715, 0.92524715)
f <- c(0.9449105 , 0.9449105)

cd <- data.frame(c,d)
ce <- data.frame(c,e)
cf <- data.frame(c,f)


plot(s2, main = "KM estimates for Central African Republic",
     xlab = "Study time", ylab = "Probability of not experiencing FGM", xlim = c(0,15))
axis(1, at = 1:15, labels = 1:15)
lines(ab)
lines(cd, lty=2)
lines(ce)
lines(cf, lty=2)
legend(0, 0.2, legend=c("KM estimate", "95% confidence interval"),
       col=c("black", "black"), lty = 1:2)
title(sub = "Data: MICS 2010", adj = 1, line = 4, font = 3)

# Standard errors

ses <- confint(s2, parm = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,15), level = 0.95)
ConInList[[2]] <- ses  
ses


# Ghana -------------------------------------------------------------------

df <- ReadSingleDTA(3) %>%
  mutate(id = 1)

# Recode variables of interest 

df$FG13 <- as.numeric(as.character(df$FG13))#Age of daughter

df$FG16 <- as.numeric(as.character(df$FG16))#Age of daughter at circumcisison

unique(df$FG15)
df[,'FG15'] <- as.character(df[,'FG15'])#Daughter circumcised or not
df <- df %>%
  mutate(fgm_status = ifelse(as.character(df$FG15)== "No", 0,
                             ifelse(as.character(df$FG15) =="Yes",1,NA)))

# Create time to event or time to censoring variable

time <- ifelse(df$fgm_status == 0, df$FG13,
               ifelse(df$fgm_status==1,df$FG16,"na"))
time <- as.integer(time)

df <- df%>%
  mutate(time = time)

#Bring in complex survey format

dhs_design <- svydesign(id = ~HH1 + HH2, #primary and secondary sampling units,
                        weights = df$wmweight, #weight expressed in 6 decimals
                        data = df)

# KM estimates
s3 <- svykm(Surv(time, fgm_status>0) ~1, design=dhs_design, se = T)

# plot 
a <- c(12, 13, 14, 15)
b <- c(0.9904473 , 0.9904473 , 0.9904473, 0.9904473)
c <- c(0.99309185, 0.99309185, 0.99309185, 0.99309185)
d <- c(0.9957364 , 0.9957364 , 0.9957364, 0.9957364)

ab <- data.frame(a,b)
ac <- data.frame(a,c)
ad <- data.frame(a,d)

plot(s3, main = "KM estimates for Ghana",
     xlab = "Study time", ylab = "Probability of not experiencing FGM", xlim = c(0,15))
axis(1, at = 1:15, labels = 1:15)
lines(ab, lty = 2)
lines(ac)
lines(ad, lty = 2)
legend(0, 0.2, legend=c("KM estimate", "95% confidence interval"),
       col=c("black", "black"), lty = 1:2)
title(sub = "Data: MICS 2011", adj = 1, line = 4, font = 3)

# Zoom

plot(s3, main = "KM estimates for Ghana",
     xlab = "Study time", ylab = "Probability of not experiencing FGM", xlim = c(0,15),
     ylim = c(0.96,1))
axis(1, at = 1:15, labels = 1:15)
lines(ab, lty = 2)
lines(ac)
lines(ad, lty = 2)
legend(0, 0.2, legend=c("KM estimate", "95% confidence interval"),
       col=c("black", "black"), lty = 1:2)
title(sub = "Data: MICS 2011", adj = 1, line = 4, font = 3)

# Standard errors

ses <- confint(s3, parm = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,15), level = 0.95)
ConInList[[3]] <- ses  
ses


