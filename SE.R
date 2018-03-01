# Copyright statement comment ---------------------------------------------

# @UNFPA


# Author comment ----------------------------------------------------------

# Due to continuing issues with RAM when calculating SE with the Survey Package, I admit defeat -
# This is a file which hardcodes Taylor linearization SE


# File description, purpose of code, inputs and output --------------------

# Uses output from DHS.R and calculates SE at each single year in survival curve


# Source codes, libraries, global options and working directory -----------

setwd("C:/Users/weny/Google Drive/2018/FGM/01 -Survival Analysis/FGM survival analysis")

source("DHS.R") 

listSe    <- list(list())
listupper <- list(list())
listlower <- list(list())

# Standard errors (manual) ------------------------------------------------

for (c in 1:length(listdta)){ 

  c <- 3
Clusters <- length(unique(dfList[[c]]$v021)) # Nr of clusters
Strata <- length(unique(dfList[[c]]$v022)) # Nr of strata

# weight df with wgt (DHS sample weight/1000000) and use df from dfList for country i
df_weight <- dfList[[c]] %>%
  filter(!is.na(time))

df_weight[, 3] <- df_weight[, 3]*df_weight[, 5]

colnames(df_weight)  <- c("n","H","g121","time","M")

# <- df_weight %>%
 # mutate(y = ifelse(time == 1 & g121 != 0, g121, 0))%>%
  #mutate(m = ifelse(time == 0 , 0, M))

#colSums(test)


for (j in 1:15){ # 0:14 years, fixed for every country

  # create empty dataframe for se calculation
  X <- nRisk[[3]][j]    # j year, takes into account censoring
  R <- nEvent[[3]][j]/X # Estimate of proportion/risk
  
  # create data frame with nr of cases
  data <- df_weight %>%
    mutate(y = ifelse(time == (j - 1) & g121 != 0, g121, 0)) %>%     # Number of positive answers
    mutate(m = ifelse(time == (j - 2), 0, M))      # At risk population by strata and cluster
   # filter(m != 0)%>%
    #filter(!is.na(m) | !is.na(y))

  # Generate z.hi matrix from DHS annex
  z.hi <- data %>%
    group_by(n) %>%
    mutate(y.hi = sum(y, na.rm=TRUE)) %>%
    mutate(m.hi = sum(m, na.rm=TRUE)) %>%
    select(c("H","n","y.hi","m.hi"))
  
  z.hi <- z.hi[-which(duplicated(z.hi)), ]
  
  z.hi  <- z.hi%>%
    dplyr::mutate(new = y.hi - R*m.hi)%>%
    select(c("H","n","new"))%>%
    spread(n,new)%>%
    select(-H)
  
  z.hi <- z.hi
  z.hi[is.na(z.hi)] <- 0
  
  ################
  
  # Generate z.h matrix from DHS annex
  
  data_test7 <- data %>%
    group_by(H) %>%
    mutate(y.h = sum(y, na.rm=TRUE)) %>%
    mutate(m.h = sum(m)) %>%
    select(c("H","y.h","m.h"))%>%
    mutate( new = y.h - R*m.h)%>%
    select(c("H","new"))
  
  data_test7 <- as.data.frame(data_test7[-which(duplicated(data_test7)), ])
  
  z.h <- data_test7[order(data_test7$H),]
  
  #####
  
  # Generate vector with number of cluster in every stratum
  
  n <- data%>%
    group_by(H)%>%
    mutate(n.h = n_distinct(n))%>%
    select(c("H","n.h"))
  
  n <- as.data.frame(n[-which(duplicated(n)), ])
  
  n <- n[order(n$H),]
  
  
# Calculation se using Taylor linearization

  temp.a = matrix(NA, Strata, Clusters) # Matrix of NAs with H number of rows and n number of columns
  temp.b = c() # Empty vector
  
  for (h in 1:Strata){ # loop through strata
    
    for (cl in 1:Clusters) { # loop through cluster
      
      temp.a[h, cl] <- (z.hi[h, cl])^2 
      
    }
    
    temp.b[h] <- (n[h,2]/(n[h,2]-1))*(rowSums(temp.a, na.rm = T)[h]-(z.h[h,2])^2/n[h,2])
    
  }
  
  var.p <- (1/X^2)*sum(temp.b, na.rm = T)
  se <- sqrt(var.p) 
  


# Store Ses and confidence intervals in lists

listSe[j]     <- se # when full loop will run, remove 1 and put c
}
listupper[[c]][j]  <- R + 1.96*se # when full loop will run, remove 1 and put c
listlower[[c]][j]  <- R - 1.96*se # when full loop will run, remove 1 and put c


} # loop over ages ends

} # loop over countries ends