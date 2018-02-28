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
  mutate(case = 1) # case = 1 if anyone sampled
df_weight[, 3] <- df_weight[, 3]*dfList[[c]][, 5]
df_weight[, 6] <- df_weight[, 6]*dfList[[c]][, 5]


# create new dataframe/matrix to store information relevant for se calculation
data <- setNames(data.frame(matrix(ncol = 4, nrow = nrow(dfList[[c]]))), c("H","n","y","m"))

###

### ONLY SE FOR 0 correct, CENSORING AND EVENTS OUT OF REST TOO

for (j in 1:15){ # 0:14 years, fixed for every country
  
  j <- 15
  # create empty dataframe for se calculation
  X <- nRisk[[3]][j]   # j year, takes into account censoring
  R <- nEvent[[3]][j]/X # Estimate
  
  # create data frame with nr of cases at time = 0 (we do not have to account for censoring at this point)
  data_test <- data %>%
    mutate(H = dfList[[c]]$v022) %>%
    mutate(n = dfList[[c]]$v021) %>%
    mutate(y = ifelse(df_weight$time == j, df_weight$g121,0)) %>% # we are nto looping yet, so ==0 
    mutate(m = df_weight$case)# %>% # TAKE INTO ACCOUNT CENSORING
  
  # Generate z.hi matrix from DHS annex
  data_test6 <- data_test %>%
    group_by(n) %>%
    mutate(y.hi = sum(y, na.rm=TRUE)) %>%
    mutate(m.hi = sum(m, na.rm=TRUE)) %>%
    select(c("H","n","y.hi","m.hi"))
  
  data_test6 <- data_test6[-which(duplicated(data_test6)), ]
  
  data_test6  <- data_test6%>%
    dplyr::mutate(new = y.hi - R*m.hi)%>%
    select(c("H","n","new"))%>%
    spread(n,new)%>%
    select(-H)
  
  z.hi <- data_test6
  z.hi[is.na(z.hi)] <- 0
  
  ################
  
  # Generate z.h matrix from DHS annex
  
  data_test7 <- data_test %>%
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
  
  n <- data_test%>%
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