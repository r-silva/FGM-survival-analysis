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

listSec    <- list()
listSe     <- list()

listupper <- list(list())
listlower <- list(list())

# Standard errors (manual) ------------------------------------------------

# for (c in 1:length(listdta)){ 

c <- 3

Clusters <- length(unique(dfList[[c]]$v021))   # Nr of clusters
Strata   <- length(unique(dfList[[c]]$v022))   # Nr of strata

# weight df with wgt (DHS sample weight/1000000) and use df from dfList for country i

df_weight <- dfList[[c]] %>%
  filter(!is.na(time))

# try to recode
df_weight <- df_weight %>%
  mutate(try = ifelse(g121 == 1, 0, 1)) %>%
  select(c("v021", "v022", "time", "wgt", "try"))

df_weight[, 5] <- df_weight[, 5]*df_weight[, 4]

colnames(df_weight)  <- c("n", "H", "time", "M", "g121")



# Generate vector with number of cluster in every stratum

n <- df_weight %>%
  group_by(H)  %>%
  mutate(n.h = n_distinct(n)) %>%
  select(c("H","n.h"))

n <- as.data.frame(n[-which(duplicated(n)), ])

n <- n[order(n$H), ]


for (j in 1:15){ # 0:14 years, fixed for every country

j <- 2
  # create empty dataframe for se calculation
  X <- nRisk[[3]][j]    # j year, takes into account censoring
  R <- nSurv[[3]][j]    # prob of surviving until age j-1
  
  # create data frame with nr of cases and 'positive answers'
  data <- df_weight %>%
    mutate(y = ifelse(time == (j - 1) & g121 != 0, g121, 0)) %>%     ##### Number of positive answers to surviving
    
    mutate(x = ifelse(  time == (j - 2)  & time == (j - 3)  &
                        time == (j - 4)  & time == (j - 5)  & time == (j - 6)  &
                        time == (j - 7)  & time == (j - 8)  & time == (j - 9)  &
                        time == (j - 10) & time == (j - 11) & time == (j - 12) &
                        time == (j - 13) & time == (j - 14), 0, M)) %>%      # At risk population by strata and cluster
    filter(x != 0)

# Generate z.hi matrix from DHS annex
  
  z.hi <- data %>%
    group_by(n) %>%
    mutate(y.hi = sum(y, na.rm=TRUE)) %>%
    mutate(x.hi = sum(x, na.rm=TRUE)) %>%
    select(c("H","n","y.hi","x.hi"))
  
  z.hi <- z.hi[-which(duplicated(z.hi)), ]

 # z.hi  <- z.hi%>%
   # dplyr::mutate(new = y.hi - R*x.hi)%>%
    #select(c("H","n","new"))%>%
   # spread(n,new)
  #z.hi[is.na(z.hi)] <- 0
  
  y.hi <- z.hi %>%
  select(c("H","n","y.hi"))  %>%
  spread(n,y.hi)
  y.hi[is.na(y.hi)] <- 0
  
 x.hi <- z.hi %>%
     select(c("H","n","x.hi"))  %>%
     spread(n,x.hi)
  x.hi[is.na(x.hi)] <- 0
  

  ################
  
  # Generate z.h matrix from DHS annex
  
  z.h <- data %>%
   group_by(H) %>%
   mutate(y.h = sum(y, na.rm=TRUE)) %>%
   mutate(x.h = sum(x)) %>%
   mutate(n.h = n_distinct(n))%>%
   select(c("H","y.h","x.h","n.h"))#%>%
   #mutate(new = (y.h - R*x.h)/n.h)%>%
   #select(c("H","new"))
   #z.h <- as.data.frame(z.h[-which(duplicated(z.h)), ])
   #z.h <- z.h[order(z.h$H),]
   
  y.h <- z.h %>%
   select(c("H","y.h"))%>%
   group_by(H)
  y.h <- as.data.frame(y.h[-which(duplicated(y.h)), ])
  y.h <- y.h[order(y.h$H),]
  
  x.h <- z.h %>%
   select(c("H","x.h"))%>%
   group_by(H)
  x.h <- as.data.frame(x.h[-which(duplicated(x.h)), ])
  x.h <- x.h[order(x.h$H),]
  
  

  
  
# Calculation se using Taylor linearization, 
  
  temp.a = matrix(NA, Strata, Clusters) # Matrix of NAs with H number of rows and n number of columns
  temp.b = c() # Empty vector
  
# as in DHS annex
  
  for (h in 1:Strata){ # loop through strata
  
    
  for (cl in 1:Clusters) { # loop through cluster
  
   temp.a[h, cl] <- (y.hi[h, cl + 1] - R*x.hi[h, cl + 1])^2 
  
  }
  
   temp.b[h] <- (n[h,2]/(n[h,2]-1))*(rowSums(temp.a, na.rm = T)[h]-(((y.h[h,2] - R*x.h[h,2]))^2)/n[h,2])
  
  }



# Pop div
  
  #for (h in 1:Strata){ # loop through strata
    
   # for (cl in 1:Clusters) { # loop through cluster
      
    # temp.a[h, cl] <- (z.hi[h, cl+1]-z.h[h,2])^2 
      
  #  }
    
   # temp.b[h] <- (n[h,2]/(n[h,2]-1))*(rowSums(temp.a, na.rm = T)[h])
    
  #}
  

  
  
  #####################

  
#  for (h in 1:Strata){ # loop through strata
    
    #for (cl in 1:Clusters) { # loop through cluster
      
     # temp.a[h, cl] <- (z.hi[h, cl])^2 
      
    #}
    
   # temp.b[h] <- (n[h,2]/(n[h,2]-1))*(rowSums(temp.a, na.rm = T)[h]-(z.h[h,2])^2/n[h,2])
    
  #}
  
  var.p <- (1/X^2)*sum(temp.b) # remove x^2
  se <- sqrt(var.p) 
  
se

# Store Ses and confidence intervals in lists

listSec[[j]]     <- se # when full loop will run, remove 1 and put c
}

listSe[[c]]   <- listSec
listupper[[c]][j]  <- R + 1.96*se # when full loop will run, remove 1 and put c
listlower[[c]][j]  <- R - 1.96*se # when full loop will run, remove 1 and put c


#} # loop over ages ends

#} # loop over countries ends