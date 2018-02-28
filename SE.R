# Copyright statement comment ---------------------------------------------

# @UNFPA


# Author comment ----------------------------------------------------------

# Due to continuing issues with RAM when calculating SE with the Survey Package
# This is a fiel which hardcodes Taylor linearization SE


# File description, purpose of code, inputs and output --------------------

# Uses output from DHS.R and calculates SE at each single year in survival curve


# Source codes, libraries, global options and working directory -----------

setwd("C:/Users/weny/Google Drive/2018/FGM/01 -Survival Analysis/FGM survival analysis")

source("DHS.R") 

# Standard errors (manual) ------------------------------------------------

#for (i in 1:length(listdta)){ This will be used for a loop through countries later

i <- 1  
# create empty dataframe for se calculation
C <- length(unique(df$v021)) # Nr of clusters
H <- length(unique(df$v022)) # Nr of strata
X <- nRisk[[i]][1]   # Sample size, for now 1 == 0(year), will be used for a loop later
R <- nEvent[[i]]/X # Estimate

# weight df with wgt (DHS sample weight/1000000) and use df from dfList for country i
df_weight <- dfList[[i]] %>%
  mutate(case = 1) # case = 1 if anyone sampled
df_weight[,3] <- df[,3]*df[,5]
df_weight[,6] <- df_weight[,6]*df[,5]

# create new dataframe/matrix to store information relevant for se calculation
data <- setNames(data.frame(matrix(ncol = 4, nrow = nrow(df))), c("H","n","y","m"))

# create data frame with nr of cases at time = 0 (we do not have to account for censoring at this point)
data_test <- data %>%
  mutate(H = df$v022) %>%
  mutate(n = df$v021) %>%
  mutate(y = ifelse(df_weight$time == 0, df_weight$g121,0)) %>% # we are nto looping yet, so ==0 
  mutate(m = df_weight$case)# %>% # no censoring yet, in future, we have to take into account

# Generate z.hi matrix from DHS annex
data_test6 <- data_test %>%
  group_by(n) %>%
  mutate(y.hi = sum(y)) %>%
  mutate(m.hi = sum(m)) %>%
  select(c("H","n","y.hi","m.hi"))

data_test6 <- data_test6[-which(duplicated(data_test6)), ]

data_test6  <- data_test6%>%
  mutate(new = y.hi - p*m.hi)%>%
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
  mutate( new = y.h - p*m.h)%>%
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

###

# Function for se calculation usign Taylor linearization

SEest <- function(data){
  
  temp.a = matrix(NA, H, C) # Matrix of NAs with H number of rows and n number of columns
  temp.b = c() # Empty vector
  
  for (h in 1:H){ # loop through strata
    
    for (i in 1:C) { # loop through cluster
      
      temp.a[h, i] <- (z.hi[h, i])^2 
      
    }
    
    temp.b[h] <- (n[h,2]/(n[h,2]-1))*(rowSums(temp.a, na.rm = T)[h]-(z.h[h,2])^2/n[h,2])
    
  }
  
  var.p <- (1/m^2)*sum(temp.b, na.rm = T)
  se <- sqrt(var.p) 
  
  return(c(se = se))
  
}

#} Until we loop through countries, commented out
