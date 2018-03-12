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


output <- list()

#data <-  setNames(data.table(matrix(nrow = 14, ncol=3 )),
 #                 c("me", "se", "confint"))

data <- data.frame(matrix(NA, nrow = 14, ncol = 4))
                   
names(data) <- c("me", "se", "confint_low", "confint_high")

#for (c in 1:length(listdta)){ # start loop through countries
 c <- 11

   df <- dfList[[c]] %>%
    filter(!is.na(time))%>%
    select(c("v021","v022","g121","time","wgt"))
  
   for (j in 1:14){ # Loop through ages
     
     # subset data frame: Throw out everybody who was censored before year j
     # Create new variable of women who have survived at age j
    df <- df %>%
      filter(time >= j | (time < j & g121 == 1)) %>% # remove censoring
     mutate(interest = ifelse((time > j) | (time == j & g121 == 0 ), 1, 0)) # survivors until time j
    # mutate(interest = ifelse(time<=j & g121 == 1, 1, 0))
    survey <- svydesign(id             = ~v021, 
                         strata         = df$v022, 
                         weight         = df$wgt,
                         data           = df)
    
    a <- svymean(~interest, survey, deff=TRUE, na.rm=FALSE)

    
    data[j,1] <- as.numeric(a[1])
    data[j,2] <- as.numeric(sqrt(attributes(a)$var))
    data[j,3] <- confint(a)[[1]]
    data[j,4] <- confint(a)[[2]]
    
   }
  data 1- 0.9918613
output[[c]] <- data   
   
   
#}


# test 

SmallSurvey <- svydesign(id             = ~v021, 
                        strata          = df$v022, 
                        weight          = df$wgt,
                         data           = df)

SmallSurvival <- svykm(Surv(time, g121 > 0) ~ 1, design = SmallSurvey, se = FALSE)
plot(SmallSurvival)
str(SmallSurvival)



