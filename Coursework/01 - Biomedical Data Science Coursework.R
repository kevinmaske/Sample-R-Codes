#install.packages("tidyverse")  # if needed
library(tidyverse)

########################
########################
###### PROBLEM 3 #######
########################
########################

## 3.1 ##

# writing the two functions

egfr.mdrd4 <- function(scr, age, sex, ethnic){
  scr <- as.numeric(scr)  # Convert from character format
  n <- length(scr)  # get length of vector input
  
  # Default Case : Male, other ethnicity; a set of vectors
  e.mult <- rep(1, n)
  s.mult <- rep(1, n)
  
  for(i in 1:n){
    # check if female
    if(sex[i] == "Female"){
      s.mult[i] <- 0.742
    }
    
    # check if black ethnicity
    if(ethnic[i] == "Black"){
      e.mult[i] <- 1.212
    }
  }
    
  eGFR <- 175 * scr ^ (-1.154) *
    age ^ (-0.203) * s.mult * e.mult
  return(eGFR)  
}

egfr.ckdepi <- function(scr, age, sex, ethnic){
  scr <- as.numeric(scr)  # Convert from character format
  n <- length(scr)  # get length of vector input

  # Default Case : Male, other ethnicity; a set of vectors
  k <- rep(0.9, n)
  a <- rep(-0.411, n)
  e.mult <- rep(1, n)
  s.mult <- rep(1, n)
  
  # Go through each patient
  for(i in 1:n){
    # check if female
    if(sex[i] == "Female"){
      k[i] <- 0.7
      a[i] <- -0.329
      s.mult[i] <- 1.018
    }
    
    # check if black ethnicity
    if(ethnic[i] == "Black"){
      e.mult[i] <- 1.159
    }
  }
  eGFR <- 141 * pmin(scr/k, 1)^a * pmax(scr/k, 1)^(-1.209) *
    0.993^age * s.mult * e.mult
  
  return(eGFR)  # Return vector of eGFRs
}



## 3.2 ##

scr <- read.csv("scr.csv", header=T)
# note that there are NA entries in all columns of the data; without knowledge
# of why these parts of the data are missing, it is difficult to deal with them
# there is no reason to believe that the value of other variables allows us to 
# determine the values of the missing column, and since none of the variables in the
# given dataset are the actual outcome we want to consider, we will just impute them
# to complete the data set.

# Initialize imputation scr dataframe
scr.imp <- scr

## imputing ethnicity
table(is.na(scr$ethnic))  # 3 NA entries
table(scr$ethnic)

# Other appears far more often than Black,
# so we replace all 3 of these values with "Others"

scr.imp[which(is.na(scr$ethnic)), "ethnic"] <- "Other"

## Imputing Gender
table(is.na(scr.imp$sex))  # 1 NA entry
table(scr$sex)  # Female is the mode, designate the one missing entry to be female

scr.imp[which(is.na(scr$sex)), "sex"] <- "Female"

## Imputing scr
# There is no reason to believe that the patients are ordered in such a way that might
# be correlated with the amount of scr they're given

plot(scr$scr)  # the plot seems mostly flat with some outliers in the first few patients
# we impute any missing entries with the overall mean

table(is.na(scr.imp$scr))  # 18 NA entry
scr.mean <- mean(scr$scr, na.rm = TRUE)

scr.imp[which(is.na(scr$scr)), "scr"] <- scr.mean

plot(scr.imp$scr)  # plot didn't change dramatically, this seems acceptable


# Imputing age

plot(scr$age)  # a totally scattered plot shows there's no relationship between
# patient order and age, so we just use means (age can be a real number anyway)

table(is.na(scr.imp$age))  # 2 NA entry
age.mean <- mean(scr$age, na.rm = TRUE)

scr.imp[which(is.na(scr$age)), "age"] <- age.mean

## Check if there are any remaining missing values ##
table(is.na(scr.imp))

# Data is now complete and ready to process
# Initialize vectors for input
s <- scr.imp$scr
a <- scr.imp$age
x <- scr.imp$sex
e <- scr.imp$ethnic

#### ESTIMATING eGFR with the 2 functions above ####

## MDRD4 Method
mdrd4 <- egfr.mdrd4(s, a, x, e)

md.mean <- mean(mdrd4) %>% round(2)
md.mean  # mean egfr from mdrd4

md.sd <- sd(mdrd4) %>% round(2)
md.sd  # st. dev. from mdrd4

## CKD-EPI method
ckdepi <- egfr.ckdepi(s, a, x, e)

ck.mean <- mean(ckdepi) %>% round(2)
ck.mean  # mean egfr from ckd-epi

ck.sd <- sd(ckdepi) %>% round(2)
ck.sd  # st. dev. from ckdepi

# Correlations
methods.cor <- cor(mdrd4, ckdepi, method = "pearson") %>% round(2)
methods.cor  # Correlation of methods; close to 1 is good since they're consistent

## 3.3 Plots ##

png("eGFR_est.png")

plot(mdrd4, ckdepi, main = "Scatter Plot of eGFR from MDRD4 and CKD-EPI methods",
     xlab = "MDRD4 estimate", ylab = "CKD-EPI estimate", pch = 20)
# Add lines for medians
abline(v=median(mdrd4), h=median(ckdepi), col = "red", lwd = 2)
# quantiles
abline(v=quantile(c(mdrd4), c(0.25, 0.75)), col = c("blue", "violet"), lwd = 2)
abline(h=quantile(c(ckdepi), c(0.25, 0.75)), col = c("blue", "violet"), lwd = 2)
legend(130, 80, c("1st Quantile", "Median", "3rd Quantile"),
       lty = 1, col = c("blue", "red", "violet"))

dev.off()
