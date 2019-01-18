## 
##
## -------------------------------------------
## LOAD PACKAGES
## -------------------------------------------
##


install.packages("lubridate")

library(tidyverse)
library(lubridate)


## ;;
## -------------------------------------------
## READING DATA
## -------------------------------------------
## ;;

load("accident.RData")
load("burlington.RData")
load("list.data.frames.RData")


## --------------------------------------------------
## --------------------------------------------------

## ;;
## ---------------------------------------------
## Q1: -- add your code below 
## ---------------------------------------------
## ;;

## 1.1

# Select only precipitation and date column, then filter accordingly
answer1.1a <- burlington %>%
  select(Date, Precipitation) %>%
  filter(Precipitation>31) 


# Creates a table that has only the location and cause variables
accident.data.loc.and.causes <- accident.data %>%
  select(-Where.and.What...Date..Period..Day, -Latitude, -Longitude, -Timehour)

# Groups the data by location, then gets percentages of each cause in those locations
# Since each cause is indicated by either 1 or 0, the percentage can also be calculated
# by getting the mean of each column
answer1.1b <- accident.data.loc.and.causes %>%
  group_by(Where.and.What...Location.list..Regular) %>%
  summarize_all(mean)

## 1.2

burlington.hour <- burlington %>%
  select(Hour, Precipitation) %>%
  group_by(Hour) %>%
  summarize(average = mean(Precipitation), sample.variance = sd(Precipitation)^2,
            minimum = min(Precipitation), maximum = max(Precipitation))

## 1.3

# Plot the averages obtained in previous item, and adds a bar
# representing the 95% two-tailed normal intervals for each point
burlington.hour %>% ggplot(aes(Hour, average)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = average - 1.96*sqrt(sample.variance),
        ymax = average + 1.96*sqrt(sample.variance)),
    colour = "red",
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FEEDBACK: *sqrt(sample.variance/ no of obs)
#        typo:  colour = "red",
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## 1.4

## Answer: A

## 1.5

# Executes rbind on each entry in the list of df's
df.rbinded <- do.call("rbind", list.dfs)


## ;;
## -------------------------------------------
## Q2: -- add your code below  
## -------------------------------------------
## ;;

## 2.1

# Inversion sampling from Pareto Distribution
rPareto <- function(n, alpha, beta) {
  u <- runif(n)  # Simulate n samples from U(0,1)
  p <- beta / ((1-u)^(1/alpha))  # Inverse of Pareto CDF
  return(p)
}

pareto1.1 <- rPareto(1000, 1, 1)  # Get 1000 samples from Pareto(1,1)
plot(density(pareto1.1))  # Evaluate via graph if dample is from pareto


## 2.2

# Inversion sampling from Rayleigh Distribution
rRayleigh <- function(n, sigma) {
  u <- runif(n)  # Get n samples from U(0,1)
  r <- sigma*sqrt(-2*log(1-u))  # Inverse of Rayleigh CDF
  return(r)
}

rayleigh1 <- rRayleigh(1000, 1) # Get 1000 samples from Rayleigh
plot(density(rayleigh1))

## 2.3

# Inversion sampling from GEV
rgev <- function(n, mu, sigma, xi, tol.xi.limit = 5*10^-2){
  u <- runif(n)  # Get n samples from U(0,1)
  if(abs(xi) < tol.xi.limit) {  # Case where xi -> 0
    g <- mu - sigma*log(-log(u))  # Inverse CDF when xi -> 0
    return(g)
  } else {
    g <- mu + (sigma * (((-log(u))^(-sigma))-1))/sigma  # Inverse CDF otherwise
    return(g)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FEEDBACK: 2nd g should be: X <- mu + (sigma/xi) * ( ( -log(U) )^(-xi) - 1 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gev1 <- rgev(1000, 0, 1, 0)  # 1000 samples from GEV(0,1,0)
gev2 <- rgev(1000, 0, 1, 1)  # 1000 samples from GEV(0,1,1)
gev3 <- rgev(1000, 0, 1, 1/2)  # 1000 samples from GEV(0,1,0.5)

## 2.4

# Define a function that allows us to get the pdf of beta given some x and parameters
pdf.beta <- function(x, alpha, beta) {
  g <- (gamma(alpha+beta)/(gamma(alpha)*gamma(beta)))*x^(alpha-1)*(1-x)^(beta-1)
  return(g)
}

# Rejection Sampling for beta distribution; knowing that rejection sampling
# won't work in the cases when beta has asymptotes, which occur when either
# parameter is less than 1.
rbeta.rs <- function(n, alpha, beta){
  if((alpha < 1) | (beta < 1)) {  # Checks parameter validity
    return("Invalid parameter values.")
  } else {
    b <- NULL  # Initialize vector of samples
    i <- 0  # Initialize counter of accepted samples
   
     # if both parameters > 1, the C chosen can be computed by maximizing the pdf
    if((alpha > 1) & (beta > 1)) {
      max.x <- (1-alpha)/(2-beta-alpha) 
      C <- pdf.beta(max.x, alpha, beta)
    } else {
      # if one of the parameters is 1, C is set to the value of the one that is > 1.
      C <- max(alpha, beta)
    }
    
    while(i < n) { # keep sampling until we have n accepted
      x <- runif(1)  # Generate one sample from proposition
      y <- runif(1, min=0, max=C)  # Generate from U(0, C)
      # accept x if y <= f(x), reject otherwise
      if(y <= pdf.beta(x, alpha, beta)){
        b <- append(b, x)  # append the accepted x into sample vector b
        i <- i + 1
      }
    }
    return(b)
  }
}
        
beta21 <- rbeta.rs(1000,2,1)  # Get 1000 beta(2,1) samples
plot(density(beta21))


## ;;
## -------------------------------------------
## Q3: -- add your code  below ;;
## -------------------------------------------
## ;;

## 3.1

# Initialize start date and time
start <- ymd_hms("2007-10-11 00:00:00", tz="Etc/GMT+0")


# Create new column for what hour accidents occurred at
# floor is used to ensure that, say, an occurence at 0:49 will be in our 0.
accidents.w.hours <- accident.data %>%
  mutate(hour = hour(start + hours(floor(Timehour))))

# Group data by hour and count how many entries are in each group
hour.frequency.data <- accidents.w.hours %>%
  group_by(hour) %>%
  summarize(frequency = n())

# Produces a pdf with a time-like plot of frequency per hour
pdf(file="s1782016.circ.hist.pdf", height=10, width=10)
p <- ggplot(data=hour.frequency.data,
            aes(x = hour , y = frequency)) +
  coord_polar(theta = "x", start = 0 , direction = 1 ) +
  geom_bar(stat = "identity")
p
dev.off()

## 3.2

llik.nhpp <- function(beta, rawdata){
  
  beta <- matrix(beta, nrow=1)  # Converts the list beta into a 1x4 matrix
  
  # Compute integral part first, spanning all hours til last accident
  T.max <- rawdata %>% floor %>% max  # Recover hour of last accident
  j <- seq(1,T.max)  # Generate sequence of hours from 1 til T
  
  j.hour <- j %% 24  # Converts j into 24-hours format
  j.wday <- ((j/24 + 5) %% 7) %>% floor  # Converts j into 7-day format, with Sunday = 1
  j.month <- month(start+hours(j)) # Converts j to 12-mo. format, with Jan = 1
  
  # Create a matrix such that each row will be the "seasonal" indicator
  j.ind <- matrix(1, ncol=T.max)  # First row is 1 for beta in intensity
  j.ind <- rbind(j.ind,  # 2nd row is 1 if hour is between 6am and 6pm
                 as.numeric((j.hour >= 6) & (j.hour < 18)))
  j.ind <- rbind(j.ind,  # 3rd row is 1 if day is monday to friday
                 as.numeric((j.wday >= 2) & (j.wday <= 6)))
  j.ind <- rbind(j.ind,  # 4th row is 1 if month is march to april
                 as.numeric((j.month >= 3) & (j.month <=4)))
  
  # Compute log intensity function for all j
  j.log.int <- beta %*% j.ind
  integral.part <- j.log.int %>% exp %>% sum  # Sum of (non-log) intensities
  
  # Compute the product part (which is also just a sum in a log framework)
 
  n <- length(rawdata)  # Get total number of accidents n
  data.hours <- hours(floor(rawdata))  # Converts timehour into actual time objects
  
  d.hour <- (rawdata %% 24) %>% floor # Retrieves hour in the day of each accident
  d.wday <- ((rawdata/24 + 5) %% 7) %>% floor  # Retrieves weekday of each accident (Sun = 1)
  d.month <- month(start + data.hours)  # Retrieves month of each accident (Jan = 1)
  
  # Create a matrix such that each row will be the seasonal indicator of actual data
  data.ind <- matrix(1, ncol=n) # These rows are the same in defn to j.ind
  data.ind <- rbind(data.ind,  # Add indicator row for hour in day
                    as.numeric((d.hour >= 6) & (d.hour < 18)))
  data.ind <- rbind(data.ind,  # Add indicator row for weekday
                    as.numeric((d.wday >= 2) & (d.wday <= 6)))
  data.ind <- rbind(data.ind,  # Add indicator for month
                    as.numeric((d.month >= 3) & (d.month <= 4)))
  
  #Compute log intensity for each accident occurence
  d.log.int <- beta %*% data.ind 
  product.part <- sum(d.log.int)  # Product part simplifies to the sum of log int.

  negllik <- integral.part - product.part  # Compute the negative log likelihood
  return(negllik)
}

# Use optim to find MLA for parameters beta
optimized.beta <- optim(c(0,0,0,0), fn = llik.nhpp,
                        rawdata = accident.data$Timehour,
                        control = c(maxit=500))

# Save the optimal parameters in beta.hat
beta.hat <- optimized.beta$par

# A function to obtain the values of the intensity function over
# input hours t, starting at start.date

intensity.function <- function(beta, start.date, timehours){
  date.info <- start.date + hours(floor(timehours))  # Convert time hours into actual date object
  n <- length(timehours)  # Gets number of observations needed
  # Obtain info necessary for indicators
  sample.hour <- hour(date.info)
  sample.wday <- wday(date.info)
  sample.month <- month(date.info)
  
  # Create the big indicator matrix
  sample.ind <- matrix(1, ncol=n) # These rows are the same in defn to j.ind
  sample.ind <- rbind(sample.ind,  # Add indicator row for hour in day
                    as.numeric((sample.hour >= 6) & (sample.hour < 18)))
  sample.ind <- rbind(sample.ind,  # Add indicator row for weekday
                    as.numeric((sample.wday >= 2) & (sample.wday <= 6)))
  sample.ind <- rbind(sample.ind,  # Add indicator for month
                    as.numeric((sample.month >= 3) & (sample.month <= 4)))
  
  # Compute Intensity
  sample.log.int <- beta %*% sample.ind  # log intensity function
  sample.intensity <- (exp(sample.log.int))  # intensity function

  return(sample.intensity)
}

# Over the last year - if today is November 5, 2017
start.date <- ymd_hms("2016-11-05 00:00:00", tz="Etc/GMT+0")
max.days <- 365*4  # Assume 365 days in a year and we want a val every 6 hours
sample.t <- 6*seq(0, max.days)  # Create a timehour like sequence, every 6 hrs

# Call function and store the intensity function into a variable
# Note that this gets the intensity function twice a day
fitted.intensity <- intensity.function(beta.hat, start.date, sample.t)

# Plot a graph on a pdf
pdf(file="s1782016.nhpp.intensity.pdf")
plot(x = (start.date + hours(sample.t)), y=fitted.intensity,
     main="Fitted Intensity over previous year, every 6 hours",
     xlab="Date", ylab="Intensity")
dev.off()


## 3.3

# Intensity Function

# Let this be the function to get probability of k incidents in a time pd
# defined by "time.start" and "time.end", both as date objects. k is also
# a vector if natural numbers to be considered.
prob.up.to.k <- function(beta, k, time.start, time.end){
  # Get total time period being considered (-1 so we don't spill over)
  total.hours <- floor((time_length(time.end - time.start))/3600)-1
  hr.seq <- seq(0:total.hours)  # Get sequence of hrs being considered
  
  date.info <- start.date + hours(hr.seq)
  n <- length(hr.seq)
  seq.k <- seq(0:k) - 1
  
  sample.hour <- hour(date.info)
  sample.wday <- wday(date.info)
  sample.month <- month(date.info)
  
  sample.ind <- matrix(1, ncol=n) # These rows are the same in defn to j.ind
  sample.ind <- rbind(sample.ind,  # Add indicator row for hour in day
                      as.numeric((sample.hour >= 6) & (sample.hour < 18)))
  sample.ind <- rbind(sample.ind,  # Add indicator row for weekday
                      as.numeric((sample.wday >= 2) & (sample.wday <= 6)))
  sample.ind <- rbind(sample.ind,  # Add indicator for month
                      as.numeric((sample.month >= 3) & (sample.month <= 4)))
  
  sample.log.int <- beta %*% sample.ind  # log intensity function
  sample.intensity <- (exp(sample.log.int))  # intensity function
  int.intensity <- sum(sample.intensity) # Integral as the sum
  
  # Evaluate the probability expression
  prob.til.k <- exp(-int.intensity) * (int.intensity)^(seq.k) / factorial(seq.k)

  return(prob.til.k)
}

# Define a list of start and end dates based on requirements
sd1 <- ymd_h("2017-12-04 09")  # Dember 4, 2017 is a monday
ed1 <- ymd_h("2017-12-04 15")
sd2 <- ymd_h("2017-12-02 09")  # December 2, 2017 is a sunday
ed2 <- ymd_h("2017-12-02 15")
sd3 <- ymd_h("2016-04-08 18")  # April 8, 2016 is a Friday
ed3 <- ymd_h("2016-04-08 20")
# April 1-7, 2018 is the first week of the month
sd4 <- ymd_h("2018-04-01 00")  
ed4 <- ymd_h("2018-04-08 00")

# Get k=1..9 estimates for each of the defined periods above
prob1 <- prob.up.to.k(beta.hat, 9, sd1, ed1)
prob2 <- prob.up.to.k(beta.hat, 9, sd2, ed2)
prob3 <- prob.up.to.k(beta.hat, 9, sd3, ed3)
prob4 <- prob.up.to.k(beta.hat, 9, sd4, ed4)

# Combine the array of probability estimates into a single list
prob.estimates <- list(prob1, prob2, prob3, prob4)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##*FEEDBACK*:
## Excellent answer and effort to exercise overall.
## One minor is that each element of the list object
## has to be in tabular form
## try e.g.,
## > prob.estimates <-
## >   prob.estimates %>% lapply(as.data.frame)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## ;;
## -------------------------------------------
## Q4: -- add your code below;;
## -------------------------------------------
## ;;


## 4.1

llik.gdp <- function(param, x, thresh, tol.xi.limit = 5*10^(-2)) {
  x.filt <- x[x > thresh]  # Define a vector of x's greater than the threshold
  n <- length(x.filt)  # Define n as the number of excesses
  # param will be a 2-d vector containing (sigma, xi)
  
  if(abs(param[2]) < tol.xi.limit){  # Check if xi is close to 0
    
    # Negative LLik when xi approaches 0; note that when xi -> 0, GDP is exponential
    nllik <- n*param[1] + (1/param[1])*sum(x.filt-thresh)
      
  } else {
    
    # Negative log Lik when xi isn't 0
    inner <- 1 + param[2]*((x.filt-thresh)/param[1])  # dummy variable
    nllik <- n*param[1] + (1+(1/param[2]))*sum(inner[inner>0])  # get only positive part
 
  }
  return(nllik)  # The function returns the negative log-likelihood
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FEEDBACK: log-likelihood is not correct:
#          nllik <- n*param[1] + (1/exp(param[1]))*sum(x.filt-thresh)
#         inner <- log(1 + param[2]*((x.filt-thresh)/param[1])) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 4.2

# Define a range of thresholds for which we will fit the model
u <- quantile(burlington$Precipitation, seq(0.85, 0.99, by=0.01))

# Procedure optimizes the parameters for each threshold defined in u
mle.u <- NULL  # initialize the matrix that will contain the parameters

# Get best parameters for each u and add it to the list
for(i in 1:length(u)){
  new.parameters <- optim(c(1,1), fn=llik.gdp,
                          x = burlington$Precipitation, thresh=u[i])
  mle.u <- rbind(mle.u, new.parameters$par)
}

# Initialize PDF file for saving
pdf(file="s1782016.stability.pdf", height=10, width=20)
param.plots <- par(mfrow=c(1,2))  # Initializes an area for two plots
# Plot sigma vs u
plot(x=u, y=mle.u[,1],  
     type="l", lwd=2, main="MLE of Sigma per Threshold",
     xlab="Threshold", ylab="Sigma")
# Plot xi vs u
plot(x=u, y=mle.u[,2],
     type="l", lwd=2, main="MLE of Xi per Threshold",
     xlab="Threshold", ylab="Xi")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FEEDBACK: produces wrong estimates due to the mistakes in defining llik.gpd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 4.3

ret.level <- function(x, thresh, ret.period){
  param <- optim(c(1,1), fn=llik.gdp,  # Gets MLE for parameters given u
                 x = x, thresh=thresh)$par
  # if zeta is an est. if P(X > threshold), we use a simple proportion est.
  zeta <- length(x[x > thresh])/length(x)
  
  # Formula for return level x_T 
  return.level <- thresh +
    (param[1]/param[2])*((ret.period*zeta)^param[2] - 1)
  
  return(return.level)
}

spring.data <- burlington %>%  # Retrieve Spring Data
  filter(Spring==1)
spring.prec.data <- spring.data$Precipitation # Retrieve Precipitation Column

T <- c(100, 500, 1000, 10000)  # Define return periods

# Initialize a matrix that will contain the returns for each row of u
return.matrix <- NULL

# Initilizize PDF for plotting
pdf(file="s1782016.ret.levels.pdf")
# We obtain returns and plots for different thresholds u
for(i in 1:length(u)){
  temp.return <- ret.level(spring.prec.data, u[i], T)
  return.matrix <- rbind(return.matrix, temp.return)
  plot(x=log(T), y=temp.return, type="l", lwd=2,
       main=paste("Spring day return levels given threshold", u[i], sep=" "),
       xlab="Log of Period", ylab="Return Level")
}
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FEEDBACK: produces wrong estimates due to the mistakes in defining llik.gpd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





