## 
##
## -------------------------------------------
## LOAD PACKAGES
## -------------------------------------------
##
## 
## Load any R package required for your code
## to successfully run!

library(tidyverse)

## ;;
setwd("")

## -------------------------------------------
## READING DATA
## -------------------------------------------
## ;;
load("burlington.RData")

## --------------------------------------------------
## --------------------------------------------------

## ;;
## ---------------------------------------------
## Q1: -- add your code below 
## ---------------------------------------------
## ;;

## 1.1

# Constructs the 95% bootstrap CI using R "replications"
CI.cor <- function(raw, R = 10000, conf = 0.95){ 
  n <- raw$n
  data <- cbind(raw$lsat,  # Compile the data into a 2 column matrix
                raw$gpa)
  rho_star <- NULL  # Initialize vector of covariances
  
  for(r in 1:R){  # For R Replications
    index_star <- 1:n %>%
      sample(size = n, replace = TRUE)  # Get n resamples from data
    data_star <- data[index_star,]  # Obtains the resamples
    rho_star[r] <- cor(data_star)[1,2]  # Gets correl and stores in vector
  }
  
  # Computing estimate of bias and variance for bootstrap rho
  bias_rho_star <- mean(rho_star - cor(data)[1,2])
  var_rho_star <- var(rho_star)
  
  # Define alpha confidence level
  alpha = 1 - conf
  
  # Obtain Normal Bootstrap CI
  CI_normal_low <- cor(data)[1,2] - bias_rho_star +
    qnorm(alpha/2) * sqrt(var_rho_star)
  CI_normal_high <- cor(data)[1,2] - bias_rho_star +
    qnorm(1-alpha/2) * sqrt(var_rho_star)
  
  # Obtain Percentile Bootstrap CI
  CI_low <- quantile(rho_star, alpha/2)
  CI_high <- quantile(rho_star, 1-alpha/2)
  
  # Tabulate the results for output
  results <- data.frame(
    "Lower Bound" = c(CI_normal_low, CI_low),
    "Upper Bound" = c(CI_normal_high, CI_high),
    row.names = c("Normal", "Percentile"))
  
  return(results)
}

# Construct Data
raw.1 <- list(lsat = c(576, 635, 558, 578, 666, 580, 555,  
                     661, 651, 605, 653, 575, 545, 572, 594),
            gpa = c(3.39, 3.30, 2.81, 3.03, 3.55, 3.07, 3.00,
                    3.43, 3.36, 3.13, 3.12, 2.74, 2.76, 2.88, 2.96),
            n = 15)

# Call Function
CI.cor(raw.1)

## 1.2

CI.var.ratio <- function(raw, R = 10000, conf = 0.95){ 
  n <- raw$n
  data <- cbind(raw$lsat,  # Compile the data into a 2 column matrix
                raw$gpa)
  vr_star <- NULL  # Initialize vector of variance ratios
  
  for(r in 1:R){  # For R Replications
    index_star <- 1:n %>%
      sample(size = n, replace = TRUE)  # Get n resamples from data
    data_star <- data[index_star,]  # Obtains the resamples
    vr_star[r] <- var(data_star[,1])/var(data_star[,2])  # Gets var ratios
  }
  
  vr_hat = var(data[,1])/var(data[,2])  # Define for easier referencing
  
  # Computing estimate of bias and variance for bootstrap var ratios
  bias_vr_star <- mean(vr_star - vr_hat)
  var_vr_star <- var(vr_star)
  
  # Define alpha confidence level
  alpha = 1 - conf
  
  # Obtain Normal Bootstrap CI
  CI_normal_low <- vr_hat - bias_vr_star +
    qnorm(alpha/2) * sqrt(var_vr_star)
  CI_normal_high <- vr_hat - bias_vr_star +
    qnorm(1-alpha/2) * sqrt(var_vr_star)
  
  # Obtain Percentile Bootstrap CI
  CI_low <- quantile(vr_star, alpha/2)
  CI_high <- quantile(vr_star, 1-alpha/2)
  
  # Tabulate the results for output
  results <- data.frame(
    "Lower Bound" = c(CI_normal_low, CI_low),
    "Upper Bound" = c(CI_normal_high, CI_high),
    row.names = c("Normal", "Percentile"))
  
  return(results)
}

CI.var.ratio(raw.1)  # Call Function

## ;;
## -------------------------------------------
## Q2: -- add your code below  
## -------------------------------------------
## ;;

# The given code for GDP
fit.gpd <- function(x, thresh, tol.xi.limit=5e-2, ...){
  llik.gpd <- function(par, x, thresh, tol.xi.limit=5e-2)
  {
    y <- x[x>thresh]
    sigma <- exp(par[1])
    xi <- par[2]
    n <- length(y)
    if(abs(xi)<=tol.xi.limit)
    {
      llik <- n*log(sigma) + sum((y-thresh)/sigma)
      return(llik)
    }
    par.log.max <- log( pmax( 1+xi*(y-thresh)/sigma, 0 ) )
    llik <- -n*log( sigma )-(1+( 1/xi ))*sum( par.log.max )
    llik <- -ifelse( llik > -Inf, llik, -1e40 )
    return(llik)
  }
  fit <- optim(par = c(0, 0), fn = llik.gpd,
               x=x, thresh=thresh,
               control=list( maxit=10000, ... ))
  sigmahat <- exp( fit$par[1] )
  xihat <- fit$par[2]
  return(c(sigmahat, xihat))
}

fit.gpd(c(0,0),
        x=burlington$Precipitation,
        thresh=quantile(burlington$Precipitation,0.8))


## 2.1

# Write an inversion sampling function for GDP
rgdp <- function(n, thresh, sigma, xi, tol.xi.limit = 5e-2){
  z <- runif(n)  # Generate uniform RV's
  if(abs(xi) < tol.xi.limit){  # Consider the limiting case
    x <- thresh - sigma*log(1-z)  # Inverse of limiting CDF
    return(x)
  } else {
    x <- thresh + (sigma/xi) * ((1-z)^(-xi) - 1)
    return(x)
  }
}

# Parametric Bootstrap Function
parboot.gpd <- function(data, thresh, R=10000){
  n <- length(data)  # Get length of data
  
  # Fit the data to GDP
  model <- fit.gpd(c(0,0),
                   x=data,
                   thresh=thresh)
  sigma <- model[1]
  xi <- model[2]
  
  # Initialize estimate of p.hat is ratio of "successes" over total trials
  p.hat <- length(data[data > thresh]) / length(data)

  mle <- c(sigma, xi, p.hat)  # Store the MLE Estimators based on data
  
  theta.star <- matrix(ncol=3, nrow=R)  # Intialize matrix to store parameter fits
  p <- p.hat  # Initialize for first replication
  for(r in 1:R){
    n.star <- rbinom(1, n, p)  # Generate binomial random variables to get exceedence counts
    x <- rgdp(n.star, thresh = thresh,
              sigma = sigma, xi = xi)  # Get GPD Samples
    mod.star <- fit.gpd(c(0,0),  # Fit the samples to a pareto
                        x=x,
                        thresh=thresh)
    theta.star[r, 1] <- mod.star[1]  # Store sigma in column 1
    theta.star[r, 2] <- mod.star[2]  # Store xi in column 2
    theta.star[r, 3] <- p
    p <- n.star/n  # Update p for next replication
  }
  
  # Compute Biases
  bias.sigma <- mean(theta.star[,1] - sigma)
  bias.xi <- mean(theta.star[,2] - xi)
  bias.p <- mean(theta.star[,3] - p.hat)
  
  biases <- c(bias.sigma, bias.xi, bias.p)  # Store biases in vector
  
  # Compute standard errors
  se.sigma <- sd(theta.star[,1])
  se.xi <- sd(theta.star[,2])
  se.p <- sd(theta.star[,3])
  
  se <- c(se.sigma, se.xi, se.p)
  
  return(list(mle = mle, bias = biases,  # Return the named list
              se = se, distn = theta.star))
}

answer2.1 <- parboot.gpd(data=burlington$Precipitation,
                         thresh=quantile(burlington$Precipitation,0.9))

## 2.2

# Non-Parametric Bootstrap Function
npboot.gpd <- function(data, thresh, R=10000){
  n <- length(data)  # Get length of data
  
  # Get MLE estimates from original data
  mle.fit <- fit.gpd(c(0,0),
                     x=data,
                     thresh=thresh)
  sigma.hat <- mle.fit[1]
  xi.hat <- mle.fit[2]
  p.hat <- length(data[data > thresh]) / length(data)
  mle <- c(sigma.hat, xi.hat, p.hat)  # Store the MLE Estimators based on data
  
  theta.star <- matrix(ncol=3, nrow=R)  # Intialize matrix to store parameter fits
  
  for(r in 1:R){
      # Get a sample from the original data
      data.star <- data %>%
        sample(size = n, replace = TRUE)
      # Fit the model to the sample data
      fit.star <- fit.gpd(c(0,0),
                          x = data.star,
                          thresh = thresh)
      theta.star[r,1] <- fit.star[1]  # Store sigma
      theta.star[r,2] <- fit.star[2]  # Store xi
      theta.star[r,3] <- length(data.star[data.star > thresh]) / n  # Store p
  }
  
  # Compute Biases
  bias.sigma <- mean(theta.star[,1] - sigma.hat)
  bias.xi <- mean(theta.star[,2] - xi.hat)
  bias.p <- mean(theta.star[,3] - p.hat)
  
  biases <- c(bias.sigma, bias.xi, bias.p)  # Store biases in vector
  
  # Compute standard errors
  se.sigma <- sd(theta.star[,1])
  se.xi <- sd(theta.star[,2])
  se.p <- sd(theta.star[,3])
  
  se <- c(se.sigma, se.xi, se.p)
  
  return(list(mle = mle, bias = biases,  # Return the named list
              se = se, distn = theta.star))
}

answer2.2 <- npboot.gpd(data=burlington$Precipitation,
                         thresh=quantile(burlington$Precipitation,0.9))

## 2.3

thresh = quantile(burlington$Precipitation,0.9)  # Define for future use

# Non-parametric bootstrap for returns and 95% CI
ret.level <- function(T, param = FALSE, thresh = quantile(burlington$Precipitation,0.9),
                   alpha = 0.05){
  
  if(param == TRUE) {
    gpd <- answer2.1  # Get parametric data
  } else {
    gpd <- answer2.2  # Get non-parametric data
  }
  
  # Obtain distribution using bootstrap distributions of GPD
  distn <- gpd$distn
  r.distn <- thresh + (distn[,1]/distn[,2]) * ((T*distn[,3])^distn[,2] - 1)

  # Compute Percentile Confidence Intervals
  CI_lower <- quantile(r.distn, alpha/2)
  CI_upper <- quantile(r.distn, 1-alpha/2)
  r.CI <- c(CI_lower, CI_upper)
  
  return(list(distn = r.distn,
              CI = r.CI))
}

# Obtain parametric bootstrap for different T's
r.p.100 <- ret.level(100, param = TRUE)
r.p.500 <- ret.level(500, param = TRUE)
r.p.1000 <- ret.level(1000, param = TRUE)
r.p.10000 <- ret.level(10000, param = TRUE)

# Obtain non-parametric bootstraps for different T's
r.np.100 <- ret.level(100)
r.np.500 <- ret.level(500)
r.np.1000 <- ret.level(1000)
r.np.10000 <- ret.level(10000)

## ;;
## -------------------------------------------
## Q3: -- add your code  below ;;
## -------------------------------------------
## ;;

## 3.1

## Posteriors in PDF file

## 3.2

# Initialize data
pumps.data <- list(t = c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5),
                   x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22),
                   n = 10)

# Define acceptance probability for alpha
accept.alpha <- function(a0, a1, V=1, n=10, beta=1, theta=rep(1,10)) {
  return( (gamma(a0)/gamma(a1))^n *
           beta^(n*(a1-a0)) *
            prod(theta) ^ (a1-a0) *
            exp(-(a1-a0)) *
            a1/a0 )
}

# Generate samples from the posterior using MCMC
# Ndraw is the number of draws, nburn is how many to discard among first entries,
# V is the variance of random walk on log a
mcmc.pumps <- function(data, ndraw=100, nburn=0, # Default initial parameters
                       theta0 = rep(0.5,10), beta0 = 1,
                       alpha0 = 0.5, V = 1) {  
  
    # Retrieve values from data
  n <- data$n
  t <- data$t
  x <- data$x
  
  # The parameters differ from each pump to pump, i think
  d.theta <- matrix(nrow = ndraw, ncol=n)  # Intialize a matrix to store draws for theta
  d.others <- matrix(nrow = ndraw, ncol=2)  # Initialize matrix for beta and alpha
  iter <- -nburn  # For burning entries
  
  # Set initial parameters as the current parameters
  theta <- theta0
  alpha <- alpha0 
  beta <- beta0
  alpha.accept.count <- 0
  
  while(iter < ndraw) {
    iter <- iter + 1  # Increase counter
    
    # Update thetas
    for(i in 1:n) {  # Update done for each pump in turn
      # Update using conditional distributions
      theta[i] <- rgamma(1, shape = x[i] + alpha0, rate = beta0 + t[i])
    }

    # Update beta
    beta <- rgamma(1, shape = n * alpha, rate = sum(theta) + 0.01)
    
    # Generate proposal alpha
    alpha.prop <- exp(log(alpha) + rnorm(1,0,V))  # Random walk on log alpha
    alpha.accept <- min(1, accept.alpha(a0 = alpha,  # Get acceptance criterion
                                        a1 = alpha.prop,
                                        V = V,
                                        n = n,
                                        beta = beta,
                                        theta = theta))
    if(runif(1) <= alpha.accept){  # If unif is below the acceptance criterion
      alpha <- alpha.prop  # accept alpha
      if(iter > 0) {
        alpha.accept.count <- alpha.accept.count + 1  # Count accepted alphas
      }
    }
    
    # Record draws if burning is done
    if(iter > 0) {
      d.theta[iter,] <- theta
      d.others[iter,1] <- beta
      d.others[iter,2] <- alpha
    }
  }
  
  # Return Values
  return( list(theta = d.theta,
               parameters = d.others[,1:2],
               acceptance = alpha.accept.count / ndraw))
}

answer3.2 <- mcmc.pumps(pumps.data)

## 3.3

# Predictive distribution for failure for ith pump
predictive.pumps <- function(i, t, pumps.data = answer3.2, max.x = 30){
  # Retrieve samples from posterior distribution
  theta <- pumps.data$theta[,i]

  x <- 1:max.x  # Maximum number of failures to compute distribution for
  
  pred.dist <- NULL
  
  for(f in 1:length(x)){  # Prob of X failures, X=0,1,2,...,30
    pred.dist[f] <- mean(dpois(f, lambda=theta*t))  # Predictive distribution
  }
  
  return(cbind("Failures" = x,  # Return a table that lists the prob for x failures
               "Probability" = pred.dist))  
}

# Call the function for pump 1 and 94.3 time units for up to 25 failures
answer3.3 <- predictive.pumps(1, 94.3, max.x = 25)
plot(answer3.3)  # Compare with expected probabilities based on data


## ;;
## -------------------------------------------
## Q4: -- add your code below;;
## -------------------------------------------
## ;;

## 4.1
accept.mu <- function(x, mu0, mu1){
  terms <- ((10+(x-mu1)^2)/(10+(x-mu0)^2))^(-11/2)
  first.part <- prod(terms)
  second.part <- exp( -(1/2) * (mu1^2 - mu0^2))
  return(first.part * second.part)
}

# MCMC for t-distribution using random walk
# V is the sd of the normal candidate generator
mcmc.t <- function(x, ndraw=1000, nburn=100, mu0=0, v=1){
  mu <- mu0  # Set initial mu as current
  iter <- -nburn  # Initialize for data discarding
  draws <- NULL  # Initialize vector of draws
  
  while(iter < ndraw){
    iter <- iter + 1  # update counter
    mu.prop <- rnorm(1, mu, v)  # Normal Proposal
    acc.mu <- min(1, accept.mu(x, mu, mu.prop))  # Compute acceptance probability
    
    if(runif(1) <= acc.mu){  # If accepted, change the current mu
      mu <- mu.prop
    }
    
    if(iter > 0){  # Once burning is done
      draws[iter] <- mu
    }
  }
  # Return the draws
  return(draws)
}

# Generate t-distribution samples
x <- rt(12, df = 10)
# Generate samples from posterior distribution for mu
answer4.1 <- mcmc.t(x)

## 4.2

mcmc.gibbs <- function(x, mu0=0, z0=1, ndraw = 1000, nburn=100){
  iter <- -nburn  # Initialize for burning entries
  n <- length(x)
  # Initialize variables
  mu <- mu0
  z <- rep(z0, n)
  
  d.mu <- NULL  # Initialize an empty vetor for mu
  d.z <- matrix(nrow = ndraw, ncol = n)  # Initialize matrix for draws of z
  
  while(iter < ndraw){
    iter <- iter + 1  # Increment counter
    # Update mu and z via conditional posteriors
    mu <- rnorm(1, mean = sum(x*z)/sum(z+1), sd = sqrt(sum(z+1)))
    for(i in 1:n) {
      z[i] <- rgamma(1, shape = 11/2, rate = (1/2)*(x[i]-mu)^2 + 5)
    }
    
    if(iter > 0){  # Record draws if done burning
      d.mu[iter] <- mu
      d.z[iter,] <- z
    }
  }
  # Return named list
  return(list(mu = d.mu,
         z = d.z))
}

# Generate data that comes from t-distribution with 10 df
x <- rt(12, df=10)
plot(density(x))
answer4.2 <- mcmc.gibbs(x)
plot(density(answer4.2$mu))   # Inspect the distribution for mu

## 4.3
# Predictive distribution of t
predictive.t <- function(mu = answer4.2$mu){
  x <- seq(-3,5, len=100)  # points at which distribution will be evaluated
  n <- length(x)
  m <- length(mu)
  dist <- NULL
  
  for(i in 1:n){  # Evaluate the integral at each point
    dist[i] <- mean(dnorm(x[i],
                          mean = mu,
                          sd = sqrt(rgamma(m, shape = 11/2,
                                           rate = (1/2)*(x[i]-mu)^2 + 5))))
  }
  return(cbind(x, dist))
}
    
answer4.3 <- predictive.t()
plot(density(answer4.3))

## 4.4
## Predictive Quantile
predictive.quantile <- function(mu = answer4.2$mu, alpha = 0.05){
  preds <- predictive.t(mu)[,2]  # Runs the prediction distribution and gets the values
  
  # Compute CI bounds
  CI_lower <- quantile(preds, alpha/2)
  CI_upper <- quantile(preds, 1-alpha/2)
  
  return(c(CI_lower, CI_upper))
}

predictive.quantile()

## ------------------------------------------------

## ;;
## -------------------------------------------
## DRAFT
## -------------------------------------------
## ;;
## foo <- rnorm(100)
##
##
## 



