library(tidyverse)
library(MASS)


data <- read.table("Bloodpressure.txt", header=T)

# Recover the data
recovery <- data$recovery_time
lqty <- data$log_drug
bp <- data$blood_pressure
lrec <- log(recovery)

data$log_bp <- log(bp)
data$qty2 <- lqty^2
data$log_rec <- log(recovery)

nd <- data.frame(log_drug = 2,
                 blood_pressure = 75,
                 log_bp = log(75),
                 qty2 = 4)
nd

nd.c <- data.frame(log_drug = 2 - mean(lqty),  # Mean corrected new data
                   blood_pressure = 75 - mean(bp),
                   log_bp = log(75)- mean(log(bp)))


# Plots #

M.bp <- lm(recovery_time ~ bp,
           data = data)
M.qty <- lm(recovery_time ~ log_drug,
            data = data)

par(mfrow=c(1,2))
plot(lqty, recovery, pch=20,  # Drug Qty
     main = "(a) Drug Qty. vs. Recovery Time",
     xlab = "log(Drug Quantity (mg))",
     ylab = "Recovery Time (mins)")
lines(lqty, M.qty$fitted.values,
      type="l", col = "red")


plot(bp, recovery, pch = 20,  # Mean BP
     main = "(b) Mean BP vs. Recovery Time",
     xlab = "Mean systolic BP during Hypotension (mm Hg)",
     ylab = "Recovery Time (mins)")
lines(bp, M.bp$fitted.values,
      type="l", col = "red")

summary(M.bp)

summary(M.qty)


# Between the two independents

lm.cov <- lm(blood_pressure ~ log_drug,
             data = data)
summary(lm.cov)


par(mfrow=c(1,1))
cor(bp, lqty)
plot(lqty, bp, pch = 20, main = "Mean Blood Pressure vs. Drug Qty.",
     ylab = "Mean Blood Pressure during Hypotension",
     xlab = "Log (Drug Dosage Volume)")
lines(lqty, lm.cov$fitted.values,
      col = "red", type = "l")

chisq.test(bp, lqty)  # We can't reject independence

# Initial Model #
# Model 1 : linear model qty + bp

M1 <- lm(recovery_time ~ log_drug + blood_pressure,
         data = data)
summary(M1)

predict(M1, newdata = nd,
        interval = "predict") %>% round(4)

hist(M1$residuals, freq = FALSE,
     col = "gray")

layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
plot(M1)

# Model 2 : Stabilized (Mean Corrected) #
data.c <- data.frame(log_drug = lqty - mean(lqty),
                     blood_pressure = bp - mean(bp),
                     recovery_time = recovery,
                     log_rec = log(recovery) - mean(log(recovery)),
                     log_bp = log(bp)-mean(log(bp)))
M2 <- lm(recovery_time ~ log_drug + blood_pressure,
         data = data.c)
summary(M2)

predict(M2, newdata = nd.c,
        interval = "pred")

# Residuals of the model M2 #
hist(M2$residuals, freq = FALSE,
     col = "gray")



# Individual Analysis #






# Exploring Qty variable
plot(lqty, recovery)  # normal plot

nqty <- sqrt(lqty)
plot(nqty, recovery)  # normal plot


# Exploring BP variable
plot(bp, recovery) 
plot(log(bp), recovery)  # try this to squeeze the ranges together

plot(bp^-1/4, recovery)

# Exploring the two independents again



plot(lqty, bp, main = "Scatter Plot of Covariates",
     pch = 20, 
     xlab = "Mean Blood Pressure",
     ylab = "log Drug Qty.")
abline(h = median(bp), v=median(lqty), lwd = 1.5, col = "green")
abline(h = quantile(bp, c(0.25,0.75)),
       v = quantile(lqty, c(0.25, 0.75)),
       col = c("blue", "violet"),
       lwd = 1)
legend(72, 1.8,
       c("25% Quantile", "Median", "75% Quantile"),
       lty = c(1,1,1),
       col = c("blue", "green", "violet"))






# Model 3 : Log all independents #
M3 <- lm(recovery_time ~ log_drug + log_bp,
         data = data)
summary(M3)

predict(M3, newdata = nd,
        interval = "pred") %>% round(4)


# Removing outliers, perhaps? #


par(mfrow = c(2,2))
boxplot(recovery, main = "(a) Boxplot of Recovery Times")
boxplot(lqty, main = "(b) Boxplot of log Drug Dosage")
boxplot(bp, main = "(c) Boxplot of Mean BP")




# Model 4: Removing outliers in recovery time #
par(mfrow=c(1,1))
rec.bp <- boxplot(recovery, main = "Boxplot of Recovery Times")  # Recovery time boxplot

data.squeeze <- filter(data, recovery_time < 60)  # Removed outliers

M4 <- lm(recovery_time ~ log_drug + blood_pressure,  # Note that nothing remained significant
         data = data.squeeze)
summary(M4)
predict(M4, newdata = nd,  # Narrower P.I., but model itself is no longer significant at 0.05 confidence
        interval = "pred") %>% round(4)

# Model 5: what about outliers in the covariates
boxplot(lqty)  # no outliers
bp.bp <- boxplot(bp)  # There are outliers, let's try to remove those points
bp.bp$out

data.sq2 <- filter(data, blood_pressure < 83)

M5 <- lm(recovery_time ~ log_drug + blood_pressure,  # Drug qty remained significant
         data = data.sq2)
summary(M5)
predict(M5, newdata = nd,  # wider P.I. than before
        interval = "pred")  %>% round(4)

# Model 6 : What about quadratic qty #
M6 <- lm(recovery_time ~ qty2 + blood_pressure,  # Drug qty 2
         data = data)

summary(M6)
predict(M6, newdata = nd,
        interval = "pred")


# Model 7 : Log Recovery with origincal covariates #
M7 <- lm(log_rec ~ log_drug + blood_pressure,
         data = data)
summary(M7)

predict(M7, newdata = nd, interval = "pred") %>% round(4)
predict(M7, newdata = nd, interval = "pred") %>% exp %>% round(4)


# Model 8 : log Everything #
M8 <- lm(log_rec ~ log_drug + log_bp,
         data = data)
summary(M7)
predict(M7, newdata = nd,
        interval = "pred") %>% exp


# Model 9 : Centered everything, log rec #
M8 <- lm(log_rec ~ log_drug + blood_pressure,
         data = data.c)
summary(M8)
predict(M8, newdata=nd.c, interval="predict") 
