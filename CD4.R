install.packages("timereg")
library("timereg")

data(cd4)
head(cd4)
tail(cd4)

# Indexes to delete ( we don't consider visit = 0.1, 5.3)
cd4[which(cd4$visit == 0.1), ]$obs # c(524, 799, 896, 1010)
cd4[which(cd4$visit == 5.3), ]$obs # c(131, 486, 1272, 1533)

cd4 <- cd4[-c(524, 799, 896, 1010, 131, 486, 1272, 1533 ), ]


# Centralising Age and PreCD4

cd4$age <- cd4$age - mean(cd4$age)
cd4$precd4 <- cd4$precd4 - mean(cd4$precd4)

dim(cd4) # 1817 10

# Vector containing the single ID's
n.ID <- unique(cd4$id)
length(n.ID) # 283

# Considering visit as time points
plot(cd4$visit[1:7], cd4$cd4[1:7], type = "o", 
     xlim = c(min(cd4$visit), 7),
     ylim = c(min(cd4$cd4), 70),
     xlab = "Time",
     ylab = "CD4",
     lwd = 1)

for(i in n.ID) {
  temp <- subset(cd4, id == i)
  lines(temp$visit, temp$cd4, type = "o")
}

# Different times
times <- sort(unique(cd4$visit))

# Creating a list. Each component, containg subset of the data
# for each different time point.

X <- list()

for(i in 1:length(times)) {
  X[[i]] <- subset(cd4, visit == times[i])
}

# Creating different regressions for each time 
models <- list()
for(i in 1:length(times)) {
  models[[i]] <- lm(X[[i]]$cd4 ~ X[[i]]$smoke + X[[i]]$age + X[[i]]$precd4)
}


# Saving coefficients for each variable along time

intercepts.coeff <- c()
smoking.coeff <- c()
age.coeff <- c()
precd4.coeff <- c()

for(i in 1:length(times)) {
  intercepts.coeff[i] <- models[[i]]$coefficients[1]
  smoking.coeff[i] <- models[[i]]$coefficients[2]
  age.coeff[i] <- models[[i]]$coefficients[3]
  precd4.coeff[i] <- models[[i]]$coefficients[4]
}

# Extracting standard errors

stderr.intercepts.coeff <- c()
stderr.smoking.coeff<- c()
stderr.age.coeff <- c()
stderr.precd4.coeff <- c()

for(i in 1:length(times)) {
  out <- summary(models[[i]])
  stderr.intercepts.coeff[i] <- out$coefficients[1, 2]
  stderr.smoking.coeff[i] <- out$coefficients[2, 2]
  stderr.age.coeff[i] <- out$coefficients[3, 2]
  stderr.precd4.coeff[i] <- out$coefficients[4, 2]
}



# Plotting coefficients varying over time

par(mfrow = c(2, 2))

# bandwith
band <- 1.5


### INTERCEPT
plot(times, intercepts.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "Intercept Coefficient")
lines(ksmooth(times, intercepts.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times, intercepts.coeff, df = 7), col = "green")
# C.I
lines(smooth.spline(times, intercepts.coeff + 2*(stderr.intercepts.coeff), df = 7), 
      col = "grey", lty = 2 )
lines(smooth.spline(times, intercepts.coeff - 2*(stderr.intercepts.coeff), df = 7), 
      col = "grey", lty = 2 )


### SMOKING COEFFICIENT
plot(times, smoking.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "Smoking Coefficient")
lines(ksmooth(times, smoking.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times, smoking.coeff, df = 7), col = "green")
# C.I
lines(smooth.spline(times, smoking.coeff + 2*(stderr.smoking.coeff), df = 7), 
      col = "grey", lty = 2 )
lines(smooth.spline(times, smoking.coeff - 2*(stderr.smoking.coeff), df = 7), 
      col = "grey", lty = 2 )

# AGE COEFFICIENT
plot(times, age.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "Age Coefficient")
lines(ksmooth(times, age.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times,  age.coeff, df = 7), col = "green")
# C.I
lines(smooth.spline(times, age.coeff + 2*(stderr.age.coeff), df = 7), 
      col = "grey", lty = 2 )
lines(smooth.spline(times, age.coeff - 2*(stderr.age.coeff), df = 7), 
      col = "grey", lty = 2 )

# PREC CD4
plot(times, precd4.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "PreCD4 Coefficient")
lines(ksmooth(times, precd4.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times, precd4.coeff, df = 7), col = "green")
# C.I
lines(smooth.spline(times, precd4.coeff + 2*(stderr.precd4.coeff), df = 7), 
      col = "grey", lty = 2 )
lines(smooth.spline(times, precd4.coeff - 2*(stderr.precd4.coeff), df = 7), 
      col = "grey", lty = 2 )




# Saving splined regression coefficient

spline.intercept <- smooth.spline(times, intercepts.coeff, df = 7)$y
spline.smoke <- smooth.spline(times, smoking.coeff, df = 7)$y
spline.age <- smooth.spline(times, age.coeff, df = 7)$y
spline.preCD4 <- smooth.spline(times, precd4.coeff, df = 7)$y

# Creating matrix prediction for each variable, for each time step

prediction.Y <- matrix(NA, nrow = length(n.ID), ncol = length(times))

for(i in 1:length(n.ID)) {
  for(t in 1:length(times)) {
    prediction.Y[i, t] <- spline.intercept[t] + spline.smoke[t] * cd4$smoke[i] + 
      spline.age[t] * cd4$age[i] + spline.preCD4[t] * cd4$precd4[i]
  }
}

# Comparing real value and predicted values

par(mfrow = c(1, 2))

# Actual
plot(cd4$visit[1:7], cd4$cd4[1:7], type = "l", 
     xlim = c(min(cd4$visit), 6),
     ylim = c(min(cd4$cd4), 70),
     xlab = "Time",
     ylab = "CD4",
     lwd = 1)

for(i in n.ID) {
  temp <- subset(cd4, id == i)
  lines(temp$visit, temp$cd4, type = "l")
}

# Predicted
plot(times, prediction.Y[1, ], type = "l", 
     xlim = c(min(cd4$visit), 6),
     ylim = c(min(cd4$cd4), 70),
     xlab = "Time",
     ylab = "CD4",
     lwd = 1)
for(i in 2:length(n.ID)) {
  lines(times, prediction.Y[i, ])
}



# Other visulisation of the same (overlapped on the same graph)

par(mfrow = c(1, 1))

plot(cd4$visit[1:7], cd4$cd4[1:7], type = "l", 
     xlim = c(min(cd4$visit), 6),
     ylim = c(min(cd4$cd4), 70),
     xlab = "Time",
     ylab = "CD4",
     lwd = 1)

for(i in n.ID) {
  temp <- subset(cd4, id == i)
  lines(temp$visit, temp$cd4, type = "l")
}

for(i in 2:length(n.ID)) {
  lines(times, prediction.Y[i, ], col = "red", lwd = 2)
}

