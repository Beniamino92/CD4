install.packages("timereg")

library("timereg")
library("lattice")

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


# Just a bit of exploratory data anlysis to check if smoking
# has some impact:



par(mfrow = c(1, 2))

# Trajectories "Non-Smoker"
plot(cd4$visit[1:7], cd4$cd4[1:7], type = "o", 
     xlim = c(min(cd4$visit), 6),
     ylim = c(min(cd4$cd4), 70),
     xlab = "Time",
     ylab = "CD4",
     lwd = 1,
     main = "Non-Smoker")

for(i in n.ID) {
  temp.non.smoker<- subset(cd4, id == i & smoke == 0)
  lines(temp.non.smoker$visit, temp.non.smoker$cd4, type = "o")
}

# Trajectories "Smoker"
plot(cd4$visit[29:33], cd4$cd4[29:33], type = "o", 
     xlim = c(min(cd4$visit), 6),
     ylim = c(min(cd4$cd4), 70),
     xlab = "Time",
     ylab = "CD4",
     lwd = 1,
     main = "Smoker"
)

for(i in n.ID) {
  temp.smoker<- subset(cd4, id == i & smoke == 1)
  lines(temp.smoker$visit, temp.smoker$cd4, type = "o")
}



# Let's also look at the box-plot
par(mfrow = c(1, 1))
boxplot(cd4 ~ smoke, data = cd4, col = "pink",
        xlab = "Smoke", ylab = "CD4")


# It doesn't look like there is a significant difference.
# Let's do a t-test to check this.

cd4.not.smokers <- subset(cd4, smoke == 0)
cd4.smokers <- subset(cd4, smoke == 1)

t.test(cd4.not.smokers, cd4.smokers) # p-value = 0.5271







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
band <- 1.5 # for kernel smoothing
df <- 2 # for splines


### INTERCEPT
plot(times, intercepts.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "Intercept Coefficient")
lines(ksmooth(times, intercepts.coeff, "normal", bandwidth = band), 
      col = "red", lwd = 2)
lines(smooth.spline(times, intercepts.coeff, df = df), col = "green")
# C.I
lines(smooth.spline(times, intercepts.coeff + 2*(stderr.intercepts.coeff), df = df), 
      col = "grey", lty = 2 )
lines(smooth.spline(times, intercepts.coeff - 2*(stderr.intercepts.coeff), df = df), 
      col = "grey", lty = 2 )


### SMOKING COEFFICIENT
plot(times, smoking.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "Smoking Coefficient")
lines(ksmooth(times, smoking.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times, smoking.coeff, df = df), col = "green")
# C.I
lines(smooth.spline(times, smoking.coeff + 2*(stderr.smoking.coeff), df = df), 
      col = "grey", lty = 2 )
lines(smooth.spline(times, smoking.coeff - 2*(stderr.smoking.coeff), df = df), 
      col = "grey", lty = 2 )

# AGE COEFFICIENT
plot(times, age.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "Age Coefficient")
lines(ksmooth(times, age.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times,  age.coeff, df = df), col = "green")
# C.I
lines(smooth.spline(times, age.coeff + 2*(stderr.age.coeff), df = df), 
      col = "grey", lty = 2 )
lines(smooth.spline(times, age.coeff - 2*(stderr.age.coeff), df = df), 
      col = "grey", lty = 2 )

# PREC CD4
plot(times, precd4.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "PreCD4 Coefficient")
lines(ksmooth(times, precd4.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times, precd4.coeff, df = df), col = "green")
# C.I
lines(smooth.spline(times, precd4.coeff + 2*(stderr.precd4.coeff), df = df), 
      col = "grey", lty = 2 )
lines(smooth.spline(times, precd4.coeff - 2*(stderr.precd4.coeff), df = df), 
      col = "grey", lty = 2 )




# Saving splined regression coefficient

spline.intercept <- smooth.spline(times, intercepts.coeff, df = df)$y
spline.smoke <- smooth.spline(times, smoking.coeff, df = df)$y
spline.age <- smooth.spline(times, age.coeff, df = df)$y
spline.preCD4 <- smooth.spline(times, precd4.coeff, df = df)$y

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






# Plotting individual trajectories to have more knowledge about
# what the fuck is going on in this fucking patients. My suggestion for them is
# to use a condom and do not inject heroin with shared niddles.

# ID: 1022 - 2334
xyplot(cd4 ~ visit | factor(id), data=cd4[1:246, ], 
       as.table=T, type = c("p", "l"))
# ID: 2385 - 3282
xyplot(cd4 ~ visit | factor(id), data=cd4[247:443, ], 
       as.table=T, type = c("p", "l"))
# ID: 3319 - 4096
xyplot(cd4 ~ visit | factor(id), data=cd4[444:668, ], 
       as.table=T, type = c("p", "l"))
# ID: 4103 - 4870
xyplot(cd4 ~ visit | factor(id), data=cd4[669:852, ], 
       as.table=T, type = c("p", "l"))
# ID: 4874 - 5856
xyplot(cd4 ~ visit | factor(id), data=cd4[853:1009, ],
       as.table=T, type = c("p", "l"))
# ID: 5863 - 7079
xyplot(cd4 ~ visit | factor(id), data=cd4[1010:1216, ], 
       as.table=T, type = c("p", "l"))
# ID: 7143 - 7980
xyplot(cd4 ~ visit | factor(id), data=cd4[1217:1405, ], 
       as.table=T, type = c("p", "l"))
# ID: 7980 - 8893
xyplot(cd4 ~ visit | factor(id), data=cd4[1406:1603, ],
       as.table=T, type = c("p", "l"))
# ID: 8942 - 9954
xyplot(cd4 ~ visit | factor(id), data=cd4[1604:1809, ], 
       as.table=T, type = c("p", "l"))


