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
  models[[i]] <- lm(X[[i]]$cd4 ~ X[[i]]$smoke + X[[i]]$age + X[[i]]$precd4 )
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

# Plotting coefficients varying over time

par(mfrow = c(2, 2))

# bandwith
band <- 1.5

plot(times, intercepts.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "Intercept Coefficient")
lines(ksmooth(times, intercepts.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times, intercepts.coeff, df = 7), col = "green")

plot(times, smoking.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "Smoking Coefficient")
lines(ksmooth(times, smoking.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times, smoking.coeff, df = 7), col = "green")

plot(times, age.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "Age Coefficient")
lines(ksmooth(times, age.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times,  age.coeff, df = 7), col = "green")

plot(times, precd4.coeff, pch = 19, col = "blue",
     xlab = "Time",
     ylab = "PreCD4 Coefficient")
lines(ksmooth(times, precd4.coeff, "normal", bandwidth = band), col = "red", lwd = 2)
lines(smooth.spline(times, precd4.coeff, df = 7), col = "green")




