require('timereg')
require('splines')
require('lattice')
require('mgcv')
require('fields')
require('plot3D')
require('eqs2lavaan')
require('ggplot2')

get_data <- function(standardise = F){
  data(cd4)
  data<-cd4
  X <- data.frame(id = double(),
                  visit = double(),
                  age = double(),
                  precd4 = double(),
                  cd4 = double())
  
  ids <- unique(data$id)
  ages <- rep(0, length(ids))
  precd4s <- rep(0, length(ids))
  
  for (i in seq_along(ids)){
    id = ids[i]
    guy <- data[data$id == id,2:7]
    ages[i] <- guy$age[1]
    precd4s[i] <- guy$precd4[1]
    if (length(unique(guy$visit)) != length(guy$visit)){
      n = length(unique(guy$visit))
      last.visit <- guy$visit[n]
      mean.last.visit <- mean(guy$cd4[(guy$visit == last.visit)])
      guy <- guy[1:n,]
      guy$cd4[n] <- mean.last.visit
    }
    X <- as.data.frame(rbind(X,guy))
  }
  
  row.names(X)<- NULL
  #X$smoke <- as.factor(X$smoke)
  X$age <- X$age - mean(ages)
  X$precd4 <- X$precd4 - mean(precd4s)
  if (standardise){
    X$age <- X$age/sd(X$age)
    X$precd4 <- X$precd4/sd(X$precd4)
    X$cd4 <- (X$cd4 - mean(X$cd4))/sd(X$cd4)
  }
  return(X)
}

spaghetti_plot<- function(X){
  plot.new()
  plot.window(xlim = range(X$visit), ylim = range(X$cd4))
  axis(side = 1, pos = -0.1)
  axis(side = 2, pos = -0.1)
  
  for (id in X$id){
    guy <- X[X$id == id,]
    color = ifelse(guy$smoke[1] == 0, 'black', 'red')
    lines(guy$visit, guy$cd4, col = color)
  }
}

find_betas <- function(X, tolerance = 0){
  times <- sort(unique(X$visit))
  times <- times[(times != 5.3) & (times != 0.1)]
  betas = matrix(NA,nrow = 0, ncol = 4)
  for (i in seq_along(times)){
    time <- times[i]
    if (tolerance == 0){
      dat <- X[X$visit == time,]
    }
    else{
      dat <- X[X$visit >= time - tolerance & X$visit <= time + tolerance,]
    }
    model.temp <- lm(cd4 ~ smoke + age + precd4, data = dat)
    betas <- rbind(betas,model.temp$coefficients)
  }
  return(list(times=times, betas=betas))
}

smoothed_fits <- function(betas, times, df = c(4,4,4,4), cv = F){
  par(mfrow = c(2,2))
  betsm = list()
  names = c('Intercept', 'Smoking', 'Age', 'PreCD4')
  for (i in 1:4){
    plot(times, betas[,i], ylab = paste(names[i], 'coefficient'))
    if (cv){
      spl <- smooth.spline(x = times, y = betas[,i],
                           all.knots = T, cv=T)
      print(paste('CV degrees of freedom: ', spl$df))
    }
    else {
      spl <- smooth.spline(x = times, y = betas[,i],
                           all.knots = T, df = df[i])
    }
    betsm[[i]]<-predict(spl,1:59/10)$y
    lines(times, predict(spl,times)$y, col = 'red')  
  }
  par(mfrow = c(1,1))
  betas_mat <- cbind(betsm[[1]], betsm[[2]], betsm[[3]], betsm[[4]])
  return(betas_mat)
}

hat.matrix <- function(X){
  X.t = t(X)
  P <- X %*% solve(X.t %*% X, X.t)
  return(P)
}

isequal <- function(A,B){
  return(ifelse(A==B,1,0))
}

gamma <- function(Xj, Xk, j.equals.k = F, accuracy = 0){
  d = 4
  if (j.equals.k){
    X.j <- cbind(1,Xj[,3:5])
    P.j <- hat.matrix(X.j)
    n.j <- nrow(X.j)
    e.j <- (diag(n.j) - P.j) %*% Xj[,6]
    return(sum(e.j * e.j)/(n.j - d))
  }
  else{
    X.j <- cbind(1,Xj[,3:5])
    X.k <- cbind(1,Xk[,3:5])
    P.j <- hat.matrix(X.j)
    P.k <- hat.matrix(X.k)
    n.j <- nrow(X.j)
    n.k <- nrow(X.k)
    M <- outer(Xj[,1], Xk[,1], isequal)
    e.j <- (diag(n.j) - P.j) %*% Xj[,6]
    e.k <- (diag(n.k) - P.k) %*% Xk[,6]
    temp1 <- sum(diag((diag(n.k) - P.k) %*% t(M) %*% t(diag(n.j)-P.j)))
    temp2 <- sum(diag(e.j %*% t(e.k)))
    if (sum(M) > accuracy & temp1 != 0) return(temp2/temp1)
    else return(NA)
  }
}

covBetas <- function(Xj, Xk, j.equals.k = F){
  if (j.equals.k){
    X.j <- cbind(1,Xj[,3:5])
    Inv.j <- solve(t(X.j) %*% X.j)
    gam <- gamma(Xj, Xk, j.equals.k)
    return(gam * diag(Inv.j))
  }
  else{
    X.j <- cbind(1,Xj[,3:5])
    X.k <- cbind(1,Xk[,3:5])
    pseudoInv.j <- solve(t(X.j) %*% X.j, t(X.j))
    pseudoInv.k <- solve(t(X.k) %*% X.k, t(X.k))
    M <- outer(Xj[,1], Xk[,1], isequal)
    gam <- gamma(Xj, Xk)
    return(gam * diag(pseudoInv.j %*% M %*% t(pseudoInv.k)))
  }
}

cov_matrices <- function(X, times, tolerance = 0){
  # order: Intercept, Smoke, Age, PreCD4
  timesExt <- c(times[1], times, times[length(times)])
  n <- length(timesExt)
  covs <- list()
  for (i in 1:4) covs[[i]]<-matrix(0, nrow = n, ncol = n)
  
  for (j in seq_along(timesExt)){
    if (tolerance == 0){
      data.j <- as.matrix(X[X$visit == timesExt[j],])
    }
    else{
      data.j <- as.matrix(X[X$visit <= timesExt[j]+tolerance & X$visit >= timesExt[j]-tolerance,])
    }
    
    for (k in seq_along(timesExt)){
      
      if (tolerance == 0){
        data.k <- as.matrix(X[X$visit == timesExt[k],])
      }
      else{
        data.k <- as.matrix(X[X$visit <= timesExt[k] + tolerance & X$visit >= timesExt[k] - tolerance,])
      }
      temp <- covBetas(data.j, data.k)
      for (i in 1:4){
        covs[[i]][j,k] <- temp[i]
      }
    }
  }
  return(covs)
}

gamma_matrix <- function(X, times, tolerance = 0, accuracy = 0){
  # order: Intercept, Smoke, Age, PreCD4
  n <- length(times)
  covar <- matrix(NA, nrow = n, ncol = n)
  
  for (j in seq_along(times)){
    if (tolerance == 0){
      data.j <- as.matrix(X[X$visit == times[j],])
    }
    else{
      data.j <- as.matrix(X[X$visit <= times[j]+tolerance & X$visit >= times[j]-tolerance,])
    }
    
    for (k in seq_along(times)){
      
      if (tolerance == 0){
        data.k <- as.matrix(X[X$visit == times[k],])
      }
      else{
        data.k <- as.matrix(X[X$visit <= times[k] + tolerance & X$visit >= times[k] - tolerance,])
      }
      covar[j,k] <- gamma(data.j, data.k, j == k, accuracy)
    }
  }
  return(covar)
}

reg_mat <- function(X, caps = c(-2,2)){
  X[is.na(X)]<-0
  X[X == Inf] <- 0
  X[X == -Inf] <- 0
  X[X >= caps[2]] <- caps[2]
  X[X <= caps[1]] <- caps[1]
  return(X)
}

plot_covariances <- function(X, times, tolerance = 0, bandwidth = .25, caps = c(-1e5,1e5), include_scatter = F){
  names <- c("intercept", "smoke", "age", "PreCD4")
  par(mfrow = c(2,2))
  covs <- cov_matrices(X, times, tolerance)
  for (i in 1:4){
    m<- prep_smooth(covs[[i]][2:58, 2:58])
    look <- smooth.2d(m$resp, x = m$coord, theta = bandwidth)
    M <- mesh(look$x, look$y)
    response <- reg_mat(look$z, caps)
    surf3D(M$x, M$y, response, bty = "b2")
    title(paste("Smoothed covariance for", names[i]))
    if (include_scatter){
      points3D(m$coord$x, m$coord$y, m$resp, add = T, colkey = list(plot=F), col = "black")
    }
  }
  par(mfrow = c(1,1))
}

correl <- function(covar){
  u <- diag(covar)^(-0.5)
  return(t(t(u * covar) * u))
}

prep_smooth <- function(X){
  z <- c()
  x <- c()
  y <- c()
  for (i in seq_len(nrow(X))){
    for (j in seq_len(ncol(X))){
      if (!is.na(X[i,j])){
        x <- c(x, i)
        y <- c(y, j)
        z <- c(z, X[i,j])
      }
    }
  }
  return(list(coord = data.frame(x = x, y = y),
              resp = z))
}

plot_gamma <- function(X, times, tolerance = 0, bandwidth = .25, caps = c(-100,100), include_scatter = F, accuracy = 0){
  m <- gamma_matrix(X, times, tolerance, accuracy)
  temp <- prep_smooth(m)
  look <- smooth.2d(temp$resp, x = temp$coord, theta = bandwidth)
  M <- mesh(look$x, look$y)
  response <- reg_mat(look$z, caps)
  surf3D(M$x, M$y, response, bty = "b2")
  title("Smoothed covariance for residuals")
  if (include_scatter){
    points3D(temp$coord$x, temp$coord$y, temp$resp, add = T, colkey = list(plot=F), col = "black")
  }
}

plot_robust_gamma <- function(X, times, tolerance = 0, bandwidth = .25, caps = c(-100,100), include_scatter = F, accuracies = c(0,2,4,6)){
  par(mfrow = c(2,2))
  for (i in 1:4){
    plot_gamma(X, times, tolerance, bandwidth, caps, include_scatter, accuracies[i])
  }
  par(mfrow = c(1,1))
}

point_var <- function(covs, times){
  vars <- list()
  for (i in 1:4) vars[[i]] <- rep(0, length(times))
  
  for (i in seq_along(times)){
    N <- matrix(bs(times[i], knots = times[2:(length(times)-1)],Boundary.knots =c(0.2,5.9), intercept = T), ncol = 1)
    for (j in 1:4){
      covs[[j]]<- reg_mat(covs[[j]], caps = c(-1e17, 1e17))
      vars[[j]][i] <- t(N) %*% covs[[j]] %*% N
    }
  }
  return(vars)
}

prep_CI_smooth <- function(times, vars){
  if (any(is.na(vars))){
    ind <- which(is.na(vars))
    return( list(t = times[-ind], v = vars[-ind]))
  }
  else{
    return(list(t = times, v = vars))
  }
}

plot_with_CI <- function(betas, times, vars, df = c(4,4,4,4), cv = F, add_kernel = F){#
  names = c('Intercept', 'Smoking', 'Age', 'PreCD4')
  par(mfrow = c(2,2))
  for (i in 1:4){
    plot(times, betas[,i], ylab = paste(names[i], 'coefficient'))
    if (cv){
      spl <- smooth.spline(x = times, y = betas[,i],
                           all.knots = T, cv = T)
    }
    else {
      spl <- smooth.spline(x = times, y = betas[,i],
                           all.knots = T, df = df[i])
      print(paste('CV degrees of freedom: ', spl$df))
    }
    predictions <- predict(spl,times)$y
    lines(times, predictions, col = 'red')
    lines(times, predictions + 2 * sqrt(vars[[i]]))
    lines(times , predictions - 2 * sqrt(vars[[i]]))
    temp <- prep_CI_smooth(times, 2 * sqrt(vars[[i]]))
    smoothed_CI <- predict(smooth.spline(x = temp$t, y = temp$v, df = 3), times)$y
    lines(times, predictions + smoothed_CI, lty = 2)
    lines(times, predictions - smoothed_CI, lty = 2)
    if (add_kernel){
      lines(ksmooth(times, betas[,i],"normal", bandwidth = 1.5), col = 'green')
    }
  }
  par(mfrow = c(1,1))
}

plot_with_CI(betas, times, vars, cv = T, add_kernel = T)

subsample_indiv <- function(X){
  ids <- unique(X$id)
  resampled_ids <- sample(ids, length(ids), replace = T)
  matrix_ids <- c()
  for (id in resampled_ids){
    matrix_ids <- c(matrix_ids, which(X$id == id))
  }
  return(X[matrix_ids,])
}

impute_NAs <- function(X){
  RM <-rowMeans(X, na.rm = T)
  n <- nrow(X)
  for (i in 1:n){
    X[i,is.na(X[i,])] <- RM[i]
  }
  return(X)
}

bootstrap_CI <- function(X, times.len, iter){
  betas <- list()
  for (i in 1:4){
    betas[[i]] <- matrix(NA, nrow = times.len, ncol = iter)
  }
  n_low <- floor(0.025 * iter)
  n_high <- floor(0.975 * iter)
  
  for (i in 1:iter){
    X_new <- subsample_indiv(X)
    temp <- find_betas(X_new, tolerance = 0)
    for (j in 1:4){
      betas[[j]][,i] <- temp$betas[,j]
    }
  }
  betas_low <- matrix(NA, nrow = times.len, ncol = 4)
  betas_high <- matrix(NA, nrow = times.len, ncol = 4)
  for (i in 1:4){
    temp <- impute_NAs(betas[[i]])
    temp <- t(apply(temp, 1, sort))
    betas_low[,i] <- temp[,n_low]
    betas_high[,i] <- temp[,n_high]
  }
  return(list(low = betas_low, high = betas_high))
}

prep_boot_CI_smooth <- function(X, caps = c(-80,80)){
  X[X < caps[1]] <- caps[1]
  X[X> caps[2]] <- caps[2]
  return (X)
}

plot_boot_CI <- function(betas, times, boot_CI, df = c(4,4,4,4), cv = F){
  names = c('Intercept', 'Smoking', 'Age', 'PreCD4')
  par(mfrow = c(2,2))
  for (i in 1:4){
    plot(times, betas[,i], ylab = paste(names[i], 'coefficient'))
    if (cv){
      spl <- smooth.spline(x = times, y = betas[,i],
                           all.knots = T, cv = T)
    }
    else {
      spl <- smooth.spline(x = times, y = betas[,i],
                           all.knots = T, df = df[i])
    }
    predictions <- predict(spl,times)$y
    lines(times, predictions, col = 'red')
    
    lines(times, boot_CI$low[,i])
    lines(times , boot_CI$high[,i])
    if (i == 3) caps = c(-10,10)
    else if (i == 4) caps = c(-2, 2)
    else caps = c(-80,80)
    smoothed_CI_low <- predict(smooth.spline(x = times,
                                             y = prep_boot_CI_smooth(boot_CI$low[,i], caps),
                                             df = 3),
                               times)$y
    smoothed_CI_high <- predict(smooth.spline(x = times,
                                              y = prep_boot_CI_smooth(boot_CI$high[,i], caps),
                                              df = 3),
                                times)$y
    lines(times, smoothed_CI_low, lty = 2)
    lines(times, smoothed_CI_high, lty = 2)
  }
  par(mfrow = c(1,1))
}

split_var <- function(X, test){
  X <- X[test(X),]
  time <- X[,2]
  y <- X[,6]
  X <- cbind(1,X[,3:5])
  return(list(X=X, y=y, time = time))
}

calculate_resid <- function(data,betas){
  time_ind <- data$time * 10
  residuals <- rep(0,nrow(data$X))
  for (i in seq_len(nrow(data$X))){
    residuals[i] <- data$y[i] - sum(data$X[i,] * betas[time_ind[i],])
  }
  return(data.frame(resid = residuals, time=data$time))
}

concatenate_groups <- function(groups, lab){
  data <- matrix(NA, ncol = 1, nrow = 1)
  counter = 1
  for (g in groups){
    g$group <- lab[counter]
    counter <- counter + 1
    if (is.na(data[1,1])) data <- g
    else{
      data <- rbind(data, g)
    }
  }
  data <- as.data.frame(data)
  #print(factor(data$group, levels = c('A', 'B')))
  #data$group <- factor(data$group, levels = 1:(counter -1),labels = lab)
  return(data)
}

plot_residuals <- function(X, betas, var_split="smoke"){
  groups <- list()
  if (var_split == "smoke"){
    func <- list(function(x) return(x$smoke == 1),
                 function(x) return(x$smoke == 0))
    n <- 2
    lab <- c('smokers', 'non-smokers')
  }
  else if (var_split == "age"){
    q <- quantile(unique(X$age), c(0.25, 0.5, 0.75))
    func <- list(function(x) return(x$age <= q[1]),
                 function(x) return(q[1] < x$age & x$age <= q[2]),
                 function(x) return(q[2] < x$age & x$age <= q[3]),
                 function(x) return(x$age > q[3]))
    n <- 4
    lab <- c(paste('age <=',round(q[1],2)),
             paste(round(q[1],2),'< age <=',round(q[2],2)),
             paste(round(q[2],2),'< age <=',round(q[3],2)),
             paste(round(q[3],2), '< age'))
  }
  else if (var_split == "preCD4"){
    q <- quantile(unique(X$precd4), c(0.25, 0.5, 0.75))
    func <- list(function(x) return(x$precd4 <= q[1]),
                 function(x) return(q[1] < x$precd4 & x$precd4 <= q[2]),
                 function(x) return(q[2] < x$precd4 & x$precd4 <= q[3]),
                 function(x) return(x$precd4 > q[3]))
    n <- 4
    lab <- c(paste('preCD4 <=',round(q[1],2)),
             paste(round(q[1],2),'< preCD4 <=',round(q[2],2)),
             paste(round(q[2],2),'< preCD4 <=',round(q[3],2)),
             paste(round(q[3],2), '< preCD4'))
    
  }
  for (i in 1:n){
    temp <- split_var(X, func[[i]])
    groups[[i]] <- calculate_resid(temp, betas)
  }
  final <- concatenate_groups(groups, lab)
  p <- ggplot(transform(final, group=factor(group,levels = lab)), aes(factor(time), resid))
  p + geom_boxplot() + geom_jitter() + facet_grid(group ~.)
}

X <- get_data(F)
spaghetti_plot(X)
temp <- find_betas(X, tolerance = 0)
betas <- temp$betas
gam <- gamma_matrix(X, times, 0)
betas_smoothed <- smoothed_fits(betas, times, cv = T)
covs <- cov_matrices(X, times, tolerance = 0)
vars <- point_var(covs, times)
plot_with_CI(betas, times, vars, df = c(4,4,4,4))
plot_with_CI(betas, times, vars, cv = T)
# next line might need to be run more than once, sometimes throws errors (depending on sampling)
boot_CI <- bootstrap_CI(X,length(times), 100)
plot_boot_CI(betas, times, boot_CI, cv = T)
plot_gamma(X, times, tolerance = 0, bandwidth = 10, caps = c(-1e5,1e5), include_scatter = T, accuracy = 4)
plot_robust_gamma(X, times, tolerance = 0, bandwidth = 10, caps = c(-1e5,1e5), include_scatter = T, accuracies = c(0,2,4,6))
plot_covariances(X, times, tolerance = 0, bandwidth = 10, include_scatter = F)
plot_residuals(X, betas_smoothed, "smoke")
plot_residuals(X, betas_smoothed, "age")
plot_residuals(X, betas_smoothed, "preCD4")

