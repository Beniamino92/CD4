#Producing graphs to examine the linearity of the relationship between the response and each predictor
par(mfrow = c(1,3))
for(i in 1:length(times)){
  data<-X[[i]]
  means<-c(0,mean(data$smoke),mean(data$age),mean(data$precd4))
  boxplot(data$cd4~data$smoke,xlab="smoke",ylab="cd4")
  points(seq(1,2,by = 0.01),models[[i]]$coefficients[1]+models[[i]]$coefficients[2]*seq(0,1,by = 0.01)+models[[i]]$coefficients[3]*means[3]+models[[i]]$coefficients[4]*means[4],type="l",col="red")
  plot(data$age,data$cd4,xlab="age",ylab="cd4")
  points(seq(min(data$age),max(data$age),by = 0.1),models[[i]]$coefficients[1]+models[[i]]$coefficients[3]*seq(min(data$age),max(data$age),by = 0.1)+models[[i]]$coefficients[2]*means[2]+models[[i]]$coefficients[4]*means[4],type="l",col="red")
  mtext(paste("t = ",i,sep=""),cex=1.5)
  plot(data$precd4,data$cd4,xlab="precd4",ylab="cd4")
  points(seq(min(data$precd4),max(data$precd4),by = 0.1),models[[i]]$coefficients[1]+models[[i]]$coefficients[4]*seq(min(data$precd4),max(data$precd4),by = 0.1)+models[[i]]$coefficients[3]*means[3]+models[[i]]$coefficients[2]*means[2],type="l",col="red")
  
}