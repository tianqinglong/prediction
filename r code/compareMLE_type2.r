setwd("~/Projects/Prediction/")
dyn.load("main2.so")

# boot=1:FRWB;boot=0:ParaBoot
getType2DataMLE <- function(seed,lower,upper,r,n,eta,beta,boot){
	DataMLE <- .Call("transferDataToR",seed,lower,upper,r,n,eta,beta,boot)
	nr <- DataMLE[1:2]
	B <- (length(DataMLE)-4-nr[1]-8)/2
	mle <- DataMLE[3:4]
	data <- DataMLE[5:(4+nr[1])]
	bootBeta <- DataMLE[(5+nr[1]):(4+B+nr[1])]
	bootEta <- DataMLE[(4+nr[1]+B+1):(4+nr[1]+2*B)]
	intervals <- DataMLE[(4+nr[1]+2*B+1):length(DataMLE)]
	
	output <- list(NR = nr,MLE = mle, Data = data,betaB = bootBeta, etaB = bootEta,
	               lowerIntervals = intervals[1:4], upperIntervals = intervals[5:8])

	return(output)
}

t1 <- Sys.time()
NOS <- 200
p1 <- matrix(nrow = NOS, ncol = 4)
p2 <- matrix(nrow = NOS, ncol = 4)
low = 0.05;up = 0.95;r=20;n=30;eta=1.2;beta=1.5

for(i in 1:NOS){
  fromC1 <- getType2DataMLE(i^2,low,up,r,n,eta,beta,0)
  lowers1 <- fromC1$lowerIntervals;uppers1 <- fromC1$upperIntervals
  fromC2 <- getType2DataMLE(i^2,low,up,r,n,eta,beta,1)
  lowers2 <- fromC2$lowerIntervals;uppers2 <- fromC2$upperIntervals
  
  p1[i,] <- pweibull(uppers1,beta,eta) - pweibull(lowers1,beta,eta)
  p2[i,] <- pweibull(uppers2,beta,eta) - pweibull(lowers2,beta,eta)
}
colMeans(p1);colMeans(p2)
t2 <- Sys.time()
t2-t1

# fromC <- getType2DataMLE(123,0.025,0.975,10,30,1,2)
# fromC$lowerIntervals;fromC$upperIntervals
# 
# # Fonseca's Formula
# fonseca <- function(q, bootBeta, bootEta, mleBeta, mleEta){
#   B <- length(bootBeta)
#   out <- mean(pweibull(qweibull(pweibull(q, mleBeta, mleEta)
#                     ,shape = bootBeta, scale = bootEta),
#            shape = mleBeta, scale = mleEta))
#   return(out)
# }
# 
# # Percentile Bootstrap
# percentileBoot <- function(q, bootBeta, bootEta){
#   out <- mean(pweibull(q, bootBeta, bootEta))
#   return(out)
# }
# 
# # GPQ
# GPQInterval <- function(q, bootBeta, bootEta){
#   muBoot <- log(bootEta)
#   sigmaBoot <- 1/bootBeta
# 
#   muhat <- log(fromC$MLE[2])
#   sigmahat <- 1/fromC$MLE[1]
# 
#   muGPQ <- muhat + (muhat-muBoot)/sigmaBoot*sigmahat
#   sigmaGPQ <- sigmahat/sigmaBoot*sigmahat
#   
#   out <- mean(pweibull(q, 1/sigmaGPQ, exp(muGPQ)))
#   return(out)
# }

#plot(x,y,pch=16,cex=0.5)
#points(x,k,pch=16,cex=0.5,col = "red")
# points(x,z,pch = 16,col="red")

# Data <- fromC$Data
# n=50;r=30
# likehood <- function(para){
# 	beta <- para[1]
# 	eta <- para[2]
# 	loglik <- 0
# 	for(i in 1:r){
# 		loglik <- loglik + dweibull(Data[i], beta, eta, log = TRUE)
# 	}
# 	for(i in (r+1):n){
# 		loglik <- loglik + pweibull(Data[i], beta, eta, lower.tail = FALSE, log.p = TRUE)
# 	}
# 	return(-loglik)
# }
# 
# library(microbenchmark)
# fromC$MLE
# optim(c(1.5,1.5), likehood)
# microbenchmark(optim(c(1.5,1.5), likehood),getType2DataMLE(), times=10000L)
