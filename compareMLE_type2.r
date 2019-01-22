setwd("~/Projects/Prediction/")
dyn.load("main2.so")

getType2DataMLE <- function(){
	DataMLE <- .Call("transferDataToR")
	nr <- DataMLE[1:2]
	mle <- DataMLE[3:4]
	data <- DataMLE[-(1:4)]
	output <- list(NR = nr,MLE = mle, Data = data)

	return(output)
}

fromC <- getType2DataMLE()

Data <- fromC$Data

n=50;r=30
likehood <- function(para){
	beta <- para[1]
	eta <- para[2]
	loglik <- 0
	for(i in 1:r){
		loglik <- loglik + dweibull(Data[i], beta, eta, log = TRUE)
	}
	for(i in (r+1):n){
		loglik <- loglik + pweibull(Data[i], beta, eta, lower.tail = FALSE, log.p = TRUE)
	}
	return(-loglik)
}

library(microbenchmark)
fromC$MLE
optim(c(1.5,1.5), likehood)
microbenchmark(optim(c(1.5,1.5), likehood),getType2DataMLE(), times=10000L)
