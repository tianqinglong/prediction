dyn.load("main1.so")
type1DataMLE <- function(){
	out <- .Call("transToR1")
	rn <- out[1:2]
	MLE <- out[3:4]
	data <- out[-(1:4)]

	return(list(RN=rn,MLEs = MLE, Data = data))
}

fromC <- type1DataMLE()

n = fromC$RN[1]
r = fromC$RN[2]

Data = fromC$Data

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

fromC$MLE
optim(c(1.5,1.5), likehood)
library(microbenchmark)
microbenchmark(optim(c(1.5,1.5), likehood),type1DataMLE(), times=1000L)