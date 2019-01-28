setwd("~/Projects/Prediction/")
dyn.load("main.so")
type1DataMLE <- function(seed, lower, upper, Er, Pt, scale, shape, bootstrap){
  out <- .Call("transToR1",seed, lower, upper, Er, Pt, scale, shape, bootstrap)
  rn <- out[1:2]
  n <- rn[1]
  MLE <- out[3:4]
  data <- out[5:(4+n)]
  B <- (length(out)-2-2-n-8)/2
  betaDraw <- out[(5+n):(4+n+B)]
  etaDraw <- out[(5+B+n):(4+2*B+n)]
  intervals <- out[(4+n+2*B+1):length(out)]
  
  return(list(RN=rn,MLEs = MLE, Data = data, 
              BetaDraws = betaDraw, EtaDraws = etaDraw,
              lowerIntervals = intervals[1:4], upperIntervals = intervals[5:8]))
}

fromC <- type1DataMLE(1, 0.05, 0.95, 10, 0.5, 1.5, 1.5, 1)

t1 <- Sys.time()
NOS <- 5000
p1 <- matrix(nrow = NOS, ncol = 4)
low = 0.05; up = 0.95; Er = 12; Pt = 0.25; eta=1.2; beta=1.5

for(i in 1:NOS){
  fromC1 <- type1DataMLE(i^2,low,up,Er,Pt,eta,beta,1)
  lowers1 <- fromC1$lowerIntervals;uppers1 <- fromC1$upperIntervals
  
  p1[i,] <- pweibull(uppers1,beta,eta) - pweibull(lowers1,beta,eta)
}
colMeans(p1)
t2 <- Sys.time()
t2-t1
