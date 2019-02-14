Weibull.mle <- function(data)
{
  r = data[1]
  n = data[2]
  time = data[-c(1,2)]	

  if(!exists("survreg")) library(survival)
  
  event <- rep(c(1,0), times = c(r,n-r))
  
  sobj <- Surv(time, event, type = "right")
  sfit <- survreg(sobj~1, dist = "weibull", weights = rexp(n,1))
  
  eta = as.numeric(exp(sfit$coef))
  beta = 1/sfit$scale
  return(c(beta, eta))
}
