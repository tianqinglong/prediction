Weibull.mle <- function(data)
{
  options(warn=-1)
  r = data[1]
  n = data[2]
  times = data[3:(2+n)]
  weights = data[(3+n):(3+n+r)]

  last_weight = weights[r+1]
  weights = c(weights[1:r], rep(last_weight/(n-r), times = n-r))

  wbfuncFRWB <- function(para){
    a = para[1]
    b = para[2]
    
    beta = exp(-a)
    eta = exp(b-exp(a)*log(-log(1-r/2/n)))
    
    value = sum(weights[1:r]*dweibull(times[1:r], beta, eta,log = 1))+
      sum(weights[(r+1):n]*pweibull(times[(r+1):n], beta, eta, lower.tail = 0, log.p = 1))

    return(-value)
  }

  init_sigma <- (log(times[r+1])-log(times[1]))/(log(-log(1-r/n))-log(-log(1-0.5/n)))
  init_a <- log(init_sigma)
  init_eta <- exp(log(times[r+1])-init_sigma*log(-log(1-r/n)))
  init_b <- log(init_eta)+init_sigma*log(-log(1-r/2/n))

  op <- optim(c(init_a, init_b), wbfuncFRWB)
  a <- op$par[1]
  b <- op$par[2]

  beta = exp(-a)
  eta = exp(b-exp(a)*log(-log(1-r/2/n)))

  if(op$convergence != 0){
    print("not coverge")
  }

  return(c(beta, eta))
}