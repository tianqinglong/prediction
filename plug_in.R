library(survival)
library(rlist)
library(lattice)

simulate_data <- function(type, r , Pt, shape, scale)
{
  n <- round( r/Pt )
  
  if(type == 1)
  {
    y_vector <- numeric(n)
    x <- 0
    
    for( i in 0 : ( n-1 ) )
    {
      u <- runif(1)
      y <- 1 - (1-x) * ( 1-u )^ ( 1.0 / (n-i) )
      
      if( y > Pt){
        break
      }
      
      y_vector[i+1] <- y
      x <- y
    }
    
    if( sum( y_vector > 0 ) < 2){
      output <- simulate_data(type, r, Pt, shape, scale)
      
      return(output)
    }
    
    num_failure <- sum( y_vector > 0 )
    failure_times <- scale * ( -log( 1 - y_vector[ 1:num_failure ] ) ) ^ ( 1/shape )
    censor_time <- scale * ( -log( 1 - Pt ) ) ^ ( 1/shape )
    
    output <- list(Number_of_Failures = num_failure, Censor_Time = censor_time,
                   Failure_Times = failure_times, Total_Number = n)
    
    return(output)
  }
  
  else if(type == 2)
  {
    
    y_vector <- numeric(r)
    x <- 0
    
    for(i in 0 : (r-1) )
    {
      u <- runif(1)
      
      y <- 1 - (1-x) * ( 1-u )^ ( 1.0 / (n-i) )
      y_vector[ i+1 ] <- y
      
      x <- y
    }
    
    failure_times <- scale * ( -log( 1 - y_vector ) ) ^ ( 1/shape )
    
    return(list(Number_of_Failures = r, Censor_Time = failure_times[r],
                Failure_Times = failure_times, Total_Number = n))
  }
}

get_Weibull_mle <- function(censor_data)
{
  n_minus_r <- censor_data[[4]] - censor_data[[1]]
  
  sobj <- Surv(time = c( censor_data[[3]],
                        rep(censor_data[[2]], n_minus_r)
                         ),
               event = c( rep(1, censor_data[[1]] ),
                          rep(0, n_minus_r)
                          ),
               type = 'right'
               )
  sfit <- survreg(sobj~1, dist = 'weibull')
  
  return( list( Shape = 1/sfit$scale,
                Scale = as.numeric( exp(sfit$coefficients) ),
                Number_of_Failures = censor_data[[1]],
                Censor_Time = censor_data[[2]],
                Failure_Times = censor_data[[3]],
                Total_Number = censor_data[[4]])
          )
}

find_binom_prob <- function(beta, eta, interval, censor_time){
  
  p_binom <- pweibull(censor_time + interval, beta, eta)-
    pweibull(censor_time, beta, eta)
  
  p_binom <- p_binom / (1 - pweibull(censor_time, beta, eta))
  
  return( p_binom )
}

get_Plugin_Coverage_Probability <- function(MLEs, Pt, inter, lower, upper, beta, eta){
  interval <- qweibull(Pt + inter, beta, eta)-
    qweibull(Pt, beta, eta)
  
  binom_prob_real <- find_binom_prob( beta, eta, interval, MLEs[[4]] )
  binom_prob_mle <- find_binom_prob( MLEs[[1]], MLEs[[2]], interval, MLEs[[4]] )
  
  n_minus_r <- MLEs[[6]] - MLEs[[3]]
  
  lowerbound <- 
    max(
      qbinom(lower, n_minus_r, binom_prob_mle) - 1,
      0
      )
  
  upperbound <- 
    qbinom(upper, n_minus_r, binom_prob_mle)
  
  conditional_probability <- 
    pbinom(upperbound, n_minus_r, binom_prob_real)-
    pbinom(lowerbound, n_minus_r, binom_prob_real)+
    dbinom(lowerbound, n_minus_r, binom_prob_real)
  
  return(list.append(MLEs,
                     Lower_Bound = lowerbound,
                     Upper_Bound = upperbound,
                     Conditional_Prob = conditional_probability,
                     Binom_Prob = binom_prob_mle)
         )
}

conditional_probability_discrete <- function(type, r, Pt, beta, eta, inter = 0.1, lower = 0, upper = 0.9)
{
  
  Monte_Carlo_Samples <- list()
  
  for(i in 1:N){
    Monte_Carlo_Samples[[i]] <- simulate_data(type, r, Pt, beta, eta)
  }
  
  Monte_Carlo_MLEs <- lapply(Monte_Carlo_Samples, get_Weibull_mle)
  Monte_Carlo_CP <- lapply(Monte_Carlo_MLEs, get_Plugin_Coverage_Probability,
                           Pt = Pt, inter = inter, lower = lower,
                           upper = upper, beta = beta, eta = eta)
  
  return(Monte_Carlo_CP)
}

