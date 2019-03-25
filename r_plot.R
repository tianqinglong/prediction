# Simulation Setting
type <- 1
beta <- 1
eta <- 1
Pt <- 0.1
N <- 5000
l_prob <- 0
u_prob <- 0.9

# Different Values of r
r_vector <- c(50, 100, 150, 200, 500, 1000, 2000)

# Store the coverage probability for each r
real_cp <- numeric( length(r_vector) )

# Store the conditional coverage probabilities for each r
coverage_probability_matrix <- matrix(ncol = length( r_vector ), nrow = N)

# Change the nominal probability below; "inter" controls the width of the interval
for(i in 1:length( r_vector ) )
{
  r <- r_vector[i]
  out <- conditional_probability_discrete(type, r, Pt, beta, eta, inter = 0.1, lower = l_prob, upper = u_prob)
  
  coverage_probability_matrix[,i] <- sapply( out, function(x) x[[9]] )
  real_cp[i] <- mean( coverage_probability_matrix[,i] )
}

plot(r_vector, real_cp, xlab = "Expected Number of Failures",
     ylab = "Coverage Probability",
     main = paste("Coverage Probability of Plug-in Method with Pt=",Pt),
     ylim = c(0.6 ,u_prob+0.01), pch = 16, cex = 0.75)
lines(r_vector, real_cp, lty = 3)
abline(h = u_prob, lty = 2, col = "red")

histogram(coverage_probability_matrix[,1],
          xlab = "Coverage Probability",
          main = paste("The Expected Number of Failure:",
                       r_vector[1])
)

histogram(coverage_probability_matrix[,10],
          xlab = "Coverage Probability",
          main = paste("The Expected Number of Failure:",
                       r_vector[10])
)

histogram(coverage_probability_matrix[,20],
          xlab = "Coverage Probability",
          main = paste("The Expected Number of Failure:",
                       r_vector[20])
)

histogram(coverage_probability_matrix[,30],
          xlab = "Coverage Probability",
          main = paste("The Expected Number of Failure:",
                       r_vector[30])
)

###

# Simulation Setting
type <- 1
beta <- 1
eta <- 1
Pt <- 0.1
N <- 5000
l_prob <- 0
u_prob <- 0.5

# Different Values of r
r_vector <- seq(5, 200, by = 5)

real_cp <- matrix(ncol = length(r_vector), nrow = 4)

for(pt in c(0.1, 0.2, 0.3, 0.4))
{
  for(i in 1:length( r_vector ) )
  {
    r <- r_vector[i]
    out <- conditional_probability_discrete(type, r, pt, beta, eta, inter = 0.1, lower = l_prob, upper = u_prob)
    
    coverage_probability_matrix <- matrix(ncol = length( r_vector ), nrow = N)
    coverage_probability_matrix[,i] <- sapply( out, function(x) x[[9]] )
    real_cp[10*pt, i] <- mean( coverage_probability_matrix[,i] )
  }
}


plot(r_vector, real_cp[4,], xlab = "Expected Number of Failures",
     ylab = "Coverage Probability",
     main = "Plug-in; Type-I Censoring; Nominal = 0.5",
     ylim = c(.4 ,0.75), pch = 16, cex = 0.75, col = 4)
lines(r_vector, real_cp[4,], lty = 3, col = 4)

points(r_vector, real_cp[3,], pch = 16, cex = 0.75, col = 3)
lines(r_vector, real_cp[3,], lty = 3, col = 3)

points(r_vector, real_cp[2,], pch = 16, cex = 0.75, col = 2)
lines(r_vector, real_cp[2,], lty = 3, col = 2)

points(r_vector, real_cp[1,], pch = 16, cex = 0.75, col = 6)
lines(r_vector, real_cp[1,], lty = 3, col = 6)

abline(h = u_prob, lty = 2, col = 8)

legend("topright", legend = c("Pt = 0.4",
                                 "Pt = 0.3",
                                 "Pt = 0.2",
                                 "Pt = 0.1"),
       col = c(4, 3, 2, 6),
       pch = rep(16,4),
       lty = rep(3,4)
)

####

# Simulation Setting
type <- 1
beta <- 1
eta <- 1
Pt <- 0.1
N <- 10000
l_prob <- 0
u_prob <- 0.9
inter <- 0.1

width <- qweibull(Pt+inter, beta, eta)-
  qweibull(Pt, beta, eta)
cenTime <- qweibull(Pt, beta, eta)

real_binom_prob <- find_binom_prob(beta, eta, width, cenTime)

r_value <- c(50, 100, 150, 200, 500, 1000, 2000)
pdf(file = "~/Desktop/BinomProb.pdf")
for(i in 1:7)
{
  out <- 
    conditional_probability_discrete(type, r = r_value[i],
                                     Pt, beta, eta,
                                     inter = inter,
                                     lower = l_prob,
                                     upper = u_prob)
  hist(sapply(out, function(x) x[[10]]),
       main = paste(
         "The Distribution of the Binomial Probability;r = ",r_value[i]
       ),
       xlab = "Probability", xlim = c(0, 0.25))
  abline(v = real_binom_prob, col = "red", lty = 2)
}
dev.off()

r_vector <- c(5, 10, 20, 50, 80, 100, 150, 200, 500, 1000, 1500, 2000)
pdf(file = "~/Desktop/pmf_of_Binomial_Distribution.pdf")
par(mfrow = c(2,2))
for( i in 1:length(r_vector) ){
  plot( dbinom( 0:(9*r_vector[i]), 9*r_vector[i], real_binom_prob),
        xlab = "Number of Failures",
        ylab = "Probability",
        main = paste("n-r=", 9*r_vector[i]),
        cex = 0.3,
        type = "l",
        xlim = c(0.5*r_vector[i], 2*r_vector[i])
        )
  lines(dbinom( 0:(9*r_vector[i]), 9*r_vector[i], real_binom_prob+0.01),
        col = "red")
}
dev.off()
