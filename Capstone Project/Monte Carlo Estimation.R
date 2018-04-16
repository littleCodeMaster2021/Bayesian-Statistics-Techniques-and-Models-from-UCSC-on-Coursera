set.seed(32) # Initializes the random number generator so we can replicate these results. 
#To get differentrandom numbers, change the seed.
m = 100 
a = 2.0 
b = 1.0 / 3.0 
#Use rgamma function to simulate m independent samples
theta = rgamma(n=m, shape=a, rate=b) 
#plot a histogram for generated data
hist(theta, freq=FALSE) 
curve(dgamma(x=x, shape=a, rate=b), col="blue", add=TRUE) 
# sample mean 
mean(theta)
a / b # true expected value E(theta)
m = 1e4 #10^4
theta = rgamma(n=m, shape=a, rate=b) 
mean(theta) 
var(theta) # sample variance 
a / b^2 # true variance of Gamma(a,b) 

#plot the distribution
hist(slope, breaks=1000, freq=F, main=main, xlab="slop value(percent)", xlim=c(0,150), ylim=c(0,.05) )
lines(density(slope, bw=1), col="green")

#Beta(5,3) posterior distribution for ??.
#the posterior distribution for ?? and use these samples 
#to approximate the posterior mean for odds of success ( E(??/1-??) ) that is greater than 1
theta = rbeta(9999, 5, 3)
alpha = theta / (1 - theta)
mean( alpha )
head(alpha) # first couple alpha
tail(alpha) # last couple alpha
mean( alpha > 1.0 )

###################approximate the probability that theta <5

ind = theta < 5.0 # set of indicators, TRUE if theta_i < 5 mean(ind) # automatically converts FALSE/TRUE to 0/1 
mean(ind) # automatically converts FALSE/TRUE to 0/1 
quantile(x=theta, probs=0.9) # 90th percentile of theta
qgamma(p=0.9, shape=a, rate=b) # true value of 0.9 quantile 
# we are reasonably confident that the Monte Carlo estimate is no more than this far from the truth
se = sd(theta) / sqrt(m) 
2.0 * se 
#Monte Carlo estimates for the probability that theta <5
# we are reasonably confident that the Monte Carlo estimate is no more than this far from the tru th 
ind = theta < 5.0 
se = sd(ind) / sqrt(m) 
2.0 * se

#Use a (large) Monte Carlo sample to approximate the 0.3 quantile of the 
#standard normal distribution (N(0,1)), the number such that the probability of being less than it is 0.3.
quantile(rnorm(9999, 0.0, 1.0), 0.3)
#check the quantile function:
qnorm(p = 0.3,0.0,1.0)

##################simulating a hierarchical model
#y ??? ?? ~Bin(10,??) and ?? ~ Beta(2,2)
m = 10e4 
y = numeric(m) # create the vectors we will fill in with simulations 
phi = numeric(m) 
for (i in 1:m) { 
  phi[i] = rbeta(n=1, shape1=2.0, shape2=2.0) 
  y[i] = rbinom(n=1, size=10, prob=phi[i]) 
  } 
# which is equivalent to the following 'vectorized' code 
phi = rbeta(n=m, shape1=2.0, shape2=2.0) 
y = rbinom(n=m, size=10, prob=phi) 
mean(y)
plot(prop.table(table(y)), ylab="P(y)", main="Marginal distribution of y") 



