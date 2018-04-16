#Before using Markov simulated chains to obtain Monte Carlo estimates,
#we should ask whether Markov chain converge to stationary distribution

#A trace plot shows the history of a parameter value across iterations of the chain. 
# It shows you precisely where the chain has been exploring. 
# If the chain is stationary, 
# it should not be showing any long-term trends.
# The average value for the chain, should be roughly flat.
#It appears that this trace plot is tracing the same distribution.
source("Metroplis-Hasting.R",echo=FALSE) # run source function
set.seed(61)
post0 = mh(n=n, ybar=ybar, n_iter=10e3, mu_init=0.0, cand_sd=0.9)
coda::traceplot(as.mcmc(post0$mu[-c(1:500)]))


# If a chain is wandering, the step size of the random walk sampler is too small, 
# and so it take many iterations for the chain to traverse across the distribution. 
set.seed(61)
post1 = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=0.04)
coda::traceplot(as.mcmc(post1$mu[-c(1:500)]))

#If run the chain for many many iterations, it can converge
set.seed(61)
post2 = mh(n=n, ybar=ybar, n_iter=100e3, mu_init=0.0, cand_sd=0.04)
coda::traceplot(as.mcmc(post2$mu))

#Autocorrelation is a number between negative 1 and positive 1 which measures how linearly 
#dependent the current value of the chain is to past values called lags. 
#Sampling 1000 iterations from a highly correlated Markov chain yields 
#less information about the stationary distribution than we would obtain from 
#1000 samples independently drawn from the stationary distribution.
library("coda")
autocorr.plot(as.mcmc(post0$mu))
# This is the same as previous 2 line code
coda::autocorr.plot(as.mcmc(post0$mu))
coda::autocorr.diag(as.mcmc(post0$mu))

coda::autocorr.plot(as.mcmc(post1$mu))
coda::autocorr.diag(as.mcmc(post1$mu))

#The Monte Carlo effective sample size is how many
#independent samples from the stationary distribution 
#you would have to draw to have equivalent information in your Markov chain. 
str(post2) # contains 100,000 iterations
#?effectiveSize
coda::effectiveSize(as.mcmc(post2$mu)) # effective sample size of ~350
## thin out the samples until autocorrelation is essentially 0. This will leave you with approximately independent samples. 
#The number of samples remaining is similar to the effective sample size.
coda::autocorr.plot(as.mcmc(post2$mu), lag.max=500)
# From the autocorr.plot, we know after 400 iterations, the autocorrelation between chain values =0
thin_interval = 400 # how far apart the iterations are for autocorrelation to be essentially 0.
thin_indx = seq(from=thin_interval, to=length(post2$mu), by=thin_interval)
head(thin_indx)

# draw two plots at the same canvas
par(mfrow = c(2,1))
traceplot(as.mcmc(post2$mu))
post2mu_thin = post2$mu[thin_indx]
traceplot(as.mcmc(post2mu_thin))


coda::autocorr.plot(as.mcmc(post2mu_thin), lag.max=10)#autocorrelation  = 0 
coda::effectiveSize(as.mcmc(post2mu_thin))#the effective sample size:250  because the values are approximately uncorrelated. 
length(post2mu_thin) #actual sample length: 250
str(post0) # contains 10,000 iterations
coda::effectiveSize(as.mcmc(post0$mu)) # effective sample size of ~2,500
#It is usually a good idea to check the Monte Carlo effective sample size of your chain when doing mcmc. 
#check the Monte Carlo effective sample size of your chain using the Raftery and Lewis diagnostic.
?raftery.diag
raftery.diag(as.mcmc(post0$mu))
raftery.diag(as.mcmc(post0$mu), q=0.005, r=0.001, s=0.95)


