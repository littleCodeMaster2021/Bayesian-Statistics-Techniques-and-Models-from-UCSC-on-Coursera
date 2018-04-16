#posterior distribution is small so we will work on log scale be more numerically stable
#estimate 
lg = function(mu, n, ybar) {
mu2 = mu^2
n * (ybar * mu - mu2 / 2.0) - log(1 + mu2)
}

#execute the Random-Walk Metropolis-Hastings sampler with normal proposals.
mh = function(n, ybar, n_iter, mu_init, cand_sd) {
  ## step 1, initialize
  mu_out = numeric(n_iter)
  accpt = 0
  mu_now = mu_init
  lg_now = lg(mu=mu_now, n=n, ybar=ybar)
  
  ## step 2, iterate
  for (i in 1:n_iter) {
    ## Draw a candidate sample from a proposal distriution mu ~ N(mu,cand_sd^2)
    mu_cand = rnorm(n=1, mean=mu_now, sd=cand_sd) 
    
    ## compute ratios:
    lg_cand = lg(mu=mu_cand, n=n, ybar=ybar) # evaluate log of g with the candidate
    lalpha = lg_cand - lg_now # log of acceptance ratio
    alpha = exp(lalpha)
    
    ## step 2c
    u = runif(1) #draw a uniform variable which will be less than alpha with probability min(1, alpha)
    # if 0<alpha <1 then mu_now = mu_cand and increase acceptance rate 
    if (u < alpha) { # then accept the candidate
      mu_now = mu_cand
      accpt = accpt + 1 # to keep track of acceptance
      lg_now = lg_cand
    }
    
    ## collect results
    mu_out[i] = mu_now # save this iteration's value of mu
  }
  
  ## return a list of output
  list(mu=mu_out, accpt=accpt/n_iter)
  #mu=mu.out. That contains all our samples from the Metropolis-Hastings sampler
  #also return the accept rate
}

#set the problem:
y = c(-0.2, -1.5, -5.3, 0.3, -0.8, -2.2)
ybar = mean(y)
n = length(y)
hist(y, freq=FALSE, xlim=c(-1.0, 3.0)) # histogram of the data
curve(dt(x=x, df=1), lty=2, add=TRUE) # prior for mu: density function
points(y, rep(0,n), pch=1) # individual data points
points(ybar, 0, pch=19) # sample mean

#run the sampler: use m=1000m=1000 iterations and proposal standard deviation
set.seed(43) # set the random seed for reproducibility
#posterior sampling
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=1)
#We usually are looking for an acceptance rate between 0.23 and 0.5 when we're working with a random walk Metropolis-Hastings algorithm. 
str(post)#what in the object
post$accpt#the acceptance rate should be in [0.23,0.5]
mean(post$mu)
#draw a trace plot which shows the history of the chain
#provides basic feedback about whether the chain has reached its stationary distribution.
#install.packages("coda")

library("coda")
traceplot(as.mcmc(post$mu))


traceplot(as.mcmc(post$mu))
# the acceptance rate is too high (above 50%). Let's try something in between
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=0.9)
post$accpt
traceplot(as.mcmc(post$mu)) #this is good!
#let's see what happens if we initialize the chain at some far-off value.
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=30.0, cand_sd=0.9)
post$accpt
traceplot(as.mcmc(post$mu))# pass the mcmc project to traceplot

#plot the posterior density against the prior to see
#how the data updated our belief about ??
post$mu_keep = post$mu[-c(1:100)] # discard the first 200 samples
str(post)
#a density estimate of the posterior distribution based on the 10 data points.
plot(density(post$mu_keep, adjust=2.0), main="", xlim=c(-1.0, 3.0), xlab=expression(mu)) # plot density estimate of the posterior
curve(dt(x=x, df=1), lty=2, add=TRUE) # prior for mu
points(ybar, 0, pch=19) # sample mean
curve(0.017*exp(lg(mu=x, n=n, ybar=ybar)), from=-1.0, to=3.0, add=TRUE, col="blue") # approximation to the true posterior in blue

#############################################################
#one iteration of an algorithm to simulate a chain 
#whose stationary distribution is p(theta)->g(theta)
# draw candidate
theta_cand = rnorm(n=1, mean=0.0, sd=10.0)

# evaluate log of g with the candidate
lg_cand = lg(theta=theta_cand)

# evaluate log of g at the current value
lg_now = lg(theta=theta_now)

# evaluate log of q at candidate
lq_cand = dnorm(theta_cand, mean=0.0, sd=10.0, log=TRUE)

# evaluate log of q at the current value
lq_now = dnorm(theta_now, mean=0.0, sd=10.0, log=TRUE)

# calculate the acceptance ratio
lalpha = lg_cand + lq_now - lg_now - lq_cand 
alpha = exp(lalpha)

# draw a uniform variable which will be less than alpha with probability min(1, alpha)
u = runif(1)

if (u < alpha) { # then accept the candidate
  theta_now = theta_cand
  accpt = accpt + 1 # to keep track of acceptance
}

#######Independent Metropolis-Hastings (q does not condition on the previous value of the chain) with normal proposal
#Candidates are always drawn from the same N(0,10^2) distribution.
