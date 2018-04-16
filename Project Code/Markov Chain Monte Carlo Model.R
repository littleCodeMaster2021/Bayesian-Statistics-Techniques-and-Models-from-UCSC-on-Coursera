#rjags for Markov Chain Monte Carlo model: y|mu~ N(mu,1) mu~t(0,1,1) mu is unknown variable

#install.packages("rjags")

#1. Specify the model
library("rjags")
mod_string = " model {
  for (i in 1:n) {
     y[i] ~ dnorm(mu, 1.0/sig2)#In JAGS this is the precision which is the reciprocal of the variance. 
   }
   mu ~ dt(0.0, 1.0/1.0, 1.0) # location, inverse scale, degrees of freedom
   sig2 = 1.0
} "

#2. Set up the model
set.seed(50) # reproduce results
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
n = length(y)
data_jags = list(y=y, n=n)
params = c("mu")
inits = function() {
  inits = list("mu"= 0.0) # set initialized value for each parameter in model, 
  #we can generate random sample: list("mu"= rnorm(1))
} 
# initialize the model
mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits)

#3.runs the MCMC sampler for 500 iterations without saving the samples anywhere.
update(mod, 500) 
# keep the simulation continue what last off
mod_sim = coda.samples(model=mod,variable.names=params,n.iter=1000)

#4.Post processing
#evaluate the Markov chains we've simulated to determine if
#they're suitable for inference.  
summary(mod_sim)
library("coda")
plot(mod_sim)
#plot the posterior distribution for mu

