#Q1

v=2
tau = sqrt(2/v)

# Q2
library("car")
data("Anscombe")
head(Anscombe)
?Anscombe

#When using an informative variable selection prior like the Laplace
#beta perimeter ~ ddexp(0,1.0)
#we typically center and scale the data:

Xc = scale(Anscombe, center=TRUE, scale=TRUE)
str(Xc)
data_jags = as.list(data.frame(Xc))



#. Fit the model in JAGS using the Laplace prior with variance 2 for each of the three coefficients, 
#and an inverse gamma prior for the observation variance with effective sample size 1 and prior guess 1.

library(rjags)
mod_string = " model {
  for (i in 1:length(education)) {
    education[i] ~ dnorm(mu[i], prec)
    mu[i] = b[1]*income[i] + b[2]*young[i] + b[3]*urban[i]
  }
  
 # prior: beta
  for (i in 1:3) {
    b[i] ~ ddexp(0.0, 1.0) # has variance 2.0
  }

  #the observation variance: prec

  prec ~ dgamma(1.0/2.0, 1.0*1.0/2.0)
  ## Initial guess of variance based on overall
  ## variance of education variable. Uses low prior
  ## effective sample size. Technically, this is not
  ## a true 'prior', but it is not very informative.
  sig2 = 1.0 / prec
  sig = sqrt(sig2)
} "

mod_q7_string = " model {
    for (i in 1:length(education)) {
        education[i] ~ dnorm(mu[i], prec)
        mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*urban[i]
    }
    
    b0 ~ dnorm(0.0, 1.0/1.0e6)
    for (i in 1:3) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)
    	## Initial guess of variance based on overall
    	## variance of education variable. Uses low prior
    	## effective sample size. Technically, this is not
    	## a true 'prior', but it is not very informative.
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
} "

params = c("b", "sig")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

modq7 = jags.model(textConnection(mod_q7_string), data=data_jags, n.chains=3)
update(modq7, 1e3)

mod_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

modq7_sim = coda.samples(model=modq7,
                       variable.names=params,
                       n.iter=5e3)
modq7_csim = as.mcmc(do.call(rbind, modq7_sim))


## convergence diagnostics
plot(mod_sim, ask=TRUE)
plot(modq7_sim, ask=TRUE)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## calculate DIC
dic = dic.samples(mod, n.iter=1e3)

dicq7= dic.samples(modq7, n.iter=1e3)

summary(mod_sim)
summary(modq7_sim)

# Q5
data("warpbreaks")

#  Re-fit the model with a separate variance for each group.
#For each variance, use an Inverse-Gamma(1/2, 1/2) prior, 
#corresponding to prior sample size 1 and prior guess 1 for each variance.

mod5_string = " model {
  for( i in 1:length(y)) {
    y[i] ~ dnorm(mu[woolGrp[i], tensGrp[i]], prec[woolGrp[i], tensGrp[i]])
  }
  
  for (j in 1:max(woolGrp)) {
    for (k in 1:max(tensGrp)) {
      mu[j,k] ~ dnorm(0.0, 1.0/1.0e6)
      prec[j,k] ~ dgamma(1.0/2.0, 1.0*1.0/2.0)
      sig[j,k] = sqrt(1.0 / prec[j,k])
    }
  }
  
} "

set.seed(83)

data5_jags = list(y=log(warpbreaks$breaks), 
                  woolGrp=as.numeric(warpbreaks$wool), 
                  tensGrp=as.numeric(warpbreaks$tension))

params5 = c("mu", "sig")

mod5 = jags.model(textConnection(mod5_string), data=data5_jags, n.chains=3)
update(mod5, 1e3)

mod5_sim = coda.samples(model=mod5,
                        variable.names=params5,
                        n.iter=5e3)

## convergence diagnostics
plot(mod5_sim)

gelman.diag(mod5_sim)
autocorr.diag(mod5_sim)
effectiveSize(mod5_sim)

summary(mod5_sim)
dic5 = dic.samples(mod5, n.iter=1e3)

