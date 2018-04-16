#fit our hierarhical model for counts of chocolate chips.
rm(list=ls())
dat = read.table(file="cookies.dat", header=TRUE)
head(dat)
# frequency of location
table(dat$location)
# visualize the distribution of chips by location
hist(dat$chips)
boxplot(chips ~ location, data=dat)

#Prior predictive checks
#we can use from the posterior distribution of mu and sigma to
#simulate the posterior predictive distribution of the mean for a new location. 
#select independent exponential priors for alpha and beta
#they are hyperparameters governing the gamma distribution lamda which is the expected number 
#of chocolate chips per cookie. alpha and beta control the distribution of these means between 
#locations. 
#If this is high, the mean number of chips will vary widely from location to location. 

set.seed(112)
n_sim = 500
alpha_pri = rexp(n_sim, rate=1.0/2.0)
beta_pri = rexp(n_sim, rate=5.0)
mu_pri = alpha_pri/beta_pri  #The mean of lamda will represent the overall mean of number of chips for all cookies.
#The variance of this gamma distribution controls the variability between locations.
sig_pri = sqrt(alpha_pri/beta_pri^2)
summary(mu_pri)
summary(sig_pri)
#use those samples of alpha and beta to simulate further down the hierarchy:
lam_pri = rgamma(n=n_sim, shape=alpha_pri, rate=beta_pri)
summary(lam_pri)
#prior predictive reconstruction of the original data set:
(lam_pri = rgamma(n=5, shape=alpha_pri[1:5], rate=beta_pri[1:5]))  
(y_pri = rpois(n=150, lambda=rep(lam_pri, each=30)))

#these priors: lamda have high variance and are somewhat noninformative, 
#they produce unrealistic predictive distributions. 
#Although enough data would overwhelm the prior, resulting in useful posterior distributions. 

###########Another approach would be to re-parameterize the lam_pri(gamma prior) to fit the model
library("rjags")
mod_string = " model {
for (i in 1:length(chips)) {
chips[i] ~ dpois(lam[location[i]])
}

for (j in 1:max(location)) {
lam[j] ~ dgamma(alpha, beta)
}

alpha = mu^2 / sig^2
beta = mu / sig^2

mu ~ dgamma(2.0, 1.0/5.0)
sig ~ dexp(1.0)

} "

set.seed(113)

data_jags = as.list(dat)

params = c("lam", "mu", "sig")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)

##########################Model checking################
#we can check the fit via residuals. 
#With a hierarhcical model, there are now two levels of residuals: the observation level 
#and the location mean level. 
#we’ll look at the residuals associated with the posterior means of the parameters.
#Check the observation residuals, based on the estimates of location means.
(pm_params = colMeans(mod_csim))
yhat = rep(pm_params[1:5], each=30) #30 samples of observation yhats
resid = dat$chips - yhat
plot(resid)
plot(jitter(yhat), resid)

var(resid[yhat<7])
var(resid[yhat>11])

## location level residuals differ from the overall mean mu. 
lam_resid = pm_params[1:5] - pm_params["mu"]
plot(lam_resid)
abline(h=0, lty=2)
# we don't see see any obvious violations of our model assumptions.
summary(mod_sim)

#Posterior predictive simulation
# use these posterior samples to get Monte Carlo estimates the posterior predictive distribution.
#draws from the posterior distribution of mu and sigma 
#to simulate the posterior predictive distribution lam_pred of the mean for a new location.
(n_sim = nrow(mod_csim))## 15000
lam_pred = rgamma(n=n_sim, shape=mod_csim[,"mu"]^2/mod_csim[,"sig"]^2, 
                  rate=mod_csim[,"mu"]/mod_csim[,"sig"]^2)
hist(lam_pred)
#the posterior probability that the number of chips at a new location would be greater than 15. 
mean(lam_pred > 15)
#Using these lamda draws, we can get the observation level and simulate y_pred
#which takes into account the uncertainty in lambda
y_pred = rpois(n=n_sim, lambda=lam_pred)
hist(y_pred)

#what is the posterior probability that the next cookie produced in 
#Location 1 will have fewer than seven chips?
y_pred1 = rpois(n=n_sim, lambda=mod_csim[,"lam[1]"])
hist(y_pred1)
#posterior probability that the next cookie contains less than seven chips. 
mean(y_pred1 < 7)

