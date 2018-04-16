#Random intercept linear model
#We’ll do this with a hierarical model, where each region has its own intercept.
library("car")
data("Leinhardt")
?Leinhardt
str(Leinhardt)
pairs(Leinhardt)
head(Leinhardt)
#worked with infant mortality and income on the logarithmic scale. 
#Recall also that we had to remove some missing data.
dat = na.omit(Leinhardt)
dat$logincome = log(dat$income)
dat$loginfant = log(dat$infant)
str(dat)
# fit the proposed model:
library("rjags")
mod_string = " model {
for (i in 1:length(y)) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = a[region[i]] + b[1]*log_income[i] + b[2]*is_oil[i]
}

for (j in 1:max(region)) {
a[j] ~ dnorm(a0, prec_a)
}

a0 ~ dnorm(0.0, 1.0/1.0e6)
prec_a ~ dgamma(1/2.0, 1*10.0/2.0)
tau = sqrt( 1.0 / prec_a )

for (j in 1:2) {
b[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*10.0/2.0)
sig = sqrt( 1.0 / prec )
} "

set.seed(116)
data_jags = list(y=dat$loginfant, log_income=dat$logincome,
                 is_oil=as.numeric(dat$oil=="yes"), region=as.numeric(dat$region))
#data_jags$is_oil
#table(data_jags$is_oil, data_jags$region)
# show parameters
params = c("a0", "a", "b", "sig", "tau")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3) # burn-in

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)

mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combine multiple chains

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)
dic = dic.samples(mod, n.iter=1e3)
#It appears that this model is an improvement over the non-hierarchical one we fit earlier. 
#Notice that the penalty term, which can be interpreted as the “effective” number of parameters, 
#is less than the actual number of parameters (nine). There are fewer “effective” parameters
#because they are “sharing” information or “borrowing strength” from each other in the hierarhical 
#structure. If we had skipped the hierarchy and fit one intercept, there would 
#have been four parameters. If we had fit separate, independent intercepts for each region, 
#there would have been seven parameters (which is close to what we ended up with).

#Check the posterior summary:
summary(mod_sim)

#the intercepts a[1],a[2],a[3],a[4],a0 do not have a real interpretation 
#because they correspond to the mean response for a country
#that does not produce oil and has $0 log-income per capita 
#a0 is the overall mean of intercepts and 
#tau as the standard deviation of intercepts across regions.

###########################Other models: #############################
#We have not investigated adding interaction terms, which might be 
#appropriate. We only considered adding hierarchy on the intercepts, 
#but in reality nothing prevents us from doing the same for other terms
#in the model, such as the coefficients for income and oil. 
#We could try any or all of these alternatives and see how 
#the DIC changes for those models. This, together with other model checking techniques 
#we have discussed could be used to identify your best model that you can use to make 
#inferences and predictions.

# Q4
#a model that assumes no hierarchy (the ANOVA cell means model). 
#We can approximate the posterior estimates for the five industries: group 1 to 5 
#means (theta) under a noninformative prior by simply calculating the sample mean growth 
#for the five industries. 
library(rjags)
mod_string = " model {
for (i in 1:length(y)) {
y[i] ~ dnorm(theta[grp[i]], prec)
}

for (j in 1:max(grp)) {
#posterior estimates 
theta[j] ~ dnorm(mu, tau_sq)
}

#noninformative prior
mu ~ dnorm(0, 1/1e6)
tau_sq ~ dgamma(1.0/2.0, 1.0*3.0/2.0)
prec ~ dgamma(2.0/2.0, 2*1/2)
sig = sqrt(1/prec)

} "

set.seed(113)

data_jags = as.list(dat)

params = c("theta", "mu", "sig")

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

## compute DIC: get DIC and effective parameters
dic = dic.samples(mod, n.iter=1e3)

pm_params = apply(mod_csim, 2, mean)
means_theta = pm_params[-c(1,2)]

means_anova = tapply(dat$y, INDEX=dat$grp, FUN=mean)
## dat is the data read from pctgrowth.csv

plot(means_anova)
points(means_theta, col="red") ## where means_theta are the posterior point estimates for the industry means.
