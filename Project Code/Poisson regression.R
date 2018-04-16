# Check current working directory
getwd()

# set directory
setwd("C:/Users/shuai/Desktop")

#Poisson regression
library("COUNT")
#badh : binary variable
data("badhealth")
?badhealth
head(badhealth)
#Check na
any(is.na(badhealth))
#visualize these data
hist(badhealth$numvisit, breaks=20)
plot(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==0, xlab="age", ylab="log(visits)")
points(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==1, col="red")
#both age and bad health are related to the number of doctor visits
#If we believe the age/visits relationship is different between healthy and non-healthy populations
#we should also include an interaction term.
library("rjags")

####################fit the full model#######################
mod_string = " model {
    for (i in 1:length(numvisit)) {
numvisit[i] ~ dpois(lam[i])
log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
}
#non-informative prior
int ~ dnorm(0.0, 1.0/1e6)
b_badh ~ dnorm(0.0, 1.0/1e4)
b_age ~ dnorm(0.0, 1.0/1e4)
b_intx ~ dnorm(0.0, 1.0/1e4)
} "

set.seed(102)

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

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
dic = dic.samples(mod, n.iter=1e3) #get DIC and effective parameters

##################Model checking##############################
#delete the first column:numvisit
X = as.matrix(badhealth[,-1])
X = cbind(X, with(badhealth, badh*age))
head(X)
#pm_params2 = colMeans(mod_csim)
(pmed_coef = apply(mod_csim, 2, median)) # take median for each column
#apply the inverse of the link function to get predictions for lambda
llam_hat = pmed_coef["int"] + X %*% pmed_coef[c("b_badh", "b_age", "b_intx")]
lam_hat = exp(llam_hat)
hist(lam_hat)
resid = badhealth$numvisit - lam_hat
plot(resid) # the data were ordered
plot(lam_hat, badhealth$numvisit)
# add zero and one line
abline(0.0, 1.0)
plot(lam_hat[which(badhealth$badh==0)], resid[which(badhealth$badh==0)], xlim=c(0, 8), ylab="residuals", xlab=expression(hat(lambda)), ylim=range(resid))
points(lam_hat[which(badhealth$badh==1)], resid[which(badhealth$badh==1)], col="red")
#the variability increases for values predicted at higher values 
#since the mean is also the variance in the Poisson distribution
#observations predicted to have about two visits should have variance about two not about 20
#observations predicted to have about six visits should have variance about six not about 30

#calculate variance

var(resid[which(badhealth$badh==0)])
var(resid[which(badhealth$badh==1)])

#The fact that we observed so much more variability than we expected indicates that
#either the model fits poorly. Meaning that the covariance don't
#explain enough of the variability in the data. 
#Or the data are over dispersed for this plus on likely that we've chosen. 

#the model fits poorly (meaning the covariates don’t explain enough of the variability in the data)
#the data are “overdispersed” for the Poisson likelihood we have chosen.
#This is a common issue with count data. 
#If the data are more variable than the Poisson likelihood would suggest,
#a good alternative is the negative binomial distribution
summary(mod_sim)

#For healthy individuals, it appears that age has a positive association with number of doctor visits.
#bad health is associated with an increase in expected number of visits. T
#Hence, for people with bad health, age is essentially unassociated with number of visits.
#What is the posterior probability that the individual with poor health will have more doctor visits? 
#This goes beyond the posterior probabilities we've calculated comparing expected responses
#in the previous lessons. Here, we need to create Monte Carlo 
#samples for the actual responses themselves. 
#We should take the Monte Carlo samples of the model parameters, and for each of those, 
#drawing a sample from the likelihood. 
#the healthy one is Person 1 and the unhealthy one is Person 2
x1 = c(0, 35, 0) # good health
x2 = c(1, 35, 35) # bad health
#beta: posterior samples of the model parameters are stored in mod_csim:
head(mod_csim) #beta parameter values
#compute the linear part of the predictor:
loglam1 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x1
loglam2 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x2
#apply the inverse link: get lambda:
lam1 = exp(loglam1)
lam2 = exp(loglam2)
# use these samples for the lambda for each individual and simulate y using the likelihood:
(n_sim = length(lam1))
y1 = rpois(n=n_sim, lambda=lam1)
y2 = rpois(n=n_sim, lambda=lam2)
plot(table(factor(y1, levels=0:18))/n_sim, ylab="posterior prob.", xlab="visits")
points(table(y2+0.1)/n_sim, col="red")
# What is the probability that the person with poor health will have more doctor visits than the person with good health?
mean(y2 > y1)
#Because we used our posterior samples for the model parameters in our simulation
#this posterior predictive distribution on the number of visits for these two new individuals
#takes into account our uncertainty in the model estimates. 
#This is a more honest/realistic distribution than we would get if we had fixed the model parameters 
#at their MLE or posterior means and simulated data for the new individuals.

#If t is the amount of time that we observe, and lamda is the rate of events per unit of time,
#then the expected number of events is tλ and the distribution of
#the number of events in this time interval is Poisson(tλ).

# Q7
#fit the model using N(0,102) priors for the intercept and both coefficients. 
#look at the residuals (don't forget to multiply lam_hat by days_active to
#obtain the model's predicted mean number of calls).
#One major advantage of the Poisson linear model is that the 
#log transformation appears in the link function. 
#We are taking the log of the mean rather than the response itself. 

mod3_string = " model {
for (i in 1:length(calls)) {
calls[i] ~ dpois(lam[i]*days_active[i])
log(lam[i]) = b[0] + b[1]*age[i] + b[2]*isgroup2[i]
}

for (j in 1:3) {
b[j] ~ dnorm(0.0, 1.0/1e2)
}
} "

data3_jags = as.list(dat)

params3 = c("b0", "b")


mod3 = jags.model(textConnection(mod3_string), data=data3_jags, n.chains=3)
update(mod3, 1e3)

mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=15e3)
mod3_csim = as.mcmc(do.call(rbind, mod3_sim))

## convergence diagnostics
plot(mod3_sim)

gelman.diag(mod3_sim)
autocorr.diag(mod3_sim)
autocorr.plot(mod3_sim)
effectiveSize(mod3_sim)

## compute DIC
dic3 = dic.samples(mod3, n.iter=1e3)

pmod3_coef = apply(mod3_csim, 2, mean)

X = as.matrix(dat[,-c(1,2)])
X = X[,c(2,1)]

#llam_hat3 = pmod3_coef['b0'] + X %*% pmod3_coef[-3] #select ['b0'] and delete the third column 
llam_hat3 =  X %*% pmod3_coef 
lam_hat3 = llam_hat3*dat$days_active
#posterior probability that beta the coefficient is greater than 0?
# beta parameter value is in second column of mod_csim:
mean(mod3_csim[,2] > 0)
