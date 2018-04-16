#two factors:wool and tension
rm(list=ls())
data("warpbreaks")
head(warpbreaks)
levels(warpbreaks$wool) # A,B
levels(warpbreaks$tension) #"L" "M" "H"
table(warpbreaks$wool, warpbreaks$tension)
#visualize the data with box plots
boxplot(breaks ~ wool + tension, data=warpbreaks)
boxplot(log(breaks) ~ wool + tension, data=warpbreaks)
#The different groups have more similar variance if we use the logarithm of breaks.
#It appears that there is a general decrease in breaks as we move from low to medium to high tension

###############one-way ANOVA model using tension only#########################
# y= mu + sig
# mu  ~ dnorm(0.0, 1.0/1.0e6)
# sig ~ dgamma(5/2.0, 5*2.0/2.0)

library("rjags")
mod1_string = " model {
    for( i in 1:length(y)) {
y[i] ~ dnorm(mu[tensGrp[i]], prec)
}

for (j in 1:3) {
mu[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*2.0/2.0)
sig = sqrt(1.0 / prec)
} "

set.seed(83)
str(warpbreaks)

data1_jags = list(y=log(warpbreaks$breaks), tensGrp=as.numeric(warpbreaks$tension))

params1 = c("mu", "sig")

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5e3)

## convergence diagnostics
plot(mod1_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)
summary(mod1_sim)

#The 95% posterior interval for the mean of mu[2] overlaps with both the mu[1] and mu[3]
# but the intervals for low and high group only slightly overlap.
# the means for low and high tension are different. 
# DIC model
dic1 = dic.samples(mod1, n.iter=1e3)

################### Two-way additive model(wool and tension have no interaction: no mu matrix) ################
#multiple factor ANOVA: With two factors, one with two levels and the other with three, we have six treatment groups,
# fit the additive model which treats the two factors separately with no interaction.no mu matrix
X = model.matrix( ~ wool + tension, data=warpbreaks)
head(X)
tail(X)
#By default, R has chosen the mean for wool A and low tension to be the intercept
mod2_string = " model {
    for( i in 1:length(y)) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = int + alpha*isWoolB[i] + beta[1]*isTensionM[i] + beta[2]*isTensionH[i]
}

int ~ dnorm(0.0, 1.0/1.0e6)
alpha ~ dnorm(0.0, 1.0/1.0e6)
for (j in 1:2) { 
beta[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(3/2.0, 3*1.0/2.0)
sig = sqrt(1.0 / prec)
} "

data2_jags = list(y=log(warpbreaks$breaks), isWoolB=X[,"woolB"], isTensionM=X[,"tensionM"], isTensionH=X[,"tensionH"])

params2 = c("int", "alpha", "beta", "sig")

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)

## convergene diagnostics
plot(mod2_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)

#summarize the results
summary(mod2_sim)
(dic2 = dic.samples(mod2, n.iter=1e3))
dic1
#This suggests there is much to be gained adding the wool factor to the model. 
# look again at the box plot with all six treatment groups.
boxplot(log(breaks) ~ wool + tension, data=warpbreaks)
lmod1 = lm(log(breaks)~ wool + tension, data=warpbreaks)
summary(lmod2)
#Our two-way model has a single effect for wool B and the estimate is negative.
#we would expect wool B to be associated with fewer breaks than its wool A counterpart on average
# This is true for low and high tension, but it's breaks are higher for wool B when there is medium tension.
# the effect for wool B is not consistent across tension levels, 
# so it may appropriate to add an interaction term. In  R, this would look like:
lmod2 = lm(log(breaks) ~ .^2, data=warpbreaks)
summary(lmod2)

#Adding the interaction, we get an effect for being in wool B and medium tension, 
#as well as for being in wool B and high tension. 
#There are now six parameters for the mean, 
#one for each treatment group, so this model is equivalent to the full cell means model. 

##############Two-way cell means model with interaction metrixs####################################
# mu will be a matrix with six entries : rows:woolGrp:2, columns: tensGrp: 3
mod3_string = " model {
    for( i in 1:length(y)) {
y[i] ~ dnorm(mu[woolGrp[i], tensGrp[i]], prec)
}

for (j in 1:max(woolGrp)) {
for (k in 1:max(tensGrp)) {
mu[j,k] ~ dnorm(0.0, 1.0/1.0e6)
}
}

prec ~ dgamma(3/2.0, 3*1.0/2.0)
sig = sqrt(1.0 / prec)
} "

str(warpbreaks)

data3_jags = list(y=log(warpbreaks$breaks), woolGrp=as.numeric(warpbreaks$wool), tensGrp=as.numeric(warpbreaks$tension))

params3 = c("mu", "sig")

mod3 = jags.model(textConnection(mod3_string), data=data3_jags, n.chains=3)
update(mod3, 1e3)

mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=5e3)
mod3_csim = as.mcmc(do.call(rbind, mod3_sim))

plot(mod3_sim, ask=TRUE)

## convergence diagnostics
gelman.diag(mod3_sim)
autocorr.diag(mod3_sim)
effectiveSize(mod3_sim)
raftery.diag(mod3_sim)
# compute the DIC and compare with our previous models.
(dic3 = dic.samples(mod3, n.iter=1e3))

#This suggests that the full model with interaction between wool and tension 
#(which is equivalent to the cell means model) is the best for explaining/predicting warp breaks.
summary(mod3_sim)
par(mfrow=c(3,2)) # arrange frame for plots for mu matrixs value distributions
densplot(mod3_csim[,1:6], xlim=c(2.0, 4.5))
#It might be tempting to look at comparisons between each combination of treatments
#Results are most reliable when we determine a relatively small number of hypotheses 
#we are interested in beforehand, collect the data, and statistically evaluate the evidence for them.

#go through our posterior samples and for each sample, find out which group has the smallest mean
#posterior probability that each of the treatment groups has the smallest mean
prop.table( table( apply(mod3_csim[,1:6], 1, which.min) ) ) # find rows with min Markov Chain Monte Carlo (MCMC) output
#The evidence supports wool B with high tension as the treatment that produces the fewest breaks.