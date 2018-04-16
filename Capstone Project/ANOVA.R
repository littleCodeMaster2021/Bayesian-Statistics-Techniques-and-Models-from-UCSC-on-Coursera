#As an example of a one-way ANOVA
data("PlantGrowth")
?PlantGrowth
head(PlantGrowth)
#group is a factor not continuous and visualize the data with box plots rather than scatter plots.
boxplot(weight ~ group, data=PlantGrowth) #It appears that treatment 2 has the highest mean yield

###########Modeling################################
#the reference analysis (with a noninformative prior) with a linear model in R
lmod = lm(weight ~ group, data=PlantGrowth)
summary(lmod)
plot(lmod) # for graphical residual analysis
anova(lmod)
#To recover the mean yield in treatment group 1, you would add the intercept term and the treatment 1 effect. 
model.matrix(lmod) # extract the X matrix.

# y= mu + sig
# mu  ~ dnorm(0.0, 1.0/1.0e6)
# sig ~ dgamma(5/2.0, 5*2.0/2.0)


# fit the cell means model: mu in JAGS
library("rjags")
mod_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dnorm(mu[grp[i]], prec[grp[i]]) #3 groups with different mu and variance
}

for (j in 1:3) {
mu[j] ~ dnorm(0.0, 1.0/1.0e6)
prec[j] ~ dgamma(5/2.0, 5*1.0/2.0)
sig[j] = sqrt( 1.0 / prec[j] )
}
} "

set.seed(82)
str(PlantGrowth)
data_jags = list(y=PlantGrowth$weight, 
                 grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
  inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(3,1.0,1.0))
  # set initialized value for each parameter in model, 
  #we can generate random sample: list("mu"= rnorm(...))
}

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains

#################Model checking: check for convergence of our MCMC################
plot(mod_sim)
gelman.diag(mod_sim)
autocorr.diag(mod_sim)
effectiveSize(mod_sim)
###########look at the residuals to see if there are any obvious problems with our model choice
(pm_params = colMeans(mod_csim))
yhat = pm_params[1:3][data_jags$grp]
resid = data_jags$y - yhat
plot(resid)
plot(yhat, resid)

#Again, it might be appropriate to have a separate variance for each group. 
#The posterior summary of the parameters.
summary(mod_sim) #check posterior for which group's mean was most affected by fitting separate group variances
HPDinterval(mod_csim)#calculates a 95% interval intervals of highest posterior density for each parameter#check if one of the treatments increases mean yield.
mean(mod_csim[,3] > mod_csim[,1]) #mu[3]>mu[1]
#There is a high posterior probability that the mean yield for treatment 2 is greater than the mean yield for the control group.
# If treatment 2 would be costly to put into production
#this treatment must increase mean yield by 10%. 
#What is the posterior probability that the increase?
mean(mod_csim[,3] > 1.1*mod_csim[,1])
#We have about 50/50 odds that adopting treatment 2 would increase mean yield by at least 10%.

HPDinterval(mod_csim[,3]-mod_csim[,1]) #calculate a 95% interval of highest posterior density (HPD) MU[3]-MU[1]
# ???1 in the model formula tells R to drop the intercept
#Because we used fairly noninformative priors for the ?? parameters in the analysis with JAGS,
mod_cm = lm(weight ~ -1 + group, data=PlantGrowth)
summary(mod_cm)