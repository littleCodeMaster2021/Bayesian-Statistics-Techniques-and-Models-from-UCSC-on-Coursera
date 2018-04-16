#An example of linear regression

library("car")
data("Leinhardt")
?Leinhardt
head(Leinhardt)
str(Leinhardt) #summary
summary(Leinhardt)
pairs(Leinhardt)
#simple linear regression model that relates infant mortality to per capita income.
plot(infant ~ income, data=Leinhardt)
#histogram plot
hist(Leinhardt$infant)
hist(Leinhardt$income)


#Since infant mortality and per capita income are positive and right-skewed quantities, 
#we consider modeling them on the logarithmic scale. 
#A linear model appears much more appropriate on this scale.

Leinhardt$loginfant = log(Leinhardt$infant)
Leinhardt$logincome = log(Leinhardt$income)
plot(loginfant ~ logincome, data=Leinhardt)
abline(lm(Leinhardt$loginfant~ Leinhardt$logincome), col='red', main='linear regression')

###################Modeling#############

#The reference Bayesian analysis (with a noninformative prior) is available directly in R
lmod = lm(loginfant ~ logincome, data=Leinhardt)
summary(lmod)
plot(lmod) # for graphical residual analysis

#fit this model in JAGS. A few countries have missing values, and we will omit those.
dat = na.omit(Leinhardt)


#correlation matrix
#Positive correlations are displayed in blue and negative correlations in red color.
#Two correlated variables will compete for the ability to predict the response variable, 
#leading to unstable estimates. 
#This is not a problem for prediction of the response, 
#if prediction is the end goal of the model. 
#But if our objective is to discover how the variables 
#relate to the response, we should avoid collinearity.
library(corrplot)
df <-data.frame(dat$loginfant,dat$logincome,dat$income,dat$infant)
M <- cor(df)
corrplot(M, type="upper", method="ellipse", tl.pos="d")
corrplot(M, type="lower", method="number", col="black", 
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")

library("rjags")
mod1_string = " model {
    for (i in 1:n) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = b[1] + b[2]*log_income[i] 
}

for (i in 1:2) {
b[i] ~ dnorm(0.0, 1.0/1.0e6)
}
## Initial guess of variance based on overall
## variance of education variable. Uses low prior
## effective sample size. Technically, this is not
## a true 'prior', but it is not very informative.

prec ~ dgamma(5/2.0, 5*10.0/2.0)
sig2 = 1.0 / prec
sig = sqrt(sig2)
} "

set.seed(72)
data1_jags = list(y=dat$loginfant, n=nrow(dat), 
                  log_income=dat$logincome)

params1 = c("b", "sig")

inits1 = function() {
  inits = list("b"=rnorm(2,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, inits=inits1, n.chains=3)
update(mod1, 1000) # burn-in

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5000)

mod1_csim = do.call(rbind, mod1_sim) # combine multiple chains


#MCMC convergence
#Before we check the inferences from the model, we should perform convergence diagnostics for our Markov chains.
plot(mod1_sim)

gelman.diag(mod1_sim)

autocorr.diag(mod1_sim)

autocorr.plot(mod1_sim)

effectiveSize(mod1_sim)

#We can get a posterior summary of the parameters in our model.
# results are for a regression model relating the logarithm of infant mortality to the logarithm of income.
summary(mod1_sim)

#these results are for a regression model: logarithm of infant mortality to the logarithm of income.
(pm_coef = colMeans(mod1_csim))

#The Raftery and Lewis diagnostic estimates how many iterations of
#the current chain would be required to reliably estimate the outer quantiles of the posterior.
raftery.diag(mod_csim)
#The dependence factor for many of the variables is large (>5.0), 
#indicating strong autocorrelation in the chains. 
#We would require a large number of iterations to reliably produce 95% probability
#intervals for the parameters.


##############Residual checks#################

#Checking residuals (observation- and the model's prediction for that value) 
#Residuals can reveal violations of the assumptions we made to specify the model. 
#Find any sign that the model is not linear, normally distributed, or that the observations are not independent (conditional on covariates).
#Fit the reference linear model to the un-transformed variables.
lmod0 = lm(infant ~ income, data=Leinhardt)
plot(resid(lmod0)) # to check independence (looks okay)
abline(h=0, col="blue") # predict and residual of linear regression
plot(predict(lmod0), resid(lmod0)) # to check for linearity, constant variance (looks bad)
abline(h=0, col="red")
qqnorm(resid(lmod0)) # to check Normality assumption (we want this to be a straight line)

#Fit to the log-transformed variables.  
#Check the residuals evaluated at the posterior mean of the parameters.
X = cbind(rep(1.0, data1_jags$n), data1_jags$log_income)
head(X)
(pm_params1 = colMeans(mod1_csim)) # posterior mean
yhat1 = drop(X %*% pm_params1[1:2])# delect b[1] and b[2]
resid1 = data1_jags$y - yhat1
plot(resid1) # against data index
abline(h=0, col="blue")
plot(yhat1, resid1) # against predicted values# residuals from the first model
abline(h=0, col="blue")
qqnorm(resid1) # checking normality of residuals
plot(predict(lmod), resid(lmod)) # to compare with reference linear model
abline(h=0, col="red")
sd(resid1) # standard deviation of residuals
rownames(dat)[order(resid1, decreasing=TRUE)[1:5]] # which countries have the largest positive residuals?

#The residuals look pretty good here (no patterns, shapes) except for two strong outliers, Saudi Arabia and Libya.
#When outliers appear, double check that they are not just errors in data entry. 
#If the values are correct,whether these data points really are representative of the data.
#If you conclude that they are not, you can justify dropping these data points from the data set.

#####################options for these outliers belong in the data set#################
#The first approach: find additional covariates that may be able to explain the outliers. 
# number of variables that about infant mortality above and beyond what income provides.
# two variables we haven't used yet: region and oil, perhaps this might explain part of the anomaly.
library("rjags")

mod2_string = " model {
for (i in 1:length(y)) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = b[1] + b[2]*log_income[i] + b[3]*is_oil[i] #add is_oil option
}

for (i in 1:3) {
b[i] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*10.0/2.0)
sig = sqrt( 1.0 / prec )
} "


set.seed(73)
data2_jags = list(y=dat$loginfant, log_income=dat$logincome,
                  is_oil=as.numeric(dat$oil=="yes"))
data2_jags$is_oil

params2 = c("b", "sig")

inits2 = function() {
    inits = list("b"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, inits=inits2, n.chains=3)
update(mod2, 1e3) # burn-in

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)

mod2_csim = as.mcmc(do.call(rbind, mod2_sim)) # combine multiple chains
#As usual, check the convergence diagnostics.
plot(mod2_sim)
gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)
#Get a posterior summary of the parameters in our model.
summary(mod2_sim)
#It looks like there is a positive relationship between oil-production and log-infant mortality. 
#Because these data are merely observational, we cannot say that oil-production causes an increase 
#in infant mortality (indeed that most certainly isn't the case)
#but we can say that they are positively correlated.
#Check the residuals:
X2 = cbind(rep(1.0, data1_jags$n), data2_jags$log_income, data2_jags$is_oil)
head(X2)
(pm_params2 = colMeans(mod2_csim)) # posterior mean
yhat2 = drop(X2 %*% pm_params2[1:3])
resid2 = data2_jags$y - yhat2
plot(resid2) # against data index
abline(h=0, col="red")
plot(yhat2, resid2) # against predicted values
abline(h=0, col="red")
mean(resid2)
sd(resid2) # standard deviation of residuals
#These look much better, although the residuals for Saudi Arabia and Libya are
#still more than three standard deviations away from the mean of the residuals.



##########################t another option when we are faced with strong outliers-t dis rather than normal
#models with the normal likelihood might be overly-influenced by outliers.
#T distribution is similar to the normal distribution, but it has thicker tails which can accommodate outliers.

#The t linear model: Notice that the t distribution has three parameters(positive "degrees of freedom" parameter)
#The smaller the degrees of freedom, the heavier the tails of the distribution. 
#We might fix the degrees of freedom to some number, or we can assign it a prior distribution.
mod3_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dt( mu[i], tau, df )
mu[i] = b[1] + b[2]*log_income[i] + b[3]*is_oil[i]
}

for (i in 1:3) {
b[i] ~ dnorm(0.0, 1.0/1.0e6)
}

df = nu + 2.0 # we want degrees of freedom > 2 to guarantee existence of mean and variance
nu ~ dexp(1.0)

tau ~ dgamma(5/2.0, 5*10.0/2.0) # tau is close to, but not equal to the precision
sig = sqrt( 1.0 / tau * df / (df - 2.0) ) # standard deviation of errors
} "

#the deviance information criterion (DIC): calculates the posterior mean of the log-likelihood and adds a penalty for model complexity.
#calculate the DIC for our first two models:
#the simple linear regression on log-income
dic.samples(mod1, n.iter=1e3)
#the second model where we add oil production.
dic.samples(mod2, n.iter=1e3)

#The better-fitting model has a lower DIC value.
# the gains we receive in deviance by adding the is_oil covariate outweigh the penalty 
#for adding an extra parameter. 

#The first number is the Monte Carlo estimated posterior mean deviance, which equals ???2 times the log-likelihood 
#???2 factor, a smaller deviance means a higher likelihood.
#a penalty for the complexity of our model:we can always 
#increase the likelihood of the model by making it more complex to fit the data exactly.
#This penalty is roughly equal to the effective number of parameters in your model. 
#With the first model, we had a variance parameter and two betas, for a total of three parameters.
#In the second model, we added one more beta for the oil effect.
#DIC is the last number reported with the title "Penalized deviance."
#find more details about the JAGS implementation through by entering ?dic.samples

###Another example###################
library("car")  # load the 'car' package
data("Anscombe")  # load the data set
?Anscombe  # read a description of the data
head(Anscombe)  # look at the first few lines of the data
pairs(Anscombe)  # scatter plots for each pair of variables

lmod = lm(education ~ income+young+urban, data=Anscombe)
summary(lmod)

library(corrplot)
df <-data.frame(Anscombe)
M <- cor(df)
corrplot(M, type="upper", method="ellipse", tl.pos="d")
corrplot(M, type="lower", method="number", col="black", 
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")




library("rjags")

mod3_string = " model {
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

data3_jags = as.list(Anscombe)

params3 = c("b", "sig")

inits3= function() {
  inits = list("b"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod3 = jags.model(textConnection(mod3_string), data=data3_jags, inits=inits3, n.chains=3)
update(mod3, 1e3) # burn-in

mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=1e5) # run many times

mod3_csim = as.mcmc(do.call(rbind, mod3_sim)) # combine multiple chains MCMC
#As usual, check the convergence diagnostics.
plot(mod3_sim)
gelman.diag(mod3_sim)
autocorr.diag(mod3_sim)
autocorr.plot(mod3_sim)
effectiveSize(mod3_sim)
#Get a posterior summary of the parameters in our model.
summary(mod3_sim)
mean(mod3_csim[,3] <mod_csim[,1])#posterior probability that b[3]<b[1]
dic.samples(mod3, n.iter=1e2)
#DIC is the last number reported with the title "Penalized deviance."