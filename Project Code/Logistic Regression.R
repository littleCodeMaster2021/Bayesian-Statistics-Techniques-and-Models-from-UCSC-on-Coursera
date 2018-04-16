#For an example of logistic regression
#The response variable is r, which takes on values of 0 or 1. 
#We will remove some rows from the data set which contain missing values.
library("boot")
data("urine")
?urine
head(urine)
#remove some rows from the data set which contain missing values.
dat = na.omit(urine)
#pairwise scatter plots of the seven variables
pairs(dat)

#Two correlated variables will compete for the ability to predict the response variable, 
#leading to unstable estimates. This is not a problem for prediction of the response,
#if prediction is the end goal of the model. 
#But if our objective is to discover how the variables relate to the response,
#we should avoid collinearity.

#estimate the correlation among these variables 
library("corrplot")
Cor = cor(dat)
corrplot(Cor, type="upper", method="ellipse", tl.pos="d")
corrplot(Cor, type="lower", method="number", col="black", 
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")

############## Variable selection###################
# find out which variables are related to the presence of goal variable:
# Method 1:  fit several models that include different sets of variables and see which one has the best DIC.
# Method 2: use a linear model where the priors for the beta coefficients favor values near 0 (indicating a weak relationship).
# This way, the burden of establishing association lies with the data.
# If there is not a strong signal, we assume it doesn't exist.
# Rather than tailoring a prior for each individual beta based on the scale its covariate takes values on
#it subtracts the mean and divide by the standard deviation for each variable.
X = scale(dat[,-1], center=TRUE, scale=TRUE)
head(X[,"gravity"])
# priors for the beta coefficients favor values near 0 (indicating a weak relationship). 
colMeans(X) # column mean
apply(X, 2, sd) # 2 is for each column sd: standard deviation

###########Orginial Full Modeling ##########################################

#Our prior for the beta (b in the model) coefficients will be the double exponential distribution,
ddexp = function(x, mu, tau) {
  0.5*tau*exp(-tau*abs(x-mu)) 
}
# double exponential distribution
curve(ddexp(x, mu=0.0, tau=1.0), from=-5.0, to=5.0, ylab="density", main="Double exponential\ndistribution") 
# normal distribution
curve(dnorm(x, mean=0.0, sd=1.0), from=-5.0, to=5.0, lty=2, add=TRUE) 
legend("topright", legend=c("double exponential", "normal"), lty=c(1,2), bty="n")

library("rjags")

mod1_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dbern(p[i])
logit(p[i]) = int + b[1]*gravity[i] + b[2]*ph[i] + b[3]*osmo[i] + b[4]*cond[i] + b[5]*urea[i] + b[6]*calc[i]
}
int ~ dnorm(0.0, 1.0/25.0)
for (j in 1:6) {
b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
}
} "

set.seed(92)
head(X)

data_jags = list(y=dat$r, gravity=X[,"gravity"], ph=X[,"ph"], osmo=X[,"osmo"], cond=X[,"cond"], urea=X[,"urea"], calc=X[,"calc"])

params = c("int", "b")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

## convergence diagnostics
plot(mod1_sim, ask=TRUE)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)
# look the result
summary(mod1_sim)
par(mfrow=c(3,2))
densplot(mod1_csim[,1:6], xlim=c(-3.0, 3.0))
colnames(X) # variable names
#It is clear that the coefficients for variables b[1], b[4], and b[6]'s means are not 0. 
# The posterior distribution for b[2] and b[3] looks like the prior, and is almost centered on 0 still
# b[3] and b[2] are not strong predictor of y
#b[5] is a borderline case. However, if we refer back to our correlations among the variables
# we see that b[5] is highly correlated with b[1]0.83, so we opt to remove it.

################Our second model with removal variables######################
library("rjags")
mod2_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dbern(p[i])#The response binary vector is r: 0 or 1 with likelihood ~ Bernoulli Distribution
logit(p[i]) = int + b[1]*gravity[i] + b[2]*cond[i] + b[3]*calc[i]
}
int ~ dnorm(0.0, 1.0/25.0)
for (j in 1:3) {
b[j] ~ dnorm(0.0, 1.0/25.0) # noninformative for logistic regression:prior for betas
}
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim, ask=TRUE)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)


## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)
dic2 = dic.samples(mod2, n.iter=1e3)
raftery.diag(mod2_csim)
# view results
summary(mod2_sim)
HPDinterval(mod2_csim)
par(mfrow=c(3,1))
densplot(mod2_csim[,1:3], xlim=c(-3.0, 3.0))
HPDinterval(mod2_csim[,3]-mod2_csim[,1]) #calculate a 95% interval of highest posterior density (HPD) b[3]-b[1]
colnames(X)[c(1,4,6)] # variable names
#The DIC is actually better for the first model. 
#Note that we did change the prior between models, 
#and we should not use the DIC to choose between priors. 
#Hence comparing DIC between these two models may not be a fair comparison. 
#Nevertheless, they both yield essentially the same conclusions.
#Higher values of b[1] and b[3] are associated with higher probabilities of r
#while higher values of b[4] are associated with lower probabilities of r.
# Ways to improve this model: transformations of variables, different priors, and interactions between the predictors...

##############Prediction from a logisic regression model: convert linear model to logistic model ######
#logit of p as a linear model:E(y)=p
#use the posterior means as estimates of the parameters.
(pm_coef = colMeans(mod2_csim))
# prediction for r when beta parameters are at their means: 
r1 <- 1/(1+exp(1)^(-pm_coef[4]))
#make a prediction for a new r whose b[1] is average, b[4] is mean-1sd, and b[6] is mean+1sd
r2<- 1/(1+exp(1)^-(pm_coef[4]+pm_coef[1]*0+pm_coef[2]*(-1)+pm_coef[3]*(1)))

#If we want to make predictions in terms of the original x variables without scale:
# Method1: For each xx variable, subtract the mean and divide by the standard deviation for that variable in the original data set used to fit the model
# Method2: Re-fit the model without centering and scaling the covariates

##################Predictive checks:calculate residuals ##################
# take the X matrix and matrix multiply it with the posterior means of the coefficients:
pm_Xb = pm_coef["int"] + X[,c(1,4,6)] %*% pm_coef[1:3] # take the 1,4,6 column from X
rows<-matrix(1:nrow(pm_Xb), nrow=nrow(pm_Xb))
pm_Xb <- cbind(rows,pm_Xb)
#These phat values are the model's predicted posterior probability of r for each data point.
#pass these linear values through the inverse of the link function 
phat = 1.0 / (1.0 + exp(-pm_Xb[,2]))
head(phat)

#plot these predicted values against the actual outcome.
plot(rows,phat,col="red",
  main = "prediction and actual points",
  xlab = "# of data: red is prediction value, blue is actual value",
 ylab = "prediction and actual values")
par(new=TRUE)
plot(rows,jitter(dat$r),col="blue",
     main = "prediction and actual points",
     xlab = "# of data: red is prediction value, blue is actual value",
     ylab = "prediction and actual values")


# If the model tells us the probability is higher than 0.5
#we will classify the observation as a 1 
#if it is less than 0.5, we will classify it as a 0.
#tabulate these classifications against the truth to see how well the model predicts the original data.
(tab0.5 = table(phat > 0.5, data_jags$y))
# The correct classification rate:
sum(diag(tab0.5)) / sum(tab0.5)

# choose to lower our threshold for classifying data points as 1s.change it to 0.3. 
#if the model says the probability is greater than 0.3, 
# we will classify it as 1
(tab0.3 = table(phat > 0.3, data_jags$y))
sum(diag(tab0.3)) / sum(tab0.3)

#appears that the accurate cases (high probability of correct responses) are well captured by the model.
#In this exercise, we obtained a point estimate of the coefficients and used that to obtain a 
#point estimate of the probabilities. If we want posterior distributions for the probabilities, 
#we could apply the inverse link transformation to each iteration of the coefficients.


# we did indeed increase our chances of detecting a true positive.
#We could repeat this exercise for many thresholds between 0 and 1,
#and each time calculate our error rates. 
#This is equivalent to calculating what is called the ROC (receiver-operating characteristic) curve,
#which is often used to evaluate classification techniques.

# We could get a less biased assessment of how well our model performs if we calculated these tables for data 
#that were not used to fit the model. For example, before fitting the model, 
#you could withhold a set of randomly selected "test" data points, and use the model fit 
#to the rest of the "training" data to make predictions on your "test" set.

  