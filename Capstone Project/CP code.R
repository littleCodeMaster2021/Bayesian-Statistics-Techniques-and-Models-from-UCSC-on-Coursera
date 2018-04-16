rm(list=ls())
df<- read.csv("Boston_house_price_data.csv",header=TRUE)
#check na
any(is.na(badhealth))
#remove some rows from the data set which contain missing values.
df = na.omit(df)
head(df)
summary(df)

#estimate the correlation among these variables 
library("corrplot")
Cor = cor(df)
corrplot(Cor, type="upper", method="ellipse", tl.pos="d")
corrplot(Cor, type="lower", method="number", col="black", 
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")

# Method 2: use a linear model where the priors for the beta coefficients favor values near 0 (indicating a weak relationship).
# This way, the burden of establishing association lies with the data.
# If there is not a strong signal, we assume it doesn't exist.
# Rather than tailoring a prior for each individual beta based on the scale its covariate takes values on
#it subtracts the mean and divide by the standard deviation for each variable.
#Should not scale category variable
X = scale(df[,c(1:3,5:8,10:14)], center=TRUE, scale=TRUE)
# priors for the beta coefficients favor values near 0 (indicating a weak relationship). 
colMeans(X) # column mean
apply(X, 2, sd) # 2 is for each column sd: standard deviation


###########Orginial Full Modeling ##########################################

#Our prior for the beta (b in the model) coefficients will be the double exponential distribution
#the Laplace distribution has fatter tails than the normal distribution.
#Remove rad which as really strong correlation with tax
library("rjags")
mod1_string = " model {
for (i in 1:length(y)) {
#likelihood of the data
y[i] ~ dnorm(mu[i], prec)
mu[i] = int + b[1]*crim[i] + b[2]*zn[i] + b[3]*indus[i] + b[4]*nox[i] + b[5]*rm[i] + b[6]*age[i]+ b[7]*dis[i]
+ b[8]*chas[i]+ b[9]*tax[i]+ b[10]*ptratio[i]+b[11]*black[i]+ b[12]*lstat[i]}

int ~ dnorm(0.0, 1.0/1.0e6) #non-informative prior with large variance
for (j in 1:12) {
b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
}

prec ~ dgamma(3/2.0, 3*10.0/2.0)
sig2 = 1.0 / prec
#gave a prior to the precision: a prior for sigma squared
sig = sqrt(sig2)
} "


set.seed(66)
head(X)

data_jags = list(y=X[,"medv"], crim=X[,"crim"], zn=X[,"zn"], indus=X[,"indus"], nox=X[,"nox"], rm=X[,"rm"], age=X[,"age"],
                 dis=X[,"dis"],chas=df[,"chas"],tax=X[,"tax"],ptratio=X[,"ptratio"],black=X[,"black"],lstat=X[,"lstat"])

params = c("int", "b","sig")

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
par(mfrow=c(4,3)) # draw 6 plots for b[1]...b[12]
densplot(mod1_csim[,1:12], xlim=c(-1.0, 1.0))
#traceplot(mod1_csim[,1:12], xlim=c(2000, 7000))
colnames(X) # variable names
#It is clear that the coefficients for variables b[1], b[4], and b[6]'s means are not 0. 
# The posterior distribution for b[2] and b[3] looks like the prior, and is almost centered on 0 still
# b[3] and b[2] are not strong predictor of y
#b[5] is a borderline case. However, if we refer back to our correlations among the variables
# we see that b[5] is highly correlated with b[1]0.83, so we opt to remove it.

# histogram:
h <-hist(df$zn, 
     main="Histogram for proportion of residential land", 
     xlab="zn", 
     border="blue", 
     col="bisque",
     xlim=c(0,100),
     ylim=c(0,400),
     las=1, 
     breaks=8)
xfit<-seq(min(df$zn),max(df$zn),length=10) 
yfit<-dnorm(xfit,mean=mean(df$zn),sd=sd(df$zn)) 
yfit <- yfit*diff(h$mids[1:2])*length(df$zn) 
lines(xfit, yfit, col="red", lwd=2)

#############Model without log and squares######################
library("rjags")
mod2_string = " model {
for (i in 1:n) {
#likelihood of the data
y[i] ~ dnorm(mu[i], prec)
#add the linear model: mu[i] is linear
mu[i] = b[1]+ b[2]*nox[i] + b[3]*rm[i] + b[4]*dis[i]
+ b[5]*chas[i]+ b[6]*ptratio[i]+ b[7]*lstat[i]  
}

for (i in 1:7) {
b[i] ~ dnorm(0.0, 1.0/1.0e6) #non-informative prior with large variance
}

prec ~ dgamma(3/2.0, 3*10.0/2.0)
sig2 = 1.0 / prec
#gave a prior to the precision: a prior for sigma squared
sig = sqrt(sig2)
} "

set.seed(72)
data2_jags = list(y=df$medv, n=nrow(df), 
                  nox=df$nox,rm=df$rm,dis=df$dis,chas=df$chas,ptratio=df$ptratio,
                  lstat=df$lstat)

params2 = c("b", "sig") #b[1] and b[2] and standard deviation: sig

inits2 = function() { #b[1] and b[2]
  inits = list("b"=rnorm(7,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

#In each chain: it'll initialize the chain with these random draws. 
#So we'll have different starting values for each chain. 

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, inits=inits2, n.chains=3)
#ran the model for 1,000 iterations but it didn't keep the samples. 
update(mod2, 1000) # burn-in

#mod2_sim: posterior simulation that we're going to keep for model inference
mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5000)

# combined simulation:stacking the matrices that contain the simulations themselves.
#Inside mod1.sim the simulations are stored as matrices. 
mod2_csim = do.call(rbind, mod2_sim) # combine multiple chains


#MCMC convergence
#Before we check the inferences from the model, we should perform convergence diagnostics for our Markov chains.
plot(mod2_sim)

gelman.diag(mod2_sim)

autocorr.diag(mod2_sim)

autocorr.plot(mod2_sim)

effectiveSize(mod2_sim)

#We can get a posterior summary of the parameters in our model.
# results are for a regression model relating the logarithm of infant mortality to the logarithm of income.
summary(mod2_sim)
dic2 = dic.samples(mod2, n.iter=1e3)
#These effective sample sizes are still okay, 
#if all we're interested in is the posterior mean for these two parameters. 
#But if we want to create something like a 95% posterior probability interval
#we want to have more effective samples and run the chains for longer. 

#these results are for a regression model: logarithm of infant mortality to the logarithm of income.
(pm_coef = colMeans(mod2_csim))


#############Model with log and squares######################
df1<-df[,c(4,6,13)]
df1$logmedv = log(df$medv)
df1$lognox = log(df$nox)
df1$logdis = log(df$dis)
df1$sqrptratio = (df$ptratio)*(df$ptratio)

library("rjags")
mod2_string = " model {
for (i in 1:n) {
#likelihood of the data
y[i] ~ dnorm(mu[i], prec)
#add the linear model: mu[i] is linear
mu[i] = b[1]+ b[2]*nox[i] + b[3]*rm[i] + b[4]*dis[i]
+ b[5]*chas[i]+ b[6]*ptratio[i]+ b[7]*lstat[i]
}

# prior of the coefficients
for (i in 1:7) {
b[i] ~ dnorm(0.0, 1.0/1.0e6) #non-informative prior with large variance
}

prec ~ dgamma(3/2.0, 3*10.0/2.0)
sig2 = 1.0 / prec
#gave a prior to the precision: a prior for sigma squared
sig = sqrt(sig2)
} "

set.seed(50)
data2_jags = list(y=df1$logmedv, n=nrow(df1), 
                  nox=df1$lognox,rm=df1$rm,dis=df1$logdis,chas=df1$chas,ptratio=df1$sqrptratio,
                  lstat=df1$lstat)

params2 = c("b", "sig") #b[1] and b[2] and standard deviation: sig

inits2 = function() { #b[1] and b[2]
  inits = list("b"=rnorm(7,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

#In each chain: it'll initialize the chain with these random draws. 
#So we'll have different starting values for each chain. 

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, inits=inits2, n.chains=3)
#ran the model for 1,000 iterations but it didn't keep the samples. 
update(mod2, 1000) # burn-in

#mod2_sim: posterior simulation that we're going to keep for model inference
mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5000)

# combined simulation:stacking the matrices that contain the simulations themselves.
#Inside mod1.sim the simulations are stored as matrices. 
mod2_csim = do.call(rbind, mod2_sim) # combine multiple chains



#MCMC convergence
#Before we check the inferences from the model, we should perform convergence diagnostics for our Markov chains.
plot(mod2_sim)

gelman.diag(mod2_sim)

autocorr.diag(mod2_sim)

autocorr.plot(mod2_sim)

effectiveSize(mod2_sim)

#We can get a posterior summary of the parameters in our model.
# results are for a regression model relating the logarithm of infant mortality to the logarithm of income.
summary(mod2_sim)

par(mfrow=c(4,2)) # draw 6 plots for b[1]...b[12]
densplot(mod2_csim[,1:7], xlim=c(-5, 5))
traceplot(mod2_csim[,1:7], xlim=c(2000, 6000))

dic2 = dic.samples(mod2, n.iter=1e3)
#These effective sample sizes are still okay, 
#if all we're interested in is the posterior mean for these two parameters. 
#But if we want to create something like a 95% posterior probability interval
#we want to have more effective samples and run the chains for longer. 

#these results are for a regression model: logarithm of infant mortality to the logarithm of income.
(pm_coef = colMeans(mod2_csim))

raftery.diag(mod2_csim)

#Check the residuals:
X2 = cbind(rep(1.0, data2_jags$n), data2_jags$nox, data2_jags$rm,data2_jags$dis,data2_jags$chas,data2_jags$ptratio,data2_jags$lstat)
head(X2)
(pm_params2 = colMeans(mod2_csim)) # posterior mean
yhat2 = drop(X2 %*% pm_params2[1:7])
resid2 = data2_jags$y - yhat2

plot(resid2,main="Residual Plot") # against data index
abline(h=0, col="red")


plot(yhat2, resid2,main="Residual Plot") # against predicted values
abline(h=0, col="red")
mean(resid2)
sd(resid2) # standard deviation of residuals
mean(abs(resid2)>mean(resid2)+2*sd(resid2))
qqnorm(resid2) # to check Normality assumption (we want this to be a straight line)

df1[rownames(df1)[order(abs(resid2), decreasing=TRUE)[1:5]],] # largest absolute residual value


d<-density(df1$logmedv)#actual
d1 <-density(yhat2) #predictive
plot(d,xlab="log transformation for House Price's median value",col="red",main="Density plot for actual and predictive values")
par(new=TRUE)
plot(d1,xlab="log transformation for House Price's median value",col="blue",main="Density plot for actual and predictive values")
labels = c("Actual values", "Predictive values")
legend("topright",legend = labels, col=c("red", "blue"), lty=c(1,1), cex=0.8)
#These look much better, although the residuals for Saudi Arabia and Libya are
#still more than three standard deviations away from the mean of the residuals.


##########plot for actual and predictive values###############
rows <-c(1:nrow(yhat2))
plot(rows,yhat2,col="red",
     main = "prediction and actual values comparsion's plot",
     xlab = "# of data: red is prediction value, blue is actual value",
     ylab = "prediction and actual values")
par(new=TRUE)
plot(rows,jitter(data2_jags$y),col="blue",
     main = "prediction and actual values comparsion's plot",
     xlab = "# of data: red is prediction value, blue is actual value",
     ylab = "prediction and actual values")


##########################t another option when we are faced with strong outliers-t dis rather than normal
#models with the normal likelihood might be overly-influenced by outliers.
#T distribution is similar to the normal distribution, but it has thicker tails which can accommodate outliers.
#The t linear model: Notice that the t distribution has three parameters(positive "degrees of freedom" parameter)
#The smaller the degrees of freedom, the heavier the tails of the distribution. 
#We might fix the degrees of freedom to some number, or we can assign it a prior distribution.

mod3_string = " model {
for (i in 1:length(y)) {
y[i] ~ dt( mu[i], tau, df )
mu[i] = b[1]+ b[2]*nox[i] + b[3]*rm[i] + b[4]*dis[i]
+ b[5]*chas[i]+ b[6]*ptratio[i]+ b[7]*lstat[i]
}

for (i in 1:7) {
b[i] ~ dnorm(0.0, 1.0/1.0e6)
}

df = nu + 2.0 # we want degrees of freedom > 2 to guarantee existence of mean and variance
nu ~ dexp(1.0)

tau ~ dgamma(3/2.0, 3*10.0/2.0) # tau is close to, but not equal to the precision
sig = sqrt( 1.0 / tau * df / (df - 2.0) ) # standard deviation of errors
} "

set.seed(50)
data3_jags = list(y=df1$logmedv, 
                  nox=df1$lognox,rm=df1$rm,dis=df1$logdis,chas=df1$chas,ptratio=df1$sqrptratio,
                  lstat=df1$lstat)

params3 = c("b", "nu", "sig") #b[1] and b[2] and standard deviation: sig

inits3 = function() { #b[1] and b[2]
  inits = list("b"=rnorm(n=7,0.0,100.0), "nu" = rexp(n=1, rate = 1), "tau"=rgamma(n=1,1.0,1.0))
}

#In each chain: it'll initialize the chain with these random draws. 
#So we'll have different starting values for each chain. 

mod3 = jags.model(textConnection(mod3_string), data=data3_jags, inits=inits3, n.chains=3)
#ran the model for 1,000 iterations but it didn't keep the samples. 
update(mod3, 1000) # burn-in

#mod2_sim: posterior simulation that we're going to keep for model inference
mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=5000)

# combined simulation:stacking the matrices that contain the simulations themselves.
#Inside mod1.sim the simulations are stored as matrices. 
mod3_csim = do.call(rbind, mod3_sim) # combine multiple chains


#MCMC convergence
#Before we check the inferences from the model, we should perform convergence diagnostics for our Markov chains.
plot(mod3_sim)

gelman.diag(mod3_sim)

autocorr.diag(mod3_sim)

autocorr.plot(mod3_sim)

effectiveSize(mod3_sim)

#We can get a posterior summary of the parameters in our model.
# results are for a regression model relating the logarithm of infant mortality to the logarithm of income.
summary(mod3_sim)




#the deviance information criterion (DIC): calculates the posterior mean of the log-likelihood and adds a penalty for model complexity.
#calculate the DIC for our first two models:
#the simple linear regression on log-income
dic3 <- dic.samples(mod3, n.iter=1e3)
pm_coef3 = colMeans(mod3_csim)

#Check the residuals:
X3 = cbind(rep(1.0, length(data3_jags$y)), data3_jags$nox, data3_jags$rm,data3_jags$dis,data3_jags$chas,data3_jags$ptratio,data3_jags$lstat)
head(X3)
(pm_params3 = colMeans(mod3_csim)) # posterior mean
yhat3 = drop(X3 %*% pm_params3[1:7])
resid3 = data3_jags$y - yhat3

plot(resid3,main="Residual Plot") # against data index
abline(h=0, col="red")


plot(yhat3, resid3,main="Residual Plot") # against predicted values
abline(h=0, col="red")
mean(resid3)
sd(resid3) # standard deviation of residuals
mean(abs(resid3)>mean(resid3)+2*sd(resid3))
#These look much better, although the residuals for Saudi Arabia and Libya are
#still more than three standard deviations away from the mean of the residuals.

rows <-c(1:nrow(yhat3))
#plot these predicted values against the actual outcome.
plot(rows,yhat3,col="green",
     main = "prediction and actual values comparsion's plot",
     xlab = "# of data: green dots are prediction values, purple dots are actual values",
     ylab = "prediction and actual values")
par(new=TRUE)
plot(rows,jitter(data3_jags$y),col="purple",
     main = "prediction and actual values comparsion's plot",
     xlab = "# of data: green dots are prediction values, purple dots are actual values",
     ylab = "prediction and actual values")

qqnorm(resid3)

df1[rownames(df1)[order(abs(resid3), decreasing=TRUE)[1:5]],] # largest absolute residual value
##########Density plot for actual and predictive values###############
d<-density(df1$logmedv)#actual
d1 <-density(yhat3) #predictive
plot(d,xlab="log transformation for House Price's median value",col="green",main="Density plot for actual and predictive values")
par(new=TRUE)
plot(d1,xlab="log transformation for House Price's median value",col="purple",main="Density plot for actual and predictive values")
labels = c("Actual values", "Predictive values")
legend("topright",legend = labels, col=c("green", "purple"), lty=c(1,1), cex=0.8)
