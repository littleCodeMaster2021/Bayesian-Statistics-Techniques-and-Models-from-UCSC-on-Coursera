# 
# ppd = function(a, b, y_star) {
# b^a/(b+y_star)^(a+1)
# }
#   
# y_star=seq(from=0,to=120,by=.1) 
# 
# plot(y_star,ppd(9,390,y_star),type="l")
# 
# 
# 
# 
# a = 3
# b = 200
# Gamma = rgamma(n=1000, shape=a, rate=b)
# IG = 1/Gamma
# mean(IG)
# Normal = rnorm(n=1000,mean = 500,sd = sqrt(IG/0.1)) 
# posterior = rnorm(n=1000,mean = Normal,sd = IG) 
# mean(posterior)

z <- rgamma(1000, shape=16.5, rate=6022.9)
sig2 <- 1/z
mu <- rnorm(1000, mean=609.3, sd=sqrt(sig2/27.1))
#95% equal-tailed credibility for mu: quantile/percentiles of the simulated value;
quantile(x=mu, probs=c(0.025, 0.975))

#http://www.randomservices.org/random/data/Challenger2.txt
# 23 previous space shuttle launches before the Challenger disaster
# T is the temperature in Fahrenheit, I is the O-ring damage index

oring=read.table("http://www.randomservices.org/random/data/Challenger2.txt",header=T)
#use the attach command to be to access the labels. 
attach(oring)# with label T and I
#note: masking T=TRUE

plot(T,I)

#lm: linear model
        #lm(y~x)
oring.lm=lm(I~T)
summary(oring.lm)

# add fitted line to scatterplot
lines(T,fitted(oring.lm))    

# 95% posterior interval for the slope
#qt gives quantile function: qt(p, df, ncp, lower.tail = TRUE, log.p = FALSE)
-0.24337 - 0.06349*qt(p=.975,df=21)
-0.24337 + 0.06349*qt(.975,21)
# note that these are the same as the frequentist confidence intervals
#when we're using the standard reference prior for Bayesian analysis. 

# the Challenger launch was at 31 degrees Fahrenheit
# how much o-ring damage would we predict?
# y-hat=slop*t+intercept
18.36508-0.24337*31  #10.82052
coef(oring.lm) #get coefficience information
coef(oring.lm)[1] + coef(oring.lm)[2]*31  
  
# posterior prediction interval (same as frequentist) for y*|y (95% credibility interval)
predict(oring.lm,data.frame(T=31),interval="predict")
10.82052-2.102*qt(.975,21)*sqrt(1+1/23+((31-mean(T))^2/22/var(T)))

# posterior probability that damage index(for slope) is greater than zero
# using the reference prior, allows us to use all the standard regression tools. 
1-pt((0-10.82052)/(2.102*sqrt(1+1/23+((31-mean(T))^2/22/var(T)))),21)

##multiple regression example. 
#http://www.randomservices.org/random/data/Galton.txt
# Galton's seminal data on predicting the height of children from the 
# heights of the parents, all in inches

heights=read.table("http://www.randomservices.org/random/data/Galton.txt",header=T)
attach(heights)# add data frame with labels
names(heights)

#lm is used to fit linear models. It can be used to carry out regression, 
#single stratum analysis of variance and analysis of covariance
pairs(heights)
summary(lm(Height~Father+Mother+Gender+Kids))
summary(lm(Height~Father+Mother+Gender))
heights.lm=lm(Height~Father+Mother+Gender)

resid(heights.lm) #List of residuals
plot(density(resid(heights.lm))) #A density plot
qqnorm(resid(heights.lm)) # A quantile normal plot - good for checking normality
qqline(resid(heights.lm))

# each extra inch taller a father is is correlated with 0.4 inch extra
#height in the child
# each extra inch taller a mother is is correlated with 0.3 inch extra
#height in the child
# a male child is on average 5.2 inches taller than a female child
# 95% posterior interval for the the difference in height by gender
5.226 - 0.144*qt(.975,894)
5.226 + 0.144*qt(.975,894)

# posterior prediction interval (same as frequentist) for y*|y (95% credibility interval)
predict(heights.lm,data.frame(Father=68,Mother=64,Gender="M"),interval="predict")
predict(heights.lm,data.frame(Father=68,Mother=64,Gender="F"),interval="predict")

#########################################################################
dat=read.table("http://www.stat.ufl.edu/~winner/data/pgalpga2008.dat",header=F)
datF =  subset(dat, V3==1, select=1:2) #"V3" is the name of the third column, ONLY has V3==1: female golfer (on the LPGA tour) 
datM =  subset(dat, V3==2, select=1:2)
attach(datF)
plot(datF)
datF.lm=lm(V2~V1)
summary(datF.lm)#posterior mean estimate of the slope parameter
#Find V2 and 95% posterior predictive interval for V2 When V1 = 260:
predict(datF.lm,data.frame(V1=260),interval="predict")

dat=read.table("http://www.stat.ufl.edu/~winner/data/pgalpga2008.dat",header=F)
dat$V3[dat$V3 == 1] <- 0
dat$V3[dat$V3 == 2] <- 1
summary(dat)
attach(dat)
dat.lm=lm(V2~V1+V3)
summary(dat.lm)
plot(fitted(dat.lm), residuals(dat.lm))
#The residuals appear to be random and lack any patterns or trends. 
#However, there is at least one outlier (extreme observation) that we may want to investigate.
#Outliers can strongly influence model estimates.

#replace B with b in column nm
junk <- data.frame(x <- rep(LETTERS[1:4], 3), y <- letters[1:12])
colnames(junk) <- c("nm", "val")
typeof(junk$val)# integer

#Easier to convert nm to characters and then make the change:
junk$nm <- as.character(junk$nm)
junk$nm[junk$nm == "B"] <- "b"
#EDIT: And if indeed you need to maintain nm as factors, add this in the end:
junk$nm <- as.factor(junk$nm)
#convert back to factor in case that was an important part of the data structure.

within(junk, levels(nm)[levels(nm) == "B"] <- "b")# method 1
# add an additional level "b" to the factor attributes.
levels(junk$nm) <- c(levels(junk$nm), "b") # method 2: levels find all letters in the column
junk$nm[junk$nm == "B"] <- "b"
