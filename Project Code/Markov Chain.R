# Markov chain:random walk(continuous)
#p(X(t+1)|X(t)=x(t)) = N(x(t),1)
# the probability distribution for the next state is Normal with variance 1 and mean = current state
set.seed(34) 
n = 100 
x = numeric(n) 
for (i in 2:n) {
  x[i] = rnorm(1, mean=x[i-1], sd=1.0) 
  } 
plot.ts(x) 

#discrete markov chain: assume the transition probability does't change
#so transition probability can be from state 1 to state 2: Q1, state1 to state3: Q2

Q = matrix(c(0.0, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0),  nrow=5, byrow=TRUE) 
Q %*% Q # Matrix multiplication in R. This is Q^2. 
#Therefore, if your secret number is currently 1, the probability that the number will be 3 two steps from now is .25. 

#p(X(t+h)|X(t)) where h is a big number:
Q5 = Q %*% Q %*% Q %*% Q %*% Q # h=5 steps in the future 
round(Q5, 3) 
Q10 = Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q # h=10 steps in the future 
round(Q10, 3)

Q30 = Q 
for (i in 2:30) {
  Q30 = Q30 %*% Q 
  } 
round(Q30, 3) # h=30 steps in the future 

#stationary distribution of a chain: 
# is the initial state distribution for which performing 
#a transition will not change the probability of ending up in any given state. 
c(0.2, 0.2, 0.2, 0.2, 0.2) %*% Q 

help("sample")
n = 5000 
x = numeric(n) # empty numeric array with size n
x[1] = 1 # fix the state as 1 for time 1 
for (i in 2:n) {
  x[i] = sample.int(5, size=1, prob=Q[x[i-1],]) # draw the next state from the intergers 1 to 5 with proba bilities from the transition matrix Q, based on the previous value of X. 
} 
#sample takes a sample of the specified size from the elements of x 
#weighted by probability matrix
table(x) / n 

# Continuous Markov Chain example:
# p(X(t+1)|X(t)=x(t)) = N(phi*x(t),1), -1<phi<1

set.seed(100) 
n = 5000
x = numeric(n) # empty numeric array with size n
phi = -0.6 
for (i in 2:n) { 
  x[i] = rnorm(1, mean=phi*x[i-1], sd=1.0) 
  } 
plot.ts(x) 

# a histogram of our chain
#compare that with the theoretical stationary distribution
hist(x, freq=FALSE) 
curve(dnorm(x, mean=0.0, sd=sqrt(1.0/(1.0-phi^2))), col="red", add=TRUE) 
legend("topright", legend="theoretical stationary\ndistribution", col="red", lty=1, bty="n") 
# This chain has reached the stationary distribution since for different n states, the distribution states the same. 

#The first row correspond to X=0 (player A not playing) 
#while the second row  correspond to X=1 (player A playing).
Q = matrix(c(0,1,0.3,0.7),nrow=2, byrow=TRUE)
#suppose that the first game is between Players B and C
#the probability that Player A will play in Game 4:
c(1,0)%*% Q %*% Q %*% Q 
# stationary distribution for X in the chess example:
Q10 = Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q
#If the chain starts in the stationary distribution, 
#the probability of Player A playing in the next game
#the game after that, and so forth, is always this stationary probability.

