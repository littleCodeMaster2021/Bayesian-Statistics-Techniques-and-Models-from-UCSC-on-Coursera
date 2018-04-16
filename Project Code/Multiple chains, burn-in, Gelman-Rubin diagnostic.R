# discard early iterations that do not appear to be coming from the stationary distribution. 
# Even if the chain appears to have converged early on, it is safer practice to discard an initial burn-in.
#Simulate multiple chains each with a different starting value:
set.seed(61)

nsim = 500
post1 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=15.0, cand_sd=0.4)
post1$accpt
post2 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=-5.0, cand_sd=0.4)
post2$accpt
post3 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=7.0, cand_sd=0.1)
post3$accpt
post4 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=23.0, cand_sd=0.5)
post4$accpt
post5 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=-17.0, cand_sd=0.4)
post5$accpt
pmc = mcmc.list(as.mcmc(post1$mu), as.mcmc(post2$mu), 
                as.mcmc(post3$mu), as.mcmc(post4$mu), as.mcmc(post5$mu))# combine all the objects using list
str(pmc)
coda::traceplot(pmc)
#This diagnostic statistic calculates the variability within chains, 
#comparing that to the variability between chains. If all chains have converged to the stationary
#distribution, the variability between chains should be relatively small, 
#and the potential scale reduction factor, should be close to one.
coda::gelman.diag(pmc)
coda::gelman.plot(pmc)
#But after about iteration 300, 
#the "shrink factor" is essentially one, indicating that by then, we have probably reached convergence. 