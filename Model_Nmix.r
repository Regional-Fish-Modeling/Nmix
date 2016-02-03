model{
  for(i in 1:nSites){
    for(t in 1:nYears){
      N[i,t] ~ dpois(lambda[i,t])
      log(lambda[i,t]) <- alpha[i] + 
        b[1,i]*fall.prcp[i,t] + b[2,i]*winter.prcp[i,t] + b[3,i]*spring.prcp[i,t] +
        b[4,i]*fall.tmean[i,t] + b[5,i]*winter.tmean[i,t] + b[6,i]*spring.tmean[i,t]      
    }
  }
  
  # Priors: intercepts
  for(i in 1:nSites){
    alpha[i] ~ dunif(-5,5)
  }
  
  # Priors: slopes
  for(h in 1:nCovs){
    for(i in 1:nSites){
      b[h,i] ~ dnorm(mu.b[h,i], tau.b[h])
      mu.b[h,i] <- g.0[h] + g.1[h]*elev[i]
    }  
    g.0[h] ~ dnorm(0, 0.01)
    g.1[h] ~ dnorm(0, 0.01)  
    tau.b[h] <- pow(sigma.b[h], -2)
    sigma.b[h] ~ dunif(0,3)
    sigma2.b[h] <- pow(sigma.b[h], 2)
  }
  
  # Detection
  for(i in 1:nSites){
    for(t in 1:nYears){
      y[i,t,1] ~ dbin(p[i,t], N[i,t])
      y[i,t,2] ~ dbin(p[i,t]*(1-p[i,t]), N[i,t])
      y[i,t,3] ~ dbin(p[i,t]*(1-p[i,t])*(1-p[i,t]), N[i,t])
      p[i,t] <- 1/(1 + exp(-lp.lim[i,t]))
      lp.lim[i,t] <- min(999, max(-999, lp[i,t]))
      lp[i,t] <- p.mu + p.b1*prcp7day[i,t] + p.b2*sampday[i,t]
    }
  }
  
  # Priors: detection
  p.mean ~ dunif(0.1, 0.9)
  p.mu <- log(p.mean/(1-p.mean))
  p.b1 ~ dnorm(0, 0.37)T(-3,3)
  p.b2 ~ dnorm(0, 0.37)T(-3,3)
}