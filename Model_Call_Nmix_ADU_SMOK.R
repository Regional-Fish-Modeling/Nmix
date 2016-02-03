### Brook trout ADULT abundance model for GSMNP (SMOK) adapted from Kanno et al. (2016) Freshwater Biology
# NP Hitt <nhitt@usgs.gov> 3 Feb 2016

# Load working directory & libraries
setwd("C:/Users/nhitt/Desktop/Nmix_SMOK")
memory.limit(size=229355) # 7x the baseline memory limit (32765) 
library(reshape2)
library(rjags)
library(plyr)
library(ggplot2)
library(arm)
library(boot)
library(coda)

# Fish count data
load("Data_SMOK_FishCountAr.RData")
# Seasonal climate data
load("Data_SMOK_SeasonalClimateStd.RData") # Standardized to site-specific mean=0 and sd=1
# Site covariate data
load("Data_SMOK_SiteCovsStd.RData") # Standardized to dataset-level mean=0 and sd=1
# Detection covariate data
load("Data_SMOK_DetectionCovsStd.RData") # Standardized to dataset-level mean=0 and sd=1

# Set up data structure
nSites=70; nYears=31; nCovs=6; nPasses=3

# Adding 1 fish to abundances to avoid log(0) problem
ADUFish <- ADUFish + 1

# Bundle data for Adult model
dat <- list(nSites=nSites, nYears=nYears, nCovs=nCovs, y=ADUFish,
            fall.prcp=FallPrcpStd, winter.prcp=WinterPrcpStd, spring.prcp=SpringPrcpStd, 
            fall.tmean=FallTmeanStd, winter.tmean=WinterTmeanStd, spring.tmean=SpringTmeanStd, 
            prcp7day=prcp7day.std, sampday=sampday.std, elev=elev.std)

# Set initial values
init <- function() list(N=array(500, dim=c(nSites, nYears)),
                        alpha=array(runif(nSites,-5,5), dim=(nSites)),
                        b=array(rnorm(nCovs*nSites,0,1), c(nCovs,nSites)),
                        g.0=array(rnorm(nCovs,0), c(nCovs)),
                        g.1=array(rnorm(nCovs,0), c(nCovs)),                                             
                        sigma.b=array(runif(nCovs,0,3), c(nCovs)),
                        p.mean=0.5, 
                        p.b1=array(rnorm(1,0), 1),
                        p.b2=array(rnorm(1,0), 1))

# Run JAGS
set.seed(1234)
Nadapt <- 1000
StageBurnin <- jags.model("Model_Nmix.R", dat, init, n.chains=3, n.adapt=Nadapt)
Niter=1000
Nthin=10
pars <- c("N","g.0","p")
out1 <- coda.samples(StageBurnin, pars, n.iter=Niter, n.thin=Nthin)
out2 <- jags.samples(StageBurnin, pars, n.iter=Niter, n.thin=Nthin)
save(out1, out2, file="ResultsRaw_ADU_SMOK_3Feb16.RData")

# Calculate estimated abundance
N.est <- p.est <- array(NA, dim=c(nSites,nYears))
for(i in 1:nSites){
  for(t in 1:nYears){
    N.est[i,t] <- median(out2$N[i,t,Nadapt,1:3])
    p.est[i,t] <- median(out2$p[i,t,Nadapt,1:3])
  }
}

y.est <- array(NA, dim=c(nSites,nYears,nPasses))
y.est[,,1] <- N.est*p.est
y.est[,,2] <- N.est*(1-p.est)*p.est
y.est[,,3] <- N.est*(1-p.est)*(1-p.est)*p.est
save(N.est, y.est, file="Result_ADUEst_SMOK_3Feb16.RData")
# gelman.diag(out1)

#########################################################
# Plot histograms for seasonal climate covariates
pdf("Result_EnvHist_ADU_SMOK_3Feb16.pdf")
varlist <- c("Fall precip","Winter precip","Spring precip","Fall temp","Winter temp","Spring temp")
for(i in 1:length(varlist)){
  hist(out2$g.0[i,,], col="black", main=paste(varlist[i]), xlab="Posterior probability")
  abline(v=quantile(out2$g.0[i,,], probs=c(0.025)), col="red")
  abline(v=quantile(out2$g.0[i,,], probs=c(0.975)), col="red")
  abline(v=quantile(out2$g.0[i,,], probs=c(0.05)), col="blue")
  abline(v=quantile(out2$g.0[i,,], probs=c(0.95)), col="blue")
}
dev.off()
