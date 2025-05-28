#########################################################################################################
######## Run Stage one of two-stage MCMC algorithm to jointly model the environmental variables #########
#########################################################################################################

set.seed(810) ## Why not...

###### will run as an array job; run grid cell K in parallel
args = Sys.getenv('SLURM_ARRAY_TASK_ID')
K = as.numeric(args[1])
outfile = paste0("StageOneCovariateOut/CovariateOut",K,".Rda")

###### load data
data <- read.csv("USDMData.csv") #### these csv files were generated when previously running ScalingValues.R.
scalingvalues <- read.csv("scalingvalues.csv")


## For this application, remove grid cell with no change in drought status
sub <- data
drop = which(sub$grid %in% c("N78","W98","GG14","WW88"))
sub <- sub[-drop,]
rm(data)

## only keep data from July 1, 2003 - 2012
sub <- sub[sub$time < 20220630,]
sub <- sub[sub$time > 20110101,]
sub <- sub[sub$grid == unique(sub$grid)[K],] #### extract just grid cell K

library(nimble) # For MCMC computation using NIMBLE.
library(coda) # For manipulation of MCMC results.

id.use = c(
  "time",
  "grid",
  "lon",
  "lat"
)

vars.use = c(
  "evp","soilm","tsoil"
)

wh.scaling <- which(scalingvalues$X %in% vars.use)
wh.data <- which(colnames(sub) %in% c(id.use, vars.use))

sub <- sub[,wh.data]
T = length(unique(sub$time)) - 13; T
I = length(unique(sub$grid)); I
J = dim(sub)[2] - 4
crdtm = data.frame("lon"=sub$lon[1:I],"lat"=sub$lat[1:I])


############ detrend the covariates by subtracting out a deterministic seasonal trend

## Define Fourier terms, 5 pairs seems like enough
x.sin1 = sin(2*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)
x.cos1 = cos(2*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)
x.sin2 = sin(4*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)
x.cos2 = cos(4*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)
x.sin3 = sin(6*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)
x.cos3 = cos(6*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)
x.sin4 = sin(8*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)
x.cos4 = cos(8*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)
x.sin5 = sin(10*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)
x.cos5 = cos(10*pi*as.numeric(as.factor(unique(sub$time)))/52.14286)

## Make design matrix
X = cbind(rep(1,I*(T+13)), x.sin1, x.cos1, x.sin2, x.cos2, x.sin3, x.cos3, x.sin4, x.cos4, x.sin5, x.cos5)
bp = length(X[1,])
colnames(X)[1]="int"

# Define response
Y = sub[(1:T),5:(4+J)]

M.Y = scalingvalues[wh.scaling,2]
SD.Y = scalingvalues[wh.scaling,3]

for (j in 1:J){
  Y[,j] = (Y[,j] - M.Y[j])/SD.Y[j]
}
head(Y)

## Fit OLS with X
## Save these OLS parameters so the trend can be added back in to forecast drought
## Keep residuals as data moving forward and model these detrended quantities.

ols = matrix(NA,J,bp)
Y2 = Y

for (j in 1:J){
  lm.hold = lm(Y[,j] ~ X[1:T,] - 1)
  ols[j,] = lm.hold$coeff
  Y2[,j] = lm.hold$residuals
}


## Stuff for NIMBLE
mod_data=list(Y2=Y2, X=X)
mod_constants=list(T=T, J=J, I=I, R=diag(J))
mod_inits=list(rho=rep(.75,J), SigInv = diag(J))

model_code=nimbleCode({
  ##prior distribution for Sigma
  SigInv[1:J,1:J] ~ dwish(R[1:J,1:J], df=J)
  
  ##prior distribution for rho, j indexes each element of beta
  for(j in 1:J){
    rho[j] ~ dunif(0,1)
  }
  
  ## Likelihood
  for (i in 1:I){
    Y2[i,1:J] ~ dmnorm(mu[i,1:J], prec = SigInv[1:J,1:J])
    for (j in 1:J){
      mu[i,j] <- 0
    }
    
    for(t in 2:T){
      Y2[((t-1)*I+i),1:J] ~ dmnorm(mu[((t-1)*I+i),1:J], prec = SigInv[1:J,1:J])
      for (j in 1:J){
        mu[((t-1)*I+i),j] <- rho[j]*(Y2[((t-1)*I+i - I),j])
      }
    }
    
    ## Prediction loop here
    for(t in (T+1):(T+1)){
      Ypred[((t-T-1)*I+i),1:J] ~ dmnorm(mu[((t-1)*I+i),1:J], prec = SigInv[1:J,1:J])
      for (j in 1:J){
        mu[((t-1)*I+i),j] <- rho[j]*(Y2[((t-1)*I+i - I),j])
      }
    }
    
    for(t in (T+2):(T+13)){
      Ypred[((t-T-1)*I+i),1:J] ~ dmnorm(mu[((t-1)*I+i),1:J], prec = SigInv[1:J,1:J])
      for (j in 1:J){
        mu[((t-1)*I+i),j] <- rho[j]*(Ypred[((t-T-1)*I+i - I),j])
      }
    }
  } ## closes location i
} ## closes nimble code
)

nimble_model <- nimbleModel(model_code, mod_constants,mod_data,mod_inits)
compiled_model <- compileNimble(nimble_model,resetFunctions = TRUE)
mcmc_conf <- configureMCMC(nimble_model,monitors=c('rho','SigInv','Ypred'),control=list(adaptive=TRUE,scale=0.1,adaptInterval=100000,sliceMaxSteps=100000,maxContractions=5000000,sliceAdaptWidthMaxIter=0,sliceAdaptFactorMaxIter=0),useConjugacy = TRUE)

### first you have to remove the default sampler and then change it
names=array()
for (j in 1:J){
  names[j] = paste("rho[",j,"]", sep="")
  mcmc_conf$removeSamplers(names[j])
  mcmc_conf$addSampler(target=names[j],type='slice')
}

nimble_mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(nimble_mcmc, project = nimble_model,resetFunctions = TRUE)

MCS=10000
samples_sub=runMCMC(compiled_mcmc,inits=mod_inits,
                    nchains = 1, nburnin=MCS/2,niter = MCS,samplesAsCodaMCMC = TRUE,thin=5,
                    summary = FALSE, WAIC = FALSE, progressBar=TRUE)

SimMCMC = list(samples_sub, X, sub, crdtm, M.Y, SD.Y, ols, Y2)

save(SimMCMC, file=outfile)




