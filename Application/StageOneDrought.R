####################################################################
##### Code Stage One of the MCMC algorithm for Drought Model #######
#####################################################################

#################################################################
##### First, process the data and get means, SDs of covariates

source('ScalingValues.R') #### this script loads in the data, transforms some of the covariates, and computes the means and SDs of covariates

data = data[-which(data$grid %in% c("N78","W98","GG14","WW88")),] ### remove these 4 locations whose drought level never changes because they are over water

data = data[which(data$time>20110101),]

training <- data[which(data$time<20220400),] ### data to model

training$timeID = as.numeric(as.factor(training$time))
gridID = training$grid[which(training$timeID==1)]

Q = length(gridID)
Tobs = nrow(training)/Q

#############################################################
####### Create the design matrix

vars <- c("evp","soilm","tsoil")
Xcov = cbind(training[,vars])
J = ncol(Xcov)

Xm = rep(0,J)
Xs = rep(0,J)
for(j in 1:J){
  Xm[j] = scalingvalues$means[which(rownames(scalingvalues)==vars[j])]
  Xs[j] = scalingvalues$sds[which(rownames(scalingvalues)==vars[j])]
}


Xfull = matrix(NA,Q*Tobs,J)
for(j in 1:(J)){
  Xfull[,j] = (Xcov[,j]-Xm[j])/Xs[j]
}

training$droughtID = factor(training$drought,levels=c("0","D0","D1","D2","D3","D4"))
yfull = as.numeric(training$droughtID)-1

rm('select','Xcov','training','holdout','scalingvalues')


#####################################################

#### now can do MCMC for grid cell q
#### run as an array job so MCMC for each q in parallel

args = Sys.getenv('SLURM_ARRAY_TASK_ID')
q = as.numeric(args[1])

II = seq(q,(Q*Tobs),by=Q)
X = cbind(rep(1,Tobs),Xfull[II,])
Y = yfull[II]
  
D = 5 ### 6 total drought levels


### set initial values for MCMC
alpha = c(-Inf,0,1,2,3,4,Inf)
tau.z = 1
Z = 1*(Y - 0.5)

fit = arima(Z,order=c(1,0,0),xreg=X[,2:4]) #### use AR(1) fit to get initial values
rho.Z = fit$coef[1]
beta = fit$coef[2:5]
if(is.na(sum(beta))){
  beta[1] = Z[1]
  beta[2:ncol(X)]=0
  rho.Z=.98
}


bp = length(beta)

#### Set up MCMC
M.iter = 100000
M.burn = 20000
M.thin = 8

library(nimble)
library(coda)


st <- Sys.time()


mod_data=list(Y=Y, X=X)
mod_constants=list(Tobs=Tobs, bp=bp, cut=c(0,1,2,3,4))
mod_inits=list(beta=beta, rho.z=rho.Z, tau.z = tau.z, Z=1*(Y-.5))


model_code=nimbleCode({
  #Drought Variable
  for(t in 1:Tobs){
    Y[t] ~ dinterval(Z[t], cut[])
  }
  
  #Latent Gaussian Variable
  mu[1] <- inprod(X[1,1:bp],beta[1:bp])
  Z[1] ~ dnorm(mu[1], tau = tau.z)
  for(t in 2:Tobs){
    mu[t] <- inprod(X[t,1:bp],beta[1:bp]) +
      rho.z*(Z[(t-1)] - (inprod(X[(t-1),1:bp], beta[1:bp])))
    Z[t] ~ dnorm(mu[t], tau = tau.z)
  }
  
  
  ##prior distribution for beta, i indexes each element of beta
  for (b in 1:bp){
    beta[b] ~ dnorm(0, tau = 1/9)
  }
  
  ##prior distribution for tau.z
  tau.z ~ dgamma(0.01, 0.01)
  
  ##prior distribution for rho
  gamma ~ dlogis(0,1)
  rho.z <- exp(gamma)/(1+exp(gamma))
  

} ## closes nimble code
)

nimble_model <- nimbleModel(model_code, mod_constants, mod_data, mod_inits)
compiled_model <- compileNimble(nimble_model,resetFunctions = TRUE)
mcmc_conf <- configureMCMC(nimble_model,monitors=c('beta','rho.z','tau.z','Z'),control=list(adaptive=TRUE,scale=0.1,adaptInterval=100,sliceMaxSteps=100000,maxContractions=100000,sliceWidth=1),useConjugacy = TRUE)

mcmc_conf$removeSamplers("Z") ##### change to a slice sampler and update vectors of length 100 jointly
mcmc_conf$addSampler(target="Z[1:100]",type='AF_slice',control=list(sliceWidth=1))
mcmc_conf$addSampler(target="Z[101:200]",type='AF_slice',control=list(sliceWidth=1))
mcmc_conf$addSampler(target="Z[201:300]",type='AF_slice',control=list(sliceWidth=1))
mcmc_conf$addSampler(target="Z[301:400]",type='AF_slice',control=list(sliceWidth=1))
mcmc_conf$addSampler(target="Z[401:500]",type='AF_slice',control=list(sliceWidth=1))
mcmc_conf$addSampler(target="Z[501:587]",type='AF_slice',control=list(sliceWidth=1))


nimble_mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(nimble_mcmc, project = nimble_model,resetFunctions = TRUE)

samples=runMCMC(compiled_mcmc,inits=mod_inits,
                    nchains = 1, nburnin=M.burn,niter = M.iter,samplesAsCodaMCMC = TRUE,thin=M.thin,
                    summary = FALSE, WAIC = FALSE, progressBar=TRUE)


time.out = Sys.time()-st


zl = which(colnames(samples)=="Z[1]")
zu = which(colnames(samples)==paste("Z[",Tobs,"]",sep=""))
bl = which(colnames(samples) == "beta[1]")
bu = which(colnames(samples)== paste("beta[",bp,"]",sep=""))
rl = which(colnames(samples)=="rho.z")
tl = which(colnames(samples)=="tau.z")

MCMCout <- list("Z"=samples[,zl:zu],"sigma.sq"=1/samples[,tl],"beta"=samples[,bl:bu],"rho.Z"=samples[,rl])


save(MCMCout,time.out, file=paste("StageOneOutput/MCMCout.",q,".Rda",sep=""))



