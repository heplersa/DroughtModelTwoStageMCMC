##################################################################################
##### Integrate Covariate and Drought model output ########
############# Forecast next 13 weeks ######################
#########################################################

######## Will need the covariates from the last observed time period to predict Z for first new time period
data = read.csv('USDMData.csv')
scalingvalues = read.csv('scalingvalues.csv')
data = data[-which(data$grid %in% c("N78","W98","GG14","WW88")),] ### remove these 4 locations whose drought level never changes because they are over water
data = data[which(data$time>20110101),]
training <- data[which(data$time<20220400),] ### data to model

training$timeID = as.numeric(as.factor(training$time))
gridID = training$grid[which(training$timeID==1)]

All = length(gridID) ## scientific runs for all locations
Tobs = nrow(training)/All
## Q = length(list.files("StageOneOutput")) ## For testing a subset
Q = All ## uncomment when fitting model to all locations instead.

#############################################################
####### Create the design matrix

vars <- c("evp","soilm","tsoil")
Xcov = cbind(training[,vars])
J = ncol(Xcov)

Xm = rep(0,J)
Xs = rep(0,J)
for(j in 1:J){
  Xm[j] = scalingvalues$means[which(scalingvalues$X==vars[j])]
  Xs[j] = scalingvalues$sds[which(scalingvalues$X==vars[j])]
}

Xfull = matrix(NA,Q*Tobs,J)
for(j in 1:(J)){
  Xfull[,j] = (Xcov[,j]-Xm[j])/Xs[j]
}

#### reshape Xfull to be array 
X = array(NA,c(Q,J,Tobs))
for(t in 1:Tobs){
  X[,,t]=Xfull[((t-1)*Q+1):(t*Q),]
}

#############################################################################################

######## Load the draws from the posterior predictive distribution of the environmental covariates

load('StageTwoCovariateOutput.Rda') #### output of MCMC for covariate model

Xpred = Ypred.out #### QxJx13x1000 draws for predicted values of the DETRENDED covariates
rm("Ypred.out")

### add back in the deterministic mean trend that was removed
### the coefficients of the deterministic model are saved in ols list 
sub <- data
rm(data)

sub <- sub[sub$time < 20220630,]
sub <- sub[sub$time > 20110101,]

## I = length(list.files("StageOneOutput")) ## For testing a subset
I = 3254 ## all locations
Tpred = 13
T = length(unique(sub$time))-Tpred

sub <- sub[which(sub$grid==sub$grid[1]),]

## Define Fourier terms to get deterministic mean trend for covariates
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

X.cov = cbind(rep(1,(T+13)), x.sin1, x.cos1, x.sin2, x.cos2, x.sin3, x.cos3, x.sin4, x.cos4, x.sin5, x.cos5)

#### add the deterministic trend back to Xpred

base <- array(NA,c(I,J,Tpred))
for(i in 1:I){
  base[i,,] <- t(X.cov[(T+1):(T+Tpred),]%*%t(ols[[i]]))
}

#### draws from the posterior predictive distribution of the standardized covariates after adding the deterministic mean trend back
M=dim(Xpred)[4]
for(i in 1:M){
  Xpred[,,,i] <- base+Xpred[,,,i]
}

#########################################################################

###### load in samples from drought model output 

######### load drought model output
load('StageTwoOutput.Rda')

#### Check alignment of the covariate and drought output based on the length of the two chains.
#### only keep a subset of particles from the longer chain.  
M.drought = dim(beta.out)[3]
if (M.drought < M) {
  keep <- round(M.drought * seq(0,1,by=1/M),0)
  beta.out = beta.out[,,keep]
  rho.Z.out = rho.Z.out[,keep]
  Z.out = Z.out[,keep] ### samples for the very last observed time period Z
  sigma.sq.out = sigma.sq.out[,keep]
}
if (M.drought > M) {
  keep <- round(M.drought * seq(0,1,by=1/M),0)
  beta.out = beta.out[,,keep]
  rho.Z.out = rho.Z.out[,keep]
  Z.out = Z.out[,keep] ### samples for the very last observed time period Z
  sigma.sq.out = sigma.sq.out[,keep]
}



#### for each posterior sample of (theta, Xpred), generate Z.pred then generate Y.pred
Z.pred = array(NA,c(Q,13,M))
Y.pred = array(NA,c(Q,13,M))

alpha = c(-Inf,0,1,2,3,4,Inf)

for(m in 1:M){
  
  #### predict for time T+1
  Z.pred[,1,m] = rnorm(Q,beta.out[,1,m]+rowSums(Xpred[,,1,m]*beta.out[,2:(J+1),m])+rho.Z.out[,m]*(Z.out[,m]-beta.out[,1,m]-rowSums(X[,,Tobs]*beta.out[,2:(J+1),m])),sqrt(sigma.sq.out[,m]))
  Y.pred[,1,m] = cut(Z.pred[,1,m],alpha,labels=FALSE,ordered_result=TRUE)-1
  
  #### predict for t=T+2,...,T+13
  for(t in 2:Tpred){
    Z.pred[,t,m] = rnorm(Q,beta.out[,1,m]+rowSums(Xpred[,,t,m]*beta.out[,2:(J+1),m])+rho.Z.out[,m]*(Z.pred[,(t-1),m]-beta.out[,1,m]-rowSums(Xpred[,,(t-1),m]*beta.out[,2:(J+1),m])),sqrt(sigma.sq.out[,m]))
    Y.pred[,t,m] = cut(Z.pred[,t,m],alpha,labels=FALSE,ordered_result=TRUE)-1
  }
  
  print(paste(m))
  
}

save(Z.pred,Y.pred,Xpred,file="ForecastOutput.Rda")





