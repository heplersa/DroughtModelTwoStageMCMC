##########################################################################
##### Run stage two of the two-stage MCMC for the covariate model ########
##### save forecasted values of the covariates for the next 13 weeks #####
#########################################################################
library(invgamma)
library(mvtnorm)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## Load all of the individual files and consolidate into array
I = 3254
names = paste0("StageOneCovariateOut/CovariateOut",1:I, ".Rda")
StageOne = list()

for (i in 1:I){
  load(names[i])
  StageOne[[i]] = SimMCMC
  rm('SimMCMC')
}

J = length(StageOne[[1]][[5]])
bp = dim(StageOne[[1]][[2]])[2]
M = length(StageOne[[1]][[1]][,1])

lon=rep(0,I)
lat=rep(0,I)
for (i in 1:I){
  lon[i]=as.numeric(StageOne[[i]][[4]][1])
  lat[i]=as.numeric(StageOne[[i]][[4]][2])
}

SimMCMC = StageOne
rm(StageOne)
rm(names)
L = 13 ## length of prediction period

## Make combined OLS estimates here, save in stage two output
ols <- list()
for (i in 1:I){
  ols[[i]] <- SimMCMC[[i]][[7]]
}


##### reorg output so in matrices instead of list
sigma.inv.samp = array(NA,c(I,J*J,M))
rho.samp = array(NA,c(I,J,M))
Ypred.samp = array(NA,c(I,J,L,M))

for(i in 1:I){
  temp = SimMCMC[[i]][[1]]
  wh.rho = which(substr(colnames(temp),1,3) == "rho")
  rho.samp[i,,]=t(temp[,wh.rho])

  wh.siginv = which(substr(colnames(temp),1,6) == "SigInv")
  sigma.inv.samp[i,,] = t(temp[,wh.siginv])

  for (j in 1:J){
    wh.Ypred.a = which(substr(colnames(temp),1,5) == "Ypred")
    wh.Ypred.b = which(substrRight(colnames(temp),2) == paste(j,"]",sep=""))
    wh.Ypred = intersect(wh.Ypred.a, wh.Ypred.b)
    Ypred.samp[i,j,,] = t(temp[,wh.Ypred]) ## Update here next
  }
}

#### convert rho.Z to gamma samples
gamma.samp = log(rho.samp/(1-rho.samp))
rm(rho.samp)

##### need adjacency matrix for grid cells and number of neighbors of each...
loc = cbind(lon, lat)
Di = as.matrix(dist(loc))
A = 1*(Di>0 & Di<=sqrt(.5))
numnns = rowSums(A)
D = matrix(0,I,I)
diag(D)=numnns

## Clean up
rm(SimMCMC)

###### set initial values to the last draw
sigma.inv = sigma.inv.samp[,,M]
gamma = gamma.samp[,,M]
Ypred = Ypred.samp[,,,M]

#### initial values ICAR precisions
tausq.gamma = rep(1,J)

##### store all final draws
M.iter = 10000
M.burn = 5000
M.thin = 5
M.out = (M.iter-M.burn)/M.thin

sigma.inv.out = array(NA,c(I,J*J,M.out))
rho.out = array(NA,c(I,J,M.out))
Ypred.out = array(NA,c(I,J,L,M.out))
tausq.gamma.out = array(NA,c(J,M.out))

progress_bar = txtProgressBar(min=0, max=M.iter, style = 3, char="=")
st<-Sys.time()
for(m in 1:M.iter){

  #### update tausq.rho
  for (j in 1:J){
    tausq.gamma[j] = rinvgamma(1,0.5+I/2,0.5+1/2*t(gamma[,j])%*%(D-A)%*%gamma[,j])
  }

  #### update individual site-specific parameters 
  for(i in 1:I){

    R.new=rep(0,J)
    R.old=rep(0,J)
    mi = sample(1:M,1)
    gamma.new = gamma.samp[i,1:J,mi]
    Ypred.new = Ypred.samp[i,1:J,,mi]
    gammam = 1/numnns[i]*A[i,]%*%gamma[,]
    sigma.inv.new = sigma.inv.samp[i,,mi]

    for (j in 1:J){ ## loops over each environmental covariate, we will sum terms for a single Metropolis update
      R.new[j] = dnorm(gamma.new[j],gammam[j],sd=sqrt(tausq.gamma[j]/numnns[i]),log=TRUE)-dlogis(gamma.new[j],0,1,log=TRUE)
      R.old[j] = dnorm(gamma[i,j],gammam[j],sd=sqrt(tausq.gamma[j]/numnns[i]),log=TRUE)-dlogis(gamma[i,j],0,1,log=TRUE)
    } ## closes variable j

    R.new.T = sum(R.new)
    R.old.T = sum(R.old)
    if(log(runif(1))<(R.new.T-R.old.T)){ ## Everything that shows up in the stage one model!!
      gamma[i,1:J] = gamma.new
      Ypred[i,1:J,] = Ypred.new ## Update Ypred here?  Yes.
      sigma.inv[i,1:(J*J)] = sigma.inv.new
    }
  } ## closes location i

  if(m>=M.burn & m/M.thin==floor(m/M.thin)){
    sigma.inv.out[,,(m-M.burn)/M.thin] = sigma.inv
    Ypred.out[,,,(m-M.burn)/M.thin] = Ypred
    rho.out[,,(m-M.burn)/M.thin] = exp(gamma)/(1+exp(gamma))
    tausq.gamma.out[,(m-M.burn)/M.thin] = tausq.gamma
  }

  setTxtProgressBar(progress_bar, value = m)
}
close(progress_bar)
Sys.time()-st
save(sigma.inv.out, Ypred.out, rho.out, tausq.gamma.out, ols, lat, lon, file="StageTwoCovariateOutput.Rda")





