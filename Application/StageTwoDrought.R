############################################################################
##### Run stage two of the 2-stage Bayesian model for drought #############
###########################################################################
library(invgamma)
library(mvtnorm)

st <- Sys.time()

###### load in all samples from the stage one model
###### consolidate all stage one output into arrays 

load(paste('StageOneOutput/MCMCout.',1,'.Rda',sep=""))

n = 3254
T = ncol(MCMCout$Z)
M = nrow(MCMCout$Z)

Z.samp = array(NA,c(n,T,M))
sigma.sq.samp = matrix(NA,n,M)
beta.samp = array(NA,c(n,dim(MCMCout$beta)[2],M))
rho.Z.samp = matrix(NA,n,M)

rm('MCMCout')

for(q in 1:n){

	load(paste('StageOneOutput/MCMCout.',q,'.Rda',sep=""))

	 Z.samp[q,,] = t(MCMCout$Z)
  sigma.sq.samp[q,]=MCMCout$sigma.sq
  beta.samp[q,,]=t(MCMCout$beta)
  rho.Z.samp[q,]=MCMCout$rho.Z
  rm('MCMCout')

}


# save(Z.samp,sigma.sq.samp,beta.samp,rho.Z.samp,file="MCMCoutputCombined.Rda") # save combined output, if needed. 

#### convert rho.Z to gamma samples
gamma.samp = log(rho.Z.samp/(1-rho.Z.samp))

##########################################################################
####### Set up MCMC for the second stage

###### set initial values
###### take one draw for each site-specific parameter
Z = Z.samp[,,M]
sigma.sq = sigma.sq.samp[,M]
beta = beta.samp[,,M]
gamma = gamma.samp[,M]
  
bp=dim(beta.samp)[2]

#### initial values for overall parameters
tausq.b = rep(1,bp)
tausq.g = 1

##### need adjacency matrix for grid cells and number of neighbors of each...
##### need adjacency matrix for grid cells and number of neighbors of each...

data <- read.csv("USDMDataAvg.csv")
data <- data[,-1] ## remove leading column
data = data[-which(data$grid %in% c("N78","W98","GG14","WW88")),] ### remove these 4 locations whose drought level never changes because they are over water


loc = data[which(data$time==data$time[1]),c("lon","lat")]
rm('data')
D = as.matrix(dist(loc))
A = 1*(D>0 & D<=sqrt(.5)) ### regular grid; queens adjacency means distances is <= sqrt(0.5)
numnns = rowSums(A)
D = matrix(0,n,n)
diag(D)=numnns



M.iter = 55000
M.burn = 5000
M.thin=10
M.out = (M.iter-M.burn)/M.thin

##### create storage matrices
### note we only need to save the last Z for forecasting purposes
Z.out = matrix(NA,n,M.out)
sigma.sq.out = matrix(NA,n,M.out)
beta.out = array(NA,c(n,bp,M.out))
rho.Z.out = matrix(NA,n,M.out)
tausq.b.out = matrix(NA,bp,M.out)
tausq.g.out = rep(NA,M.out)

progress_bar = txtProgressBar(min=0, max=M.iter, style = 3, char="=")

for(m in 1:M.iter){
  
  #### update tausq.b for each covariate
  for(p in 1:bp){
    tausq.b[p] = rinvgamma(1,0.5+n/2,0.5+1/2*t(beta[,p])%*%(D-A)%*%beta[,p])   
  }
  
  #### update tausq.g
  tausq.g = rinvgamma(1,0.5+n/2,0.5+1/2*t(gamma)%*%(D-A)%*%gamma)

  
  #### update the site-specific parameters for each location
  #### update full vector of site-specific parameters jointly
  for(i in 1:n){
   mi = sample(1:M,1)
    gamma.new = gamma.samp[i,mi]
    gammam = 1/numnns[i]*A[i,]%*%gamma
    betap.new = beta.samp[i,,mi]
    betam = 1/numnns[i]*A[i,]%*%beta
    ### make gamma prior ICAR
    R.new = dnorm(gamma.new,gammam,sd=sqrt(tausq.g/numnns[i]),log=TRUE)+sum(dnorm(betap.new,betam,sqrt(tausq.b/numnns[i]),log=TRUE))-dlogis(gamma.new,0,1,log=TRUE)-sum(dnorm(betap.new,0,sd=3,log=TRUE))
    R.old = dnorm(gamma[i],gammam,sd=sqrt(tausq.g/numnns[i]),log=TRUE)+sum(dnorm(beta[i,],betam,sqrt(tausq.b/numnns[i]),log=TRUE ))-dlogis(gamma[i],0,1,log=TRUE)-sum(dnorm(beta[i,],0,sd=3,log=TRUE))
    if(log(runif(1))<(R.new-R.old)){
      gamma[i] = gamma.new
      sigma.sq[i] = sigma.sq.samp[i,mi]
      Z[i,] = Z.samp[i,,mi]
      beta[i,] = betap.new
    }
  }
  
  
  
  if(m>=M.burn & m/M.thin==floor(m/M.thin)){
    Z.out[,(m-M.burn)/M.thin] = Z[,T]
    sigma.sq.out[,(m-M.burn)/M.thin] = sigma.sq
    beta.out[,,(m-M.burn)/M.thin] = beta
    rho.Z.out[,(m-M.burn)/M.thin] = exp(gamma)/(1+exp(gamma))
    tausq.b.out[,(m-M.burn)/M.thin] = tausq.b
    tausq.g.out[(m-M.burn)/M.thin] = tausq.g
  }
  
  setTxtProgressBar(progress_bar, value = m)
  
}
close(progress_bar)
time.out <- Sys.time()-st

save(Z.out,sigma.sq.out,beta.out,rho.Z.out,tausq.b.out,tausq.g.out,time.out,file="Stage2Output.Rda")
