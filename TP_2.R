


#Function covariance matrix
matern1s2 = function(x,y,sigma2,l) {
  M = outer(x,y,"-")
  M = sigma2*exp(-(abs(M)/l))
  return(M)
}

matern3s2 = function(x,y,sigma2,l) {
  M = outer(x,y,"-")
  M = sigma2*(1+sqrt(6)*(abs(M)/l))*exp(-sqrt(6)*(abs(M)/l))
  return(M)
}

matern5s2 = function(x,y,sigma2,l) {
  M = outer(x,y,"-")
  M = sigma2*(1+sqrt(10)*(abs(M)/l)+(10/3)*(abs(M)/l)^2)*exp(-sqrt(10)*(abs(M)/l))
  return(M)
}

gaussien = function(x,y,sigma2,l) {
  M = outer(x,y,"-")
  M = sigma2*exp(-(abs(M)/l)^2)
  return(M)
}

#covariance parameters we use
sigma2=0.2
l=0.4

#The unknown function
f = function(x) {
  2 +cos(4*x) + x + x^2 - exp(x) + 0.3*sin(12*x)
}

#the large grid for the plots
grid = seq(from=0,to=1,length=1000)
fgrid = f(grid)
plot(grid,fgrid,type="l")

#The observation vector
xobs = c(0,0.2,0.4,0.6,0.8,1)
yobs = f(xobs)
points(xobs,yobs,col="red")


#The conditional mean function
pred = function(x,xobs,yobs,sigma2,l)  {
  yobs = matrix(nrow=length(yobs),ncol=1,data=yobs)
  Kobs=matern3s2(xobs,xobs,sigma2,l)
  KobsPred=matern3s2(xobs,x,sigma2,l)
  pred = t(KobsPred)%*%solve(Kobs)%*%yobs
  return(pred)
}

#add plot of the prediction on the grid
predgrid = pred(grid,xobs,yobs,sigma2,l)
lines(grid,predgrid,col="blue")

#The conditional variance function
condVar = function(x,xobs,yobs,sigma2,l)  {
  yobs = matrix(nrow=length(yobs),ncol=1,data=yobs)
  Kobs=matern3s2(xobs,xobs,sigma2,l)
  KobsPred=matern3s2(xobs,x,sigma2,l)
  Kpred = matern3s2(x,x,sigma2,l)
  condVar = Kpred - t(KobsPred)%*%solve(Kobs)%*%KobsPred
  return(diag(condVar))
}

#add plot of the 95 % confidence intervals on the grid
condVarGrid = condVar(grid,xobs,yobs,sigma2,l)
q975 = qnorm(p=0.975,mean=0,sd=1)
lines(grid, predgrid - q975*sqrt(condVarGrid) ,col="green")
lines(grid, predgrid + q975*sqrt(condVarGrid) ,col="green")

#The true needed length of the cable
trueLength = sum(
  sqrt(
    (grid[2:length(grid)] - grid[1:length(grid)-1])^2
                       + (fgrid[2:length(grid)] - fgrid[1:length(grid)-1])^2
    ) 
  ) 

#The length estimated by the metamodel mean
estLengthMean = sum(
  sqrt(
    (grid[2:length(grid)] - grid[1:length(grid)-1])^2
    + (predgrid[2:length(grid)] - predgrid[1:length(grid)-1])^2
  ) 
) 

#The conditional realization function
condReal = function(x,xobs,yobs,sigma2,l)  {
  yobs = matrix(nrow=length(yobs),ncol=1,data=yobs)
  Kobs=matern3s2(xobs,xobs,sigma2,l)
  KobsPred=matern3s2(xobs,x,sigma2,l)
  Kpred = matern3s2(x,x,sigma2,l)
  pred = t(KobsPred)%*%solve(Kobs)%*%yobs
  condCov = Kpred - t(KobsPred)%*%solve(Kobs)%*%KobsPred + 1e-8*diag(length(x))
  CcondCov = chol(condCov)
  pred + t(CcondCov)%*%matrix(nrow=length(x),ncol=1,data=rnorm(length(x)))
}

#add plot of a conditional realization
for (i in 1:5) {
  lines(grid,condReal(grid,xobs,yobs,sigma2,l),col="red")
}


#length evaluated by Monte Carlo
nmc = 200
vLength = seq(from=0,to=0,length=nmc)

for (i in 1:nmc) {
  ycond = condReal(grid,xobs,yobs,sigma2,l)
  estLengthCond = sum(
    sqrt(
      (grid[2:length(grid)] - grid[1:length(grid)-1])^2
      + (ycond[2:length(grid)] - ycond[1:length(grid)-1])^2
    ) 
  )
  vLength[i] = estLengthCond
}

lengthEst = mean(vLength)
lengthQsupMC = quantile(vLength,probs=0.95)
lengthQinfMC = quantile(vLength,probs=0.05)
