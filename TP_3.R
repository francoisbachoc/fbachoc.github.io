


#Function covariance matrix
matern1s2 = function(x1,x2,y1,y2,sigma2,l) {
  M1 = outer(x1,y1,"-")
  M1 = exp(-(abs(M1)/l))
  M2 = outer(x2,y2,"-")
  M2 = exp(-(abs(M2)/l))
  sigma2*M1*M2
}

#The conditional mean function
pred = function(x,y,xobs,yobs,fobs,sigma2,l)  {
  fobs = matrix(nrow=length(fobs),ncol=1,data=fobs)
  Kobs=matern1s2(xobs,yobs,xobs,yobs,sigma2,l) + 1e-8*diag(length(fobs))
  KobsPred=matern1s2(xobs,yobs,x,y,sigma2,l)
  pred = t(KobsPred)%*%solve(Kobs)%*%fobs
  return(pred)
}

#covariance parameters we use
sigma2=0.2
l=0.4

#The unknown function
f = function(x,y) {
  (x*y)
}

##############################################
#Les vrais indices
##############################################
I1 = 3/7
I12 = 4/7
I2 = 3/7
TI1=4/7


##############################################
#Evaluation of indices by Monte Carlo
##############################################
nmc=100000
xA=runif(nmc)
yA=runif(nmc)
xB=runif(nmc)
yB=runif(nmc)

#I1
fA = f(xA,yA)
fB = f(xA,yB)
I1hatMC = (mean(fA*fB) - mean(fA)^2) / ( mean( fA^2)-mean(fA)^2 )

#I2
fA = f(xA,yA)
fB = f(xB,yA)
I2hatMC = (mean(fA*fB) - mean(fA)^2) / ( mean( fA^2)-mean(fA)^2 )

#TI1
fA = f(xA,yA)
fB = f(xB,yA)
TI1hatMC = (1/2)*mean( (fA-fB)^2  ) / ( mean( fA^2)-mean(fA)^2 )

#TI2
fA = f(xA,yA)
fB = f(xA,yB)
TI2hatMC = (1/2)*mean( (fA-fB)^2  ) / ( mean( fA^2)-mean(fA)^2 )

##############################################
#Evaluation of indices by Kriging
##############################################
#The observation vector for kriging
xobs = seq(from=0,to=0,length=4)
yobs = seq(from=0,to=0,length=4)
cpt=0
for (x in seq(from=0,to=1,length=2)) {
  for (y in seq(from=0,to=1,length=2)) {
    cpt=cpt+1
    xobs[cpt] = x
    yobs[cpt] = y
  }
}
fobs = f(xobs,yobs)

#plot prediction errors
x=runif(100)
y=runif(100)
plot(f(x,y),pred(x,y,xobs,yobs,fobs,sigma2,l))

#for Monte Carlo
nmc=100000
xA=runif(nmc)
yA=runif(nmc)
xB=runif(nmc)
yB=runif(nmc)

#I1
fA = pred(xA,yA,xobs,yobs,fobs,sigma2,l)
fB = pred(xA,yB,xobs,yobs,fobs,sigma2,l)
I1krig = (mean(fA*fB) - mean(fA)^2) / ( mean( fA^2)-mean(fA)^2 )

#I2
fA = pred(xA,yA,xobs,yobs,fobs,sigma2,l)
fB = pred(xB,yA,xobs,yobs,fobs,sigma2,l)
I2krig = (mean(fA*fB) - mean(fA)^2) / ( mean( fA^2)-mean(fA)^2 )

#TI1
fA = pred(xA,yA,xobs,yobs,fobs,sigma2,l)
fB = pred(xB,yA,xobs,yobs,fobs,sigma2,l)
TI1krig = (1/2)*mean( (fA-fB)^2  ) / ( mean( fA^2)-mean(fA)^2 )

#TI2
fA = pred(xA,yA,xobs,yobs,fobs,sigma2,l)
fB = pred(xA,yB,xobs,yobs,fobs,sigma2,l)
TI2krig = (1/2)*mean( (fA-fB)^2  ) / ( mean( fA^2)-mean(fA)^2 )


