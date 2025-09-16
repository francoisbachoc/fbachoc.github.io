


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
sigma2=1
l=0.5

#failure of the gaussien



# Condition number
vn = c(10,100,500,1000)
vcondexp = vn
vcond3s2 = vn
vcond5s2 = vn
vcondgauss = vn
for (i in 1:4) {
  x = seq(from=0,to=1,length=vn[i])
  R = matern1s2(x,x,sigma2,l)
  vcondexp[i] = rcond(R)
  R = matern3s2(x,x,sigma2,l)
  vcond3s2[i] = rcond(R)
  R = matern5s2(x,x,sigma2,l)
  vcond5s2[i] = rcond(R)
  R = gaussien(x,x,sigma2,l)
  vcondgauss[i] = rcond(R)
}

plot(vn,log10(vcondexp),ylim=c(-30,0))
points(vn,log10(vcond3s2),col='blue')
points(vn,log10(vcond5s2),col='green')
points(vn,log10(vcondgauss),col='red')

#comparaison regularite
nugget = 10^(-8)
n = 1000
x = seq(from=0,to=1,length=n)
R = matern1s2(x,x,sigma2,l) + nugget*diag(n)
CR = chol(R)
y = t(CR)%*%matrix(nrow=n,ncol=1,data=rnorm(n))
plot(x,y,ylim=c(-2,2),type='l')
R = matern3s2(x,x,sigma2,l) + nugget*diag(n)
CR = chol(R)
y = t(CR)%*%matrix(nrow=n,ncol=1,data=rnorm(n))
points(x,y,type='l',col='blue')
R = gaussien(x,x,sigma2,l) + nugget*diag(n)
CR = chol(R)
y = t(CR)%*%matrix(nrow=n,ncol=1,data=rnorm(n))
points(x,y,type='l',col='green')



#comparaison ell
nugget = 10^(-8)
n = 1000
x = seq(from=0,to=1,length=n)
R = matern3s2(x,x,sigma2,0.2) + nugget*diag(n)
CR = chol(R)
y = t(CR)%*%matrix(nrow=n,ncol=1,data=rnorm(n))
plot(x,y,ylim=c(-2,2),type='l')
R = matern3s2(x,x,sigma2,0.5) + nugget*diag(n)
CR = chol(R)
y = t(CR)%*%matrix(nrow=n,ncol=1,data=rnorm(n))
points(x,y,type='l',col='blue')
R = matern3s2(x,x,sigma2,1) + nugget*diag(n)
CR = chol(R)
y = t(CR)%*%matrix(nrow=n,ncol=1,data=rnorm(n))
points(x,y,type='l',col='green')



