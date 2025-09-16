#Author: Francois BACHOC
#
Balpha = function(q,N,alpha,I=1000) {
  #compute the upper bound for quantile of max of absolute value of inner products between standard Gaussian vector and unit vectors
  #cf K4 in Bachoc Leep Poetscher 2014 and Balpha in Bachoc Preinerstorfer Steinberger 2016 and in Bachoc Blanchard Neuvial 2017
  #q...........: dimension of the Gaussian vector
  #N...........: number of unit vectors
  #alpha.......: confidence level
  #I...........: numerical precision
  #return......: upper bound, numerically approximated
  #
  vC = sqrt(qbeta(p=seq(from=0,to=1/N,length=I),shape1=1/2,shape2=(q-1)/2, lower.tail=FALSE)) #vector of quantiles of Beta distribution
  fconfidence = function(K){mean(pchisq(q=(K/vC)^2,df=q,lower.tail=FALSE))-(1-alpha)} #MC evaluation of confidence level for a constant K
  uniroot(fconfidence,interval=c(1,100))$root
}

