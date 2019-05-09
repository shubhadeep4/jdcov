## JdCov.test : performs a bootstrap based test for joint independence of d random vectors 
## X_1, .. , X_d with dimensions p_1, .. , p_d, given in Section 4 in Chakraborty and Zhang (2019).

## Inputs :
## x : A list of d elements, the i^th element being the numeric n*p_i data matrix for the random vector X_i,
##     i=1,..,d.
## cc : A choice for the tuning parameter c (Section 2.2, Chakraborty and Zhang (2019)), default is 1.
## B : number of bootstrap resamples, default is 100.
## stat.type : the type of the estimator of JdCov to be used. Must be one of "V" (V-statistic type), "U" 
## (U-statistic type), "US" (scale invariant U-statistic type), "Rank V" (rank based V-statistic type) or
## "Rank U" (rank based U-statistic type), default is "U".
##  alpha : level of the test, default is 0.05.


JdCov.test <- function(x, cc=1, B=100, stat.type="U", alpha=0.05){
  
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  
  if(stat.type=="V") {F <- Jdcov.sq.V ; stat <- F(x,cc)}
  if(stat.type=="U") {F <- Jdcov.sq.U ; stat <- F(x,cc)}
  if(stat.type=="US"){F <- Jdcov.sq.US ; stat <- F(x,cc)}
  
  if(stat.type=="Rank V"){
    F <- Jdcov.sq.V
    x.r <- x
    for(j in 1:d) {
      for(l in 1:ncol(x[[j]])){
        f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
      }
    }
    stat <- F(x.r,cc)
  } 
  
  if(stat.type=="Rank U"){
    F <- Jdcov.sq.U
    x.r <- x
    for(j in 1:d) {
      for(l in 1:ncol(x[[j]])){
        f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
      }
    }
    stat <- F(x.r,cc)
  }
  
  stat.p <- rep(0, B)
  
  for(i in 1:B)
  {
    x.p <- x
    for(j in 1:d) x.p[[j]] <- as.matrix(x[[j]][sample(1:n,n,replace=TRUE),])
    
    if(stat.type=="V") stat.p[i] <- F(x.p,cc)
    if(stat.type=="U") stat.p[i] <- F(x.p,cc)
    if(stat.type=="US") stat.p[i] <- F(x.p,cc)
    
    if(stat.type=="Rank V"){
      for(j in 1:d) {
        x.pr <-  x
        for(l in 1:ncol(x[[j]])){
          f.cdf <- ecdf(x.p[[j]][,l]); x.pr[[j]][,l] <- f.cdf(x.p[[j]][,l])
        }
      }
      stat.p[i] <- F(x.pr,cc)
    }
    
    if(stat.type=="Rank U"){
      for(j in 1:d) {
        x.pr <-  x
        for(l in 1:ncol(x[[j]])){
          f.cdf <- ecdf(x.p[[j]][,l]); x.pr[[j]][,l] <- f.cdf(x.p[[j]][,l])
        }
      }
      stat.p[i] <- F(x.pr,cc)
    }
    
  }
  
  crit <- as.numeric(quantile(stat.p, 1-alpha))
  pvalue <- (1 + sum(stat.p>=stat))/(B+1)
  
  result=list(statistic = stat.type, value=stat, crit.value=crit, p.value= pvalue)
  return(result)
}

