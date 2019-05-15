#' Test for joint independence based on joint distance covariance
#' 
#' Performs a bootstrap based test for joint independence among d(>=2) random vectors X_1, .. , X_d  of (arbitrary) dimensions p_1, .. , p_d, described in Section 4 in Chakraborty and Zhang (2019). The test is based on the joint distance covariance (JdCov) among X_1, .. , X_d, which is shown to completely characterize joint independence among the random vectors. The null hypothesis (H_0) is that all the random vectors are jointly independent.
#' @param x a list of d (>=2) elements, the i^th element being a numeric n*p_i data matrix for the random vector X_i,  i = 1,..,d.  n denotes the sample size and p_i  denotes the dimension of the i^th random vector X_i.
#' @param cc a numeric value specifying the choice of the tuning parameter c, default is 1.
#' @param B an integer value specifying the number of bootstrap replicates to be considered, default is 100.
#' @param stat.type a character string specifying the type of the estimator of joint distance covariance to be computed. The available options are : "V" (the V-statistic type estimator), "U" (the U-statistic type estimator), "US" (the scale invariant U-statistic type estimator), "Rank V" (the rank based V-statistic type estimator) and "Rank U" (the rank based U-statistic type estimator). Default is "U".
#' @param alpha a numeric value specifying the level of the test, default is 0.05.
#'
#' @return A list containing the following components :
#' 
#' statistic       :       the observed value of the test statistic
#' 
#' crit.value      :       the critical value of the bootstrap-based hypothesis test. The null hypothesis (H_0: joint independence) is rejected if the observed value of the statistic is greater than crit.value
#' 
#' p.value         :       the p-value of the bootstrap-based hypothesis test. The null hypothesis (H_0: joint independence) is rejected if the p.value is smaller than alpha.
#'
#' @references
#' Chakraborty, S. and Zhang, X. (2019). Distance Metrics for Measuring Joint Dependence with Application to Causal Inference, Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1513364.
#'
#' @examples
#' ## (X_1, .. , X_d) are generated as follows : X=Z^3, where
#' ## Z~N(0,diag(d)). 
#' 
#' library(mvtnorm)
#' n=100; d=5
#' set.seed(10)
#' z=rmvnorm(n,mean=rep(0,d),sigma=diag(d)) ; x=z^3
#' X <- lapply(seq_len(ncol(x)), function(i) as.matrix(x[,i]))
#' jdcov.test(X, cc=1, B=100, stat.type = "U", alpha=0.05) 
#' 
#' 
#' ## X, Y and Z are p-dimensional random vectors, where X, Y  
#' ## are i.i.d N(0, diag(p)), Z_1 = sign(X_1 * Y_1) * W and 
#' ## Z_{2:p} ~ N(0, diag(p-1)), W ~ exponential(mean=sqrt(2)). 
#' 
#' library(mvtnorm)
#' n=100 ; d=3; p=5
#' x=list() ; x[[1]]=x[[2]]=x[[3]]=matrix(0,n,p)
#' set.seed(1)
#' x[[1]]=rmvnorm(n,rep(0,p),diag(p)) 
#' set.seed(2)
#' x[[2]]=rmvnorm(n,rep(0,p),diag(p))
#' set.seed(3)
#' W=rexp(n,1/sqrt(2))
#' x[[3]][,1]=(sign(x[[1]][,1] * x[[2]][,1])) * W
#' set.seed(4)
#' x[[3]][,2:p]=rmvnorm(n,rep(0,(p-1)),diag(p-1))
#' jdcov.test(x, cc=1, B=100, stat.type = "U", alpha=0.05)

 
jdcov.test <- function(x, cc=1, B=100, stat.type="U", alpha=0.05){
  
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
  
  result=list(statistic=stat, crit.value=crit, p.value= pvalue)
  return(result)
}

