
## JdCov : computes the joint distance covariance among random vectors X_1, .. , X_d with 
## dimensions p_1, .. , p_d, defined in Section 2.2 in Chakraborty and Zhang (2019).

## Inputs :
## x : A list of d elements, the i^th element being the numeric n*p_i data matrix for the random vector X_i,
##     i=1,..,d.
## cc : A choice for the tuning parameter c (Section 2.2, Chakraborty and Zhang (2019)), default is 1.
## stat.type : the type of the estimator of JdCov to be used. Must be one of "V" (V-statistic type), "U" 
## (U-statistic type), "US" (scale invariant U-statistic type), "Rank V" (rank based V-statistic type) or
## "Rank U" (rank based U-statistic type), default is "U".

#' Title
#'
#' @param x 
#' @param cc 
#' @param stat.type 
#'
#' @return
#' @export
#'
#' @examples
JdCov <- function(x, cc=1, stat.type="U"){
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  
  if(stat.type=="V") stat <- Jdcov.sq.V(x,cc)
  if(stat.type=="U") stat <- Jdcov.sq.U(x,cc)
  if(stat.type=="US") stat <- Jdcov.sq.US(x,cc)
  if(stat.type=="Rank U") {
    x.r <- x
    for(j in 1:d) {
      for(l in 1:ncol(x[[j]])){
        f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
      }
    }
    stat <- Jdcov.sq.U(x.r,cc)
  }
  if(stat.type=="Rank V") {
    x.r <- x
    for(j in 1:d) {
      for(l in 1:ncol(x[[j]])){
        f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
      }
    }
    stat <- Jdcov.sq.V(x.r,cc)
  }
  result=list(statistic=stat.type, value=stat)
  return(result)
}

