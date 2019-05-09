## high.order.dCov : computes the d^th order distance covariance among random vectors 
## X_1, .. , X_d with dimensions p_1, .. , p_d, defined in Section 2.1 in Chakraborty and Zhang (2019).

## Inputs :
## x : A list of d elements, the i^th element being the numeric n*p_i data matrix for the random vector X_i,
##     i=1,..,d.
## type : the type of the estimator of dCov to be used. Must be one of "V" (V-statistic type) or "U" 
## (U-statistic type), default is "U".

#' Title
#'
#' @param x 
#' @param type 
#'
#' @return
#' @export
#'
#' @examples
high.order.dCov <- function(x, type="U"){
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  if(type=="U"){
    A <- u.center(x[[1]])
    for(i in 2:d) A <- A* u.center(x[[i]])
    return( sum(A)/n/(n-3))
  }
  if(type=="V"){
    A <- v.center(x[[1]])
    for(i in 2:d) A <- A* v.center(x[[i]])
    return( sum(A)/n/n)
  }
}
