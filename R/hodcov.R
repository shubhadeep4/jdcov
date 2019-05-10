#' High order distance covariance
#'
#' @param x  a list of d (>=2) elements, the i^th element being a numeric n*p_i data matrix for the random vector X_i, i=1,..,d. n denotes the sample size and p_i denotes the dimension of the i^th random vector X_i.
#' @param type a character, either "U" (computes the U-statistic type estimator of the squared d^th order dCov) or "V" (computes the V-statistic type estimator of the squared d^th order dCov)
#'
#' @return the value of the U-statistic type or the V-statistic type estimator of the squared d^th order dCov, based on whether type is "U" or "V".
#'
#' @examples 
#' ## (X_1, .. , X_d) are generated as follows : X=Z^3, where Z~N(0,diag(d)).
#' library(mvtnorm)
#' n=100; d=5
#' set.seed(10)
#' z=rmvnorm(n,mean=rep(0,d),sigma=diag(d)) ; x=z^3
#' X <- lapply(seq_len(ncol(x)), function(i) as.matrix(x[,i]))
#' hodcov(X, type="U")
#' ## X, Y and Z are p-dimensional random vectors, where X, Y are i.i.d N(0, diag(p)), 
#' ## Z_1 = sign(X_1 * Y_1) * W and Z_{2:p} ~ N(0, diag(p-1)), W ~ exponential(mean=sqrt(2)).
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
#' hodcov(x, type="V")

hodcov <- function(x, type="U"){
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
