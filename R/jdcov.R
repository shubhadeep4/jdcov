

#' Title  Joint distance covariance
#' 
#' computes the joint distance covariance among random vectors X_1, .. , X_d with 
#' dimensions p_1, .. , p_d, defined in Section 2.2 in Chakraborty and Zhang (2019).
#'
#' @param x a list of d (>=2) elements, the i^th element being a numeric n*p_i data matrix for the random vector 
#'    X_i, i=1,..,d. n denotes the sample size and p_i denotes the dimension of the i^th random vector X_i.
#' @param cc a choice for the tuning parameter c, default is 1.
#' @param stat.type  a character string specifying the type of the estimator of joint distance covariance 
#'                   to be computed. The available options  are : "V" (the V-statistic type estimator), 
#'                   "U" (the U-statistic type estimator), "US" (the scale invariant U-statistic type estimator),
#'                  "Rank V" (the rank based V-statistic type estimator) and "Rank U" (the rank based U-statistic
#'                   type estimator). Default is "U".
#'
#' @return  the value of the squared sample joint distance covariance depending on the type of estimator 
#'          specified
#' @export
#'
#' @examples
#' 
#'  (X_1, .. , X_d) are generated as follows : X=Z^3, where Z~N(0,diag(d)).
#' library(mvtnorm)
#' n=100; d=5
#' set.seed(10)
#' z=rmvnorm(n,mean=rep(0,d),sigma=diag(d)) ; x=z^3
#' X <- lapply(seq_len(ncol(x)), function(i) as.matrix(x[,i]))
#' jdcov(X, cc=1, stat.type = "U")



#'  X, Y and Z are p-dimensional random vectors, where X, Y are i.i.d N(0, diag(p)), 
#' Z_1 = sign(X_1 * Y_1) * W and Z_{2:p} ~ N(0, diag(p-1)), W ~ exponential(mean=sqrt(2)).
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
#' jdcov(x, cc=1, stat.type = "US")

 
jdcov <- function(x, cc=1, stat.type="U"){
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
  return(stat)
}

