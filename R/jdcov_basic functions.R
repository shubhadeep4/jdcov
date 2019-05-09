## --------------------------- Basic functions------------------------------------------ ##


## v.center : computes the matrix U^hat for a random vector, defined in Section 3.1 of Chakraborty and 
## Zhang (2019).

## Input : 
## x : a numeric n*p data matrix (n=sample size, p= dimension of the random vector).

v.center <- function(x){
  n=nrow(x)
  A=as.matrix(dist(x))
  R=rowSums(A)
  C=colSums(A)
  T=sum(A)                             
  r=matrix(rep(R,n),n,n)/n
  c=t(matrix(rep(C,n),n,n))/n
  t=matrix(T/n^2,n,n)
  UA=-(A-r-c+t)
  return(UA)
}


## Jdcov.sq.V : computes the V-statistic type estimator of JdCov for d random vectors X_1, .. , X_d 
## with dimensions p_1, .. , p_d, given in equation (11) in Chakraborty and Zhang (2019).

## Inputs :
## x : A list of d elements, the i^th element being the numeric n*p_i data matrix for the random vector X_i,
##     i=1,..,d.
## cc : A choice for the tuning parameter c (Section 2.2, Chakraborty and Zhang (2019)), default is 1.

Jdcov.sq.V <- function(x,cc=1){
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  
  A <- v.center(x[[1]]) + cc 
  for(i in 2:d) A <- A*( v.center(x[[i]]) + cc )
  return( sum(A)/n^2 - cc^d )
}



## u.center : computes the U-centered matrix U^tilde for a random vector, defined in Section 3.2 in 
## Chakraborty and Zhang (2019).

## Input : 
## x : a numeric n*p data matrix (n=sample size, p= dimension of the random vector).

u.center <- function(x){
  n=nrow(x)
  A=as.matrix(dist(x))
  R=rowSums(A)
  C=colSums(A)
  T=sum(A)
  r=matrix(rep(R,n),n,n)/(n-2)
  c=t(matrix(rep(C,n),n,n))/(n-2)
  t=matrix(T/(n-1)/(n-2),n,n)
  UA=-(A-r-c+t)
  diag(UA)=0
  return(UA)
}



## Jdcov.sq.U : computes the U-statistic type estimator of JdCov for d random vectors X_1, .. , X_d 
## with dimensions p_1, .. , p_d, given in Section 3.2 in Chakraborty and Zhang (2019).

## Inputs :
## x : A list of d elements, the i^th element being the numeric n*p_i data matrix for the random vector X_i,
##     i=1,..,d.
## cc : A choice for the tuning parameter c (Section 2.2, Chakraborty and Zhang (2019)), default is 1.

Jdcov.sq.U <- function(x,cc=1){
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  A <- u.center(x[[1]]) + cc 
  for(i in 2:d) A <- A*( u.center(x[[i]]) + cc )
  return( sum(A)/n/(n-3)-cc^d*n/(n-3) )
}



## dCov.sq.U : given a random vector X, it computes the U-centered estimator of dCov^2(X,X), given in
## Section 3.2, Chakraborty and Zhang (2019)

## Input : 
## x : a numeric n*p data matrix (n=sample size, p= dimension of the random vector).

dCov.sq.U <- function(x,n) {
  A <- u.center(x)
  return( sum(A^2)/n/(n-3) )
}



## Jdcov.sq.US : computes the scale-invariant U-statistic type estimator of JdCov for d random vectors 
## X_1, .. , X_d with dimensions p_1, .. , p_d, defined in Section 3.2 in Chakraborty and Zhang (2019).

## Inputs :
## x : A list of d elements, the i^th element being the numeric n*p_i data matrix for the random vector X_i,
##     i=1,..,d.
## cc : A choice for the tuning parameter c (Section 2.2, Chakraborty and Zhang (2019)), default is 1.


Jdcov.sq.US <- function(x,cc=1){
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  A <- u.center(x[[1]])/sqrt( dCov.sq.U(x[[1]],n) ) + cc         
  for(i in 2:d) A <- A*( u.center(x[[i]])/sqrt( dCov.sq.U(x[[i]],n) ) + cc )
  return( sum(A)/n/(n-3)-cc^d*n/(n-3) )
}