
ParallelBootstrap_UVUS <- function(x, B, stat.type, cl, F, cc){
  
  d = length(x)
  n.sample=as.matrix(sapply(x,dim))[1,]
  n=n.sample[1]
  p=as.matrix(sapply(x,dim))[2,]
  
  registerDoParallel(cl)
  
  output <- foreach(i = 1:B) %dopar% {
    source("R/jdcov_basic_functions.R")
    x.p <- x
    stat.p = 0
    
    for(j in 1:d) {x.p[[j]] <- as.matrix(x[[j]][sample(1:n,n,replace=TRUE),])}
    stat.p <- F(x.p,cc)
    
    stat.p
  }  
  stopCluster(cl)
  return(output)
}



ParallelBootstrap_RankUV <- function(x, B, stat.type, cl, F, cc){
  
  d = length(x)
  n.sample=as.matrix(sapply(x,dim))[1,]
  n=n.sample[1]
  p=as.matrix(sapply(x,dim))[2,]
  
  registerDoParallel(cl)
  
  output <- foreach(i = 1:B) %dopar% {
    source("R/jdcov_basic_functions.R")
    x.p <- x
    stat.p = 0
    for(j in 1:d) { x.p[[j]] <- as.matrix(x[[j]][sample(1:n,n,replace=TRUE),])}
    x.pr <- rank_list(x.p)  ;  stat.p <- F(x.pr,cc)
    stat.p
  }
  stopCluster(cl)
  return(output)
}



jdcov.test_parallel_edited <- function(x, cc, B, stat.type, alpha){
  
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  
  nworkers <- detectCores()
  cl = makeCluster(nworkers)
  
  if(stat.type=="V") {F <- Jdcov.sq.V ; stat <- F(x,cc)
  stat.p = unlist(ParallelBootstrap_UVUS(x,B,stat.type,cl,F,cc))
  }
  if(stat.type=="U") {F <- Jdcov.sq.U ; stat <- F(x,cc)
  stat.p = unlist(ParallelBootstrap_UVUS(x,B,stat.type,cl,F,cc))
  }
  if(stat.type=="US"){F <- Jdcov.sq.US ; stat <- F(x,cc)
  stat.p = unlist(ParallelBootstrap_UVUS(x,B,stat.type,cl,F,cc))
  }
  
  if(stat.type=="Rank V"){ F <- Jdcov.sq.V ; x.r <- rank_list(x) ; stat <- F(x.r,cc)
    stat.p = unlist(ParallelBootstrap_RankUV(x,B,stat.type,cl,F,cc))
  } 
  
  if(stat.type=="Rank U"){F <- Jdcov.sq.U ; x.r <- rank_list(x) ; stat <- F(x.r,cc)
    stat.p = unlist(ParallelBootstrap_RankUV(x,B,stat.type,cl,F,cc))
  }
  
  crit <- as.numeric(quantile(stat.p, 1-alpha))
  pvalue <- (1 + sum(stat.p>=stat))/(B+1)
  
  result=list(statistic=stat, crit.value=crit, p.value= pvalue)
  return(result)
}


library(mvtnorm)
n=100; d=5
set.seed(10)
z=rmvnorm(n,mean=rep(0,d),sigma=diag(d)) ; x=z^3
X <- lapply(seq_len(ncol(x)), function(i) as.matrix(x[,i]))
source("R/jdcov-test.R")
source("R/jdcov_basic_functions.R")
t1<-proc.time()
jdcov.test_parallel(X, cc=1, B=100, stat.type = "US", alpha=0.05)
proc.time()-t1

t2<-proc.time()
jdcov.test(X, cc=1, B=100, stat.type = "US", alpha=0.05)
proc.time()-t2

t3<-proc.time()
jdcov.test_parallel_edited(X, cc=1, B=100, stat.type = "U", alpha=0.05)
proc.time()-t3


library(mvtnorm)
n=100 ; d=3; p=5
x=list() ; x[[1]]=x[[2]]=x[[3]]=matrix(0,n,p)
set.seed(1)
x[[1]]=rmvnorm(n,rep(0,p),diag(p)) 
set.seed(2)
x[[2]]=rmvnorm(n,rep(0,p),diag(p))
set.seed(3)
W=rexp(n,1/sqrt(2))
x[[3]][,1]=(sign(x[[1]][,1] * x[[2]][,1])) * W
set.seed(4)
x[[3]][,2:p]=rmvnorm(n,rep(0,(p-1)),diag(p-1))

t1<-proc.time()
jdcov.test_parallel(X, cc=1, B=100, stat.type = "US", alpha=0.05)
proc.time()-t1

t2<-proc.time()
jdcov.test(X, cc=1, B=100, stat.type = "US", alpha=0.05)
proc.time()-t2

t3<-proc.time()
jdcov.test_parallel_edited(X, cc=1, B=100, stat.type = "US", alpha=0.05)
proc.time()-t3




