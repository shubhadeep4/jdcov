nworkers <- detectCores()
cl = makeCluster(nworkers)
 
ParallelBootstrap <- function(x, B, stat.type, cl, F, cc){
  d = length(x)
  n.sample=as.matrix(sapply(x,dim))[1,]
  n=n.sample[1]
  p=as.matrix(sapply(x,dim))[2,]
  
  registerDoParallel(cl)
  
  output <- foreach(i = 1:B) %dopar% {
    source("R/jdcov_basic_functions.R")
    x.p <- x
    stat.p = 0
    for(j in 1:d) x.p[[j]] <- as.matrix(x[[j]][sample(1:n,n,replace=TRUE),])
    
    if(stat.type=="V") stat.p <- F(x.p,cc)
    if(stat.type=="U") stat.p <- F(x.p,cc)
    if(stat.type=="US") stat.p <- F(x.p,cc)
    
    if(stat.type=="Rank V"){
      for(j in 1:d) {
        x.pr <-  x
        for(l in 1:ncol(x[[j]])){
          f.cdf <- ecdf(x.p[[j]][,l]); x.pr[[j]][,l] <- f.cdf(x.p[[j]][,l])
        }
      }
      stat.p <- F(x.pr,cc)
    }
    
    if(stat.type=="Rank U"){
      for(j in 1:d) {
        x.pr <-  x
        for(l in 1:ncol(x[[j]])){
          f.cdf <- ecdf(x.p[[j]][,l]); x.pr[[j]][,l] <- f.cdf(x.p[[j]][,l])
        }
      }
      stat.p <- F(x.pr,cc)
    }
    stat.p
  }
  stopCluster(cl)
  return(output)
}

jdcov.test_parallel <- function(x, cc=1, B=100, stat.type="U", alpha=0.05){
  
  n.sample=as.matrix(sapply(x,dim))[1,]
  if(length(unique(n.sample))!=1) stop("Unequal Sample Sizes")
  n=n.sample[1]
  d=length(x)
  p=as.matrix(sapply(x,dim))[2,]
  
  if(stat.type=="V") {F <- Jdcov.sq.V ; stat <- F(x,cc)
  stat.p = unlist(ParallelBootstrap(x,B,stat.type,cl,F,cc))
  }
  if(stat.type=="U") {F <- Jdcov.sq.U ; stat <- F(x,cc)
  stat.p = unlist(ParallelBootstrap(x,B,stat.type,cl,F,cc))
  }
  if(stat.type=="US"){F <- Jdcov.sq.US ; stat <- F(x,cc)
  stat.p = unlist(ParallelBootstrap(x,B,stat.type,cl,F,cc))
  }
  
  if(stat.type=="Rank V"){
    F <- Jdcov.sq.V
    x.r <- x
    for(j in 1:d) {
      for(l in 1:ncol(x[[j]])){
        f.cdf <- ecdf(x[[j]][,l]); x.r[[j]][,l] <- f.cdf(x[[j]][,l])
      }
    }
    stat <- F(x.r,cc)
    stat.p = unlist(ParallelBootstrap(x,B,stat.type,cl,F,cc))
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
    stat.p = unlist(ParallelBootstrap(x,B,stat.type,cl,F,cc))
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
jdcov.test_parallel(X, cc=1, B=100, stat.type = "U", alpha=0.05)
source("R/jdcov-test.R")
jdcov.test_parallel(X, cc=1, B=100, stat.type = "U", alpha=0.05)
proc.time()

jdcov.test(X, cc=1, B=100, stat.type = "U", alpha=0.05)
proc.time()


