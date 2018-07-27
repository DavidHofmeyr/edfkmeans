BICkmeans <- function(X, maxk, sig, nstart = 10){
  SS <- c()
  dfs <- c()
  X <- t(t(X)-colMeans(X))
  sds <- sqrt(diag(cov(X)))
  if(min(sds)<1e-10) X <- X[,-which(sds<1e-10)]
  n <- nrow(X)
  for(k in 2:maxk){
    km <- kmeans(X, k, nstart = nstart)
    SS <- c(SS, km$tot.withinss)
    dfs <- c(dfs, edfkmeans(X, km$centers, km$cluster-1, km$size, km$tot.withinss, k, nrow(X), ncol(X), c(sig), 1))
  }
  SS/sig^2+log(n)*dfs
}

get_k <- function(X, maxk, sig, nstart = 10){
  bic <- BICkmeans(X, maxk, sig, nstart)
  which.min(bic)+1
}

get_k_elbow <- function(X, maxk, nstart = 10){
  bic <- BICkmeans(X, maxk, get_sig_elbow(X), nstart)
  which.min(bic)+1
}

get_k_grid <- function(X, maxk, nstart = 10, ngrid = 15){
  X <- t(t(X)-colMeans(X))
  sds <- sqrt(diag(cov(X)))
  if(min(sds)<1e-10) X <- X[,-which(sds<1e-10)]
  n <- nrow(X)
  d <- ncol(X)
  BICS <- matrix(0, maxk, ngrid)
  for(k in 2:maxk){
    km <- kmeans(X, k, nstart = nstart)
    BICS[k,] <- km$tot.withinss/(seq(.31, 1.75, length = ngrid)^2*km$totss/n/d) + log(n)*edfkmeans_cpp(X, km$centers, km$cluster-1, km$size, km$tot.withinss, k, n, d, seq(.31, 1.75, length = ngrid)*sqrt(km$totss/n/d), ngrid)
  }
  BICS[1,] <- km$totss/(seq(.31, 1.75, length = ngrid)^2*km$totss/n/d) + log(n)*d
  mins <- unlist(apply(BICS,2,function(x){
    ret <- c()
    for(i in 2:(length(x)-1)) if((x[i]<min(x[i-1],x[i+1])) && x[i]==min(x[-c(1,length(x))])) ret <- c(ret, i)
    ret
  }))
  which.max(sapply(1:maxk, function(k) sum(mins==k)))
}


get_sig_elbow <- function(X){
  yvals <- eigen(cov(X))$values
  d <- ncol(X)
  angles <- sapply(1:(d-1), function(i){
    a1 <- atan((i-1)/(d-1)*(yvals[1]-yvals[d])/(yvals[1]-yvals[i]))
    a2 <- atan((yvals[i]-yvals[d])/(yvals[1]-yvals[d])*(d-1)/(d-i))
    return(a1+a2)
  })
  sqrt(yvals[1+which.min(angles)])
}
