kmeans_edf = function (X, maxk, nstart = 10, ngrid = 30, bw = 5) 
{
  X <- t(t(X) - colMeans(X))
  sds <- sqrt(diag(cov(X)))
  if (min(sds) < 1e-10) 
    X <- X[, -which(sds < 1e-10)]
  n <- nrow(X)
  d <- ncol(X)
  clusters <- matrix(1, maxk, n)
  sols <- list()
  for(k in 2:(maxk+1)) sols[[k]] <- kmeans(X, k, nstart = nstart)
  sols[[1]]$cluster <- numeric(n) + 1
  sols[[1]]$centers <- matrix(colMeans(X), nrow = 1)
  SS <- c(sols[[2]]$totss, unlist(lapply(sols, function(l) l$tot.withinss)))
  c0 <- elbow(SS)
  if(ngrid==1) sigs <- c(sqrt(SS[c0]/nrow(X)/ncol(X)))
  else sigs <- sqrt(seq(SS[c0]/2.5, 2.5*SS[c0], length = ngrid)/nrow(X)/ncol(X))
  EDFS <- matrix(0, maxk, ngrid)
  EDFS[1,] <- d
  for(k in 2:maxk){
    EDFS[k,] <- edfkmeans(X, sols[[k]]$centers, sapply(1:(k+1), function(jj) sols[[k+1]]$cluster==jj)%*%sols[[k+1]]$centers, sols[[k]]$cluster - 1, sols[[k]]$size, sols[[k]]$tot.withinss, k, n, d, 
                                                   sigs, ngrid)
    clusters[k,] <- sols[[k]]$cluster
  }
  edfs <- EDFS
  for(i in 1:ngrid) EDFS[,i] <- smoothed(EDFS[,i], bw)
  BICS <- SS[1:maxk]%*%t(1/sigs^2) + EDFS*log(n*d)
  ks <- apply(BICS, 2, fmin)
  k <- ifelse(sum(ks==1)>=(ngrid*2/3), 1, ifelse(sum(ks==maxk)>=(ngrid*2/3), maxk, which.max(sapply(2:(maxk-1), function(i) sum(ks==i)))+1))
  list(ss = SS, bic = BICS, cluster = clusters, sigs = sigs, k = k, edfs = edfs)
}






elbow <- function(yvals){
  d <- length(yvals)
  angles <- sapply(1:(d - 1), function(i) {
    a1 <- atan((i - 1)/(d - 1) * (yvals[1] - yvals[d])/(yvals[1] - 
                                                          yvals[i]))
    a2 <- atan((yvals[i] - yvals[d])/(yvals[1] - yvals[d]) * 
                 (d - 1)/(d - i))
    return(a1 + a2)
  })
  which.min(angles)
}

fmin = function(x){
  if(x[1]<min(x[-1])) 1
  else{
    mins = which(sapply(2:(length(x)-1), function(i) x[i]<=(min(x[i+1], x[i-1]))))
    if(length(mins)>0) min(mins) + 1
    else length(x)
  } 
}



smoothed <- function(x, bw){
  KernSmooth::locpoly(1:length(x), x, range.x = c(1, length(x)), gridsize = length(x), bandwidth = bw)$y
 }
