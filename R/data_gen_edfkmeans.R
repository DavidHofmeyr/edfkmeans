datagen_edfkmeans <- function (n, d, k, var = 0, v = 0, scale = 0, noisedim = 0, props = 0) 
{
  X0 <- matrix(runif(k * d^2), k * d, d) * 0.8^var * sqrt(10)/sqrt(d)
  MU <- kmeans(X0, k, iter.max = 1)$centers
  X <- matrix(0, n, d)
  scl <- runif(k)^scale
  scl <- scl/mean(scl)/2
  S <- matrix(rnorm(d^2), ncol = d)
  SS <- t(S) %*% S
  V <- t(eigen(SS)$vectors[, 1:(d - noisedim)])
  p <- runif(k)^props
  p <- p/sum(p)
  ns <- rowSums(rmultinom(n, 1, p))
  while (sum(ns == 0) > 0) ns <- rowSums(rmultinom(n, 1, p))
  done <- 0
  labels <- c()
  for (i in 1:k) {
    X[(done + 1):(done + ns[i]), ] <- t(matrix(rnorm(ns[i] * 
                                                       (d - noisedim)), nrow = d - noisedim) * runif(d - 
                                                                                                       noisedim)^v/3 * (v + 1) * scl[i] + MU[i, 1:(d - noisedim)]) %*% 
      V
    done <- done + ns[i]
    labels <- c(labels, rep(i, ns[i]))
  }
  if (noisedim > 0) {
    X <- X + matrix(rnorm(nrow(X) * noisedim)/2, nrow = nrow(X)) %*% 
      t(eigen(SS)$vectors[, (d - noisedim + 1):d])
  }
  list(x = X, c = labels)
}



smoothed <- function(x){
  bw <- .5
  sm <- KernSmooth::locpoly(1:length(x), x, bandwidth = bw, gridsize = length(x))$y
  mins <- sapply(2:(length(x)-1), function(i) sm[i] <= min(sm[(i-1):(i+1)]))
  while(sum(mins) > 1){
    bw <- bw + .1
    sm <- KernSmooth::locpoly(1:length(x), x, bandwidth = bw, gridsize = length(x))$y
    mins <- sapply(2:(length(x)-1), function(i) sm[i] <= min(sm[(i-1):(i+1)]))
  }
  sm
}




kmeans_edf = function (X, kmax, nstart = 10, ntest = 1, extra = 1) 
{
  X <- t(t(X) - colMeans(X))
  sds <- sqrt(diag(cov(X)))
  if (min(sds) < 1e-10) 
    X <- X[, -which(sds < 1e-10)]
  n <- nrow(X)
  d <- ncol(X)
  clusters <- matrix(1, kmax, n)
  sols <- list()
  for(k in 2:kmax) sols[[k]] <- kmeans(X, k, nstart = nstart)
  sols[[1]]$cluster <- numeric(n) + 1
  sols[[1]]$centers <- matrix(colMeans(X), nrow = 1)
  SS <- c(sols[[2]]$totss, unlist(lapply(sols, function(l) l$tot.withinss)))
  MM <- matrix(0, n, d*ntest)
  sigs <- numeric(ntest)
  for(k in 1:ntest){
    km <- kmeans(X, kmax + extra, nstart = nstart)
    MM[,((k-1)*d+1):(k*d)] <- sapply(1:(kmax+extra), function(jj) km$cluster==jj)%*%km$centers
    sigs[k] <- sqrt(km$tot.withinss/n/d)
  }
  EDFS <- matrix(0, kmax, ntest)
  EDFS[1,] <- d
  for(k in 2:kmax){
    EDFS[k,] <- edfkmeans_all(X, sols[[k]]$centers, MM, sols[[k]]$cluster - 1, sols[[k]]$size, sols[[k]]$tot.withinss, k, n, d, 
                              sigs, ntest)
    clusters[k,] <- sols[[k]]$cluster
  }
  edfsraw <- rowMeans(EDFS)
  edfs <- smoothed(edfsraw)
  edfs[1] <- d
  bic <- n*d*log(SS) + log(n*d)*edfs
  k <- fmin(bic)
  list(ss = SS, bic = bic, cluster = clusters, edfs = edfs, edfs0 = edfsraw, k = k)
}
                                 
fmin = function(x){
  if(x[1]<min(x[-1])) 1
  else{
    mins = which(sapply(2:(length(x)-1), function(i) x[i]<=(min(x[i+1], x[i-1]))))
    if(length(mins)>0) min(mins) + 1
    else length(x)
  } 
}



