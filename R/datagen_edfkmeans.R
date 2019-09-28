datagen_edfkmeans <- function(n, d, k, var = 1.5, v = 0, scale = 0, noisedim = 0, props = 0){
  X0 <- matrix(rnorm(n*d), n, d)*.8^var*sqrt(10)/sqrt(d)
  MU <- kmeans(X0, k, iter.max = 1)$centers
  X <- matrix(0, n, d)
  
  scl <- runif(k)^scale
  scl <- scl/mean(scl)/2
  
  S <- matrix(rnorm(d^2), ncol = d)
  SS <- t(S)%*%S
  V <- t(eigen(SS)$vectors[,1:(d-noisedim)])
  
  p <- runif(k)^props
  p <- p/sum(p)
  ns <- rowSums(rmultinom(n, 1, p))
  while(sum(ns==0)>0) ns <- rowSums(rmultinom(n, 1, p))
  
  done <- 0
  labels <- c()
  for(i in 1:k){
    X[(done+1):(done+ns[i]),] <- t(matrix(rnorm(ns[i]*(d-noisedim)), nrow = d-noisedim)*runif(d-noisedim)^v/3*(v+1)*scl[i] + MU[i,1:(d-noisedim)])%*%V
    done <- done + ns[i]
    labels <- c(labels, rep(i, ns[i]))
  }
  if(noisedim>0){
    X <- X + matrix(rnorm(nrow(X)*noisedim)/2, nrow = nrow(X))%*%t(eigen(SS)$vectors[,(d-noisedim+1):d])
  } 
  list(x = X, c = labels)
}
