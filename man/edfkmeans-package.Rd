\name{edfkmeans-package}
\alias{edfkmeans-package}
\alias{edfkmeans}
\docType{package}
\title{
  Degrees of freedom and model selection for kmeans clustering
}
\description{
  Computes estimates of the effective degrees of freedom in the kmeans model, and uses these within BIC to perform model selection.
}
\details{
  The function kmeans_edf(X, maxk, nstart = 10, ngrid = 30) returns a list containing an estimate of the number of clusters, k, as well as BIC curves for ngrid equally spaced values of sig^2 in place of the within cluster variance. Effective degrees of freedom are estimated, and subsequent model selection is performed as described in Hofmeyr, D. "Degrees of freedom and model selection for kmeans clustering", ArXiv preprint [https://arxiv.org/pdf/1806.02034]. In particular, the most frequent first minimiser of the multiple smoothed BIC curves is given as the estimated number of clusters.
}
\examples{
## Generate data with 10 clusters in 10 dimensions
data <- datagen_edfkmeans(1000, 10, 10, var = 4)
  
## The first three arguments are the number of data, dimensions and clusters respectively. The cluster overlap can be increased
## by increasing the var argument. The mixing proportions can be varied by setting argument props to a positive value.
## The individual cluster variances can be varied by setting the scale argument to a positive value. 
## The cluster shapes can be varied by setting the argument v to a positive value.
## Noise dimensions can be included using the noisedim argument.
  
## Compute BIC estimates using effective degrees of freedom
sol <- kmeans_edf(data$x, 30)

## The estimated number of clusters is
sol$k

## Plots of the BIC curves can be inspected for alternative appropriate values of k
plot(sol$bic[,1], ylim = c(min(sol$bic), max(sol$bic)), col = rgb(.5, .5, .5, .5), type = 'l')
for(i in 2:30) lines(sol$bic[,i], col = rgb(.5, .5, .5, .5))


}
\references{
Hofmeyr, D. (2018) Degrees of freedom and model selection for kmeans clustering, \emph{ArXiv preprint}, ArXiv 1806.02034.
}
\author{
David P. Hofmeyr
}
