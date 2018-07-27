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
  The function get_k(X, maxk, sig, nstart) returns the value of k between 2 and maxk with minimum BIC. Effective degrees of freedom are estimated using sig^2 in place of the residual variance. The functions get_k_elbow(X, maxk, nstart) and get_k_grid(X, maxk, nstart, ngrid) return the BIC minimiser using the two approaches for selecting sig^2 discussed in Hofmeyr, D. "Degrees of freedom and model selection for kmeans clustering", ArXiv preprint [https://arxiv.org/pdf/1806.02034].
}
\examples{
## Generate data with 10 clusters in 10 dimensions
data <- datagen_edfkmeans(1000, 10, 10, var = 3)

## Select the number of clusters to extract using the effective degrees of freedom based BIC estimates
get_k_elbow(data$x, 30)
get_k_grid(data$x, 30)
}
\author{
David P. Hofmeyr
}
