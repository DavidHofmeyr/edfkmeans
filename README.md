# edfkmeans
Computes effective degrees of freedom for kmeans clustering. This is used within the BIC to perform model selection

To install from within R console:

if(!("devtools"%in%installed.packages())) install.packages("devtools")

library(devtools)

install_github("DavidHofmeyr/edfkmeans")

library(edfkmeans)

help('edfkmeans-package')
