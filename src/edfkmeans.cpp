#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]

arma::mat Pdist(arma::mat A, arma::mat B) {
  
  arma::colvec An =  sum(square(A),1);
  arma::colvec Bn =  sum(square(B),1);
  
  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  
  return(C); 
}

// [[Rcpp::export]]

NumericVector edfkmeans(arma::mat X, arma::mat M, IntegerVector asg, IntegerVector ns, double SS, int nclust, int n, int d, NumericVector sig, int nsig){
  arma::mat ds = Pdist(X, M);
  NumericVector jumps(n*d*(nclust-1));
  NumericVector offsets(n*d*(nclust-1));
  int count = 0;
  double a, b, c;
  double del1, del2, del;
  for(int i=0; i<n; i++){
    for(int j=0; j<d; j++){
      for(int cl=0; cl<nclust; cl++){
        if(cl!=asg[i]){
          a = 1.0-pow((ns[asg[i]]-1.0)/ns[asg[i]],2);
          b = 2.0*(X(i,j)-M(cl,j)-(X(i,j)-M(asg[i],j))*(ns[asg[i]]-1.0)/ns[asg[i]]);
          c = ds(i,cl)-ds(i,asg[i]);
          if((b*b-4*a*c)>=0){
            del1 = (-b+pow(b*b-4*a*c,0.5))/2.0/a;
            del2 = (-b-pow(b*b-4*a*c,0.5))/2.0/a;
            if(fabs(del1)<fabs(del2)) del = del1;
            else del = del2;
            if(del<0){
              jumps[count] = M(asg[i],j)-M(cl,j)*ns[cl]/(ns[cl]+1.0)-X(i,j)/(ns[cl]+1.0)+del*(ns[cl]-ns[asg[i]]+1.0)/(ns[asg[i]]*(ns[cl]+1.0));
              offsets[count] = pow(X(i,j)+del-M(asg[i],j),2);
              count+=1;
            }
            else{
              jumps[count] = -(M(asg[i],j)-M(cl,j)*ns[cl]/(ns[cl]+1.0)-X(i,j)/(ns[cl]+1.0)+del*(ns[cl]-ns[asg[i]]+1.0)/(ns[asg[i]]*(ns[cl]+1.0)));
              offsets[count] = pow(X(i,j)+del-M(asg[i],j),2);
              count+=1;
            }
          }
        }
      }
    }
  }
  NumericVector edfs(nsig);
  double sum;
  for(int i=0; i<nsig; i++){
    sum = 0;
    for(int j=0; j<count; j++) sum += jumps[j]*exp(-offsets[j]/sig[i]/sig[i]/2.0);
    sum *= 1.0/sig[i]/pow(2*3.1416, 0.5);
    edfs[i] = nclust*d+sum;
  }
  return(edfs);
}

