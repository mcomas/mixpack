#include <Rcpp.h>
#include <limits>

using namespace Rcpp;


double xlog(double x){
  return( x * log(x) );
}

// [[Rcpp::export]]
double lambda_entropy(NumericVector v_tau, int a, int b) {
  return( xlog(v_tau[a] + v_tau[b]) - xlog(v_tau[a]) - xlog(v_tau[b]) );
}

// [[Rcpp::export]]
double lambda_demp(NumericVector v_tau, int a, int b) {
  bool maximum = true;
  for(int i=0; maximum && i<v_tau.size(); maximum = (v_tau[i] <= v_tau[b] ? true : false), i++);
  return( (double) maximum);
}

// [[Rcpp::export]]
double lambda_dempMod(NumericVector v_tau, int a, int b) {
  return( v_tau[b] / (v_tau[a] + v_tau[b]) );
}

// [[Rcpp::export]]
double lambda_coda(NumericVector v_tau, int a, int b) {
  return( log(v_tau[b] / v_tau[a]) );
}

// [[Rcpp::export]]
double lambda_codaNorm(NumericVector v_tau, int a, int b) {
  return( -pow(log(v_tau[b] / v_tau[a]), 2) );
}

// [[Rcpp::export]]
double lambda_prop(NumericVector v_tau, int a, int b) {
  return( v_tau[b] );
}

// [[Rcpp::export]]
double omega_const(NumericVector v_tau, int a, int b) {
  return( 1 );
}

// [[Rcpp::export]]
double omega_prop(NumericVector v_tau, int a, int b) {
  return( v_tau[a] );
}

// [[Rcpp::export]]
double omega_dich(NumericVector v_tau, int a, int b) {
  bool maximum = true;
  for(int i=0; maximum && i<v_tau.size(); maximum = (v_tau[i] <= v_tau[a] ? true : false), i++);
  return( (double) maximum);
}

double confusion(NumericMatrix post, int a, int b, 
                 double (*omega)(NumericVector, int, int), 
                 double (*lambda)(NumericVector, int, int)){
  double numerator = 0, denominator = 0;
  for(int i = 0;i<post.rows();i++){
    NumericVector row = post.row(i);
    double v_omega = omega(row, a, b);
    double v_lambda = lambda(row, a, b);
    numerator +=  v_omega * v_lambda;
    denominator += v_omega;
  }
  return( numerator/denominator );
}

NumericVector optimum(NumericMatrix post, 
                          double (*omega)(NumericVector, int, int), 
                          double (*lambda)(NumericVector, int, int)){
  int m = post.cols();
  double maximum = -std::numeric_limits<double>::max();
  int I = -1, J = -1;
  NumericVector res(2);
  for(int i=0;i<m; i++){
    for(int j=0;j<m;j++){
      double cur = confusion(post, i, j, omega, lambda);
      if(i != j && maximum < cur){
        maximum = cur;
        I = i;
        J = j;
      }
    }
  }
  res(0) = I;
  res(1) = J;
  return(res);
}

// [[Rcpp::export]]
NumericMatrix mergeComponents(NumericMatrix post, int a, int b){
  int m = post.cols();
  int n = post.rows();
  int min_ab = std::min(a,b);
  int max_ab = std::max(a,b);
  NumericMatrix mergeM(post.rows(), post.cols()-1);
  if(a < m && b < m){
    int I = 0;
    for(int i=0; i<min_ab; i++, I++){
      for(int j=0; j<n; j++) mergeM(j,I) = post(j,i);
    }
    for(int j=0; j<n; j++) mergeM(j,I) = post(j,a) + post(j,b);
    I++;
    for(int i=min_ab+1; i<max_ab; i++, I++){
      for(int j=0; j<n; j++) mergeM(j,I) = post(j,i);
    }
    for(int i=max_ab+1; i<m; i++, I++){
      for(int j=0; j<n; j++) mergeM(j,I) = post(j,i);
    }
  }
  return(mergeM);
}

NumericMatrix mergingMatrix(int m, int a, int b){
  int min_ab = std::min(a,b);
  int max_ab = std::max(a,b);
  NumericMatrix mergingM(m, m-1);
  if(a < m && b < m){
    int I = 0;
    for(int i=0; i<min_ab; i++, I++){
      mergingM(i,I) = 1;
    }
    mergingM(a,I) = mergingM(b,I) = 1;
    I++;
    for(int i=min_ab+1; i<max_ab; i++, I++){
      mergingM(i,I) = 1;
    }
    for(int i=max_ab+1; i<m; i++, I++){
      mergingM(i,I) = 1;
    }
  }
  return(mergingM);
}

List mergeStep(NumericMatrix post, 
               double (*omega)(NumericVector, int, int), 
               double (*lambda)(NumericVector, int, int)){
  List out(3);
  NumericVector v = optimum(post, omega_prop, lambda_codaNorm);
  out[0] = mergingMatrix(post.cols(), v(0), v(1));
  out[1] = mergeComponents(post, v(0), v(1));
  return(out);
}


//' Merging components step
//' 
//' @param post Matrix with the posterior probabilities
//' @return partition using prop and codaNorm
//' @export
// [[Rcpp::export]]
List mergeStep_prop_codaNorm(NumericMatrix post){
  return( mergeStep(post, omega_prop, lambda_codaNorm) );
}

//' Merging components step
//' 
//' @param post Matrix with the posterior probabilities
//' @return partition using const and entropy
//' @export
// [[Rcpp::export]]
List mergeStep_const_entropy(NumericMatrix post){
  return( mergeStep(post, omega_const, lambda_entropy) );
}

// [[Rcpp::export]]
double confusion_prop_codaNorm(NumericMatrix post, int a, int b){
  return( confusion(post, a, b, omega_prop, lambda_codaNorm) );
}
  
// [[Rcpp::export]]
double confusion_const_entropy(NumericMatrix post, int a, int b){
  return( confusion(post, a, b, omega_const, lambda_entropy) );
}

// [[Rcpp::export]]
double confusion_prop_demp(NumericMatrix post, int a, int b){
  return( confusion(post, a, b, omega_prop, lambda_demp) );
}

