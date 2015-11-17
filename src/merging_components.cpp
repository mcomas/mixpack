// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <limits>

using namespace Rcpp;


double xlog(double x){
  if(x == 0) return(0);
  else return( x * log(x) );
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

List _mergeStep_(NumericMatrix post, 
               double (*omega)(NumericVector, int, int), 
               double (*lambda)(NumericVector, int, int)){
  List out(3);
  NumericVector v = optimum(post, omega, lambda);
  out[0] = mergingMatrix(post.cols(), v(0), v(1));
  out[1] = mergeComponents(post, v(0), v(1));
  return(out);
}

double (*get_omega(String omega))(NumericVector, int, int) {
  double (*fomega)(NumericVector, int, int);
  
  if( omega == "const" ){
    fomega = omega_const;
  }else if(omega == "prop"){
    fomega = omega_prop;
  }else if(omega == "dich"){
    fomega = omega_dich;
  }
  return(fomega);
}

double (*get_lambda(String lambda))(NumericVector, int, int) {
  double (*flambda)(NumericVector, int, int);
  
  if( lambda == "entropy" ){
    flambda = lambda_entropy;
  }else if(lambda == "demp"){
    flambda = lambda_demp;
  }else if(lambda == "dempMod"){
    flambda = lambda_dempMod;
  }else if(lambda == "coda"){
    flambda = lambda_coda;
  }else if(lambda == "codaNorm"){
    flambda = lambda_codaNorm;
  }else if(lambda == "prop"){
    flambda = lambda_prop;
  }
  return(flambda);
}

//' Merging components step
//' 
//' @param post Matrix with the posterior probabilities
//' @param omega omega function name
//' @param lambda lambda function name
//' @return partition using prop and codaNorm
//' @export
// [[Rcpp::export]]
List mergeStep(NumericMatrix post, String omega = "prop", String lambda = "coda"){
  return( _mergeStep_(post, get_omega(omega), get_lambda(lambda)) );
}


//' Build a hierchical partition from posterior probabilities
//' 
//' This function applies the methodology described in [citar article]
//' to build a hierarchy of classes using the weights or probabilities 
//' that an element belongs to each class
//' @param tau dataframe of probabilities/weights (\code{tau} must be strictly positive)
//' 
//' @param omega function with two parameters (\code{v_tau}, \code{a}). Parameter 
//' \code{v_tau} is a vector of probabilities, parameter \code{a} is the a selected class.
//'\code{omega}(\code{v_tau}, \code{a}) gives the representativeness of element with
//' probabities \code{v_tau} to class \code{a}
//' 
//' @param lambda function with three parameters (\code{v_tau}, \code{a}, \code{b}).
//' Parameter \code{v_tau} is a vector of probabilities, parameters \code{a} and \code{b}
//' are classes to be combined.
//' @export
// [[Rcpp::export]]
List get_hierarchical_partition_fast(NumericMatrix post, String omega = "prop", String lambda = "coda"){
  int LEVEL = post.cols();
  List hp(LEVEL);
  
  arma::mat comb_prev = arma::eye(LEVEL, LEVEL);
  for(int lvl=LEVEL,i=0;0<lvl;lvl--,i++){
    NumericVector v = optimum(post, get_omega(omega), get_lambda(lambda));
    arma::mat comb = Rcpp::as<arma::mat>(mergingMatrix(post.cols(), v(0), v(1)));
    post = mergeComponents(post, v(0), v(1));
    List l_lvl(lvl);
    for(int j=0;j<lvl;j++){
      l_lvl(j) = j;
    }
    hp(i) = l_lvl;
    comb_prev *= comb;
  }
  return(hp);
}

