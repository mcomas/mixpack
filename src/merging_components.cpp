// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <limits>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

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
  double logA = log(v_tau[a]);
  double logB = log(v_tau[b]);
  return(2 * logA * logB -logB*logB - logA * logA);
  //return( -pow(log(v_tau[b] / v_tau[a]), 2) );
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

List optimum(NumericMatrix post, 
                          double (*omega)(NumericVector, int, int), 
                          double (*lambda)(NumericVector, int, int)){
  int m = post.cols();
  double maximum = -std::numeric_limits<double>::max();
  int I = -1, J = -1;
  List out(2);
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
  out(0) = res;
  out(1) = maximum;
  return(out);
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

// [[Rcpp::export]]
NumericMatrix mergeComponents_mult(NumericMatrix post, int a, int b){
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
    for(int j=0; j<n; j++) mergeM(j,I) = post(j,a) * post(j,b);
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

arma::mat _merge_step_(NumericMatrix post, 
               double (*omega)(NumericVector, int, int), 
               double (*lambda)(NumericVector, int, int)){
  int k = post.cols();
  arma::mat conf_matrix = arma::eye(k, k);
  for(int i=0;i<k;i++){
    for(int j=0;j<k;j++){
      conf_matrix(i,j) = confusion(post, i, j, omega, lambda);
    }
  }
  return(conf_matrix);
}

double (*get_omega(String omega))(NumericVector, int, int) {
  double (*fomega)(NumericVector, int, int) = NULL;
  if( omega == "cnst" ){
    fomega = omega_const;
  }else if(omega == "prop"){
    fomega = omega_prop;
  }else if(omega == "dich"){
    fomega = omega_dich;
  }
  return(fomega);
}

double (*get_lambda(String lambda))(NumericVector, int, int) {
  double (*flambda)(NumericVector, int, int) = NULL;
  
  if( lambda == "entr" ){
    flambda = lambda_entropy;
  }else if(lambda == "demp"){
    flambda = lambda_demp;
  }else if(lambda == "demp.mod"){
    flambda = lambda_dempMod;
  }else if(lambda == "coda"){
    flambda = lambda_coda;
  }else if(lambda == "coda.norm"){
    flambda = lambda_codaNorm;
  }else if(lambda == "prop"){
    flambda = lambda_prop;
  }
  return(flambda);
}

// [[Rcpp::export]]
arma::mat merge_step_cpp(NumericMatrix post, String omega, String lambda){
  return( _merge_step_(post, get_omega(omega), get_lambda(lambda)) );
}


// [[Rcpp::export]]
List get_hierarchical_partition_cpp(NumericMatrix post, String omega = "prop", String lambda = "coda"){
  int LEVEL = post.cols();
  List hp(LEVEL);
  
  double (*flambda)(NumericVector, int, int) = get_lambda(lambda);
  double (*fomega)(NumericVector, int, int) = get_omega(omega);
  
  arma::mat comb_level = arma::eye(LEVEL, LEVEL);
  NumericMatrix post_level = NumericMatrix(post);
  NumericVector values = NumericVector(LEVEL);
  NumericVector v;
  List out;
  for(int lvl=LEVEL,l=0;1<lvl;lvl--,l++){
    List l_lvl(lvl);
    for(int j=0;j<lvl;j++){
      std::vector<double> vec;
      for(int i=0;i<LEVEL;i++) if( comb_level(i,j) == 1 ) vec.push_back(i+1);
      l_lvl(j) = wrap(vec);
    }
    int cur = LEVEL-l-1;
    hp(cur) = l_lvl;
    if(cur == 1) break;
    out = optimum( post_level, fomega, flambda );
    v =  out(0);
    values[cur-1] = out(1);
    post_level = mergeComponents(post_level, v[0], v[1]);
    comb_level = comb_level * Rcpp::as<arma::mat>( mergingMatrix(lvl, v[0], v[1]) );
    Rcpp::checkUserInterrupt();
  }
  out = optimum( post_level, fomega, flambda );
  values[0] = out(1);
  std::vector<double> vec(0);
  for(int i=0;i<LEVEL;i++) vec.push_back(i+1);
  List l_lvl(1);
  l_lvl(0) = wrap(vec);
  hp(0) = l_lvl;
  hp(1) = values;
  return(hp);
}


// [[Rcpp::export]]
List get_hierarchical_partition_mult_fast(NumericMatrix post, String omega = "prop", String lambda = "coda"){
  int LEVEL = post.cols();
  List hp(LEVEL);
  
  double (*flambda)(NumericVector, int, int) = get_lambda(lambda);
  double (*fomega)(NumericVector, int, int) = get_omega(omega);
  
  arma::mat comb_level = arma::eye(LEVEL, LEVEL);
  NumericMatrix post_level = NumericMatrix(post);
  
  for(int lvl=LEVEL,l=0;1<lvl;lvl--,l++){
    List l_lvl(lvl);
    for(int j=0;j<lvl;j++){
      std::vector<int> vec;
      for(int i=0;i<LEVEL;i++) if( comb_level(i,j) == 1 ) vec.push_back(i+1);
      l_lvl(j) = wrap(vec);
    }
    hp(l) = l_lvl;
    
    List opt = optimum( post_level, fomega, flambda );
    NumericVector v = opt[0];
    
    post_level = mergeComponents_mult(post_level, v(0), v(1));
    comb_level *= Rcpp::as<arma::mat>( mergingMatrix(lvl, v(0), v(1)) );
  }
  
  return(hp);
}

