#include <algorithm>
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
NumericVector mv_mult(NumericMatrix X, NumericVector beta){
  int n = X.nrow();
  NumericVector Xbeta (n);
  for(int i = 0; i<n ; i++){
    Xbeta[i] = sum(X(i,_) * beta);
  }
  return Xbeta;
}


// [[Rcpp::export]]
NumericVector expect_cum_weibull_cpp(Nullable<NumericVector> t_,
                                     double k = 1,
                                     double b = 1){
  if(t_.isNull()){
    return  NumericVector::create();
  }else{
    NumericVector t;
    t = t_;
    return b * pow(t,k);
  }
}


// [[Rcpp::export]]
NumericVector expect_cum_weibull_inverse_cpp(Nullable<NumericVector> m_,
                                             double k = 1,
                                             double b = 1){
  if(m_.isNull()){
    return  NumericVector::create();
  }else{
    NumericVector m;
    m = m_;
    return pow(m/b,1/k);
  }
}



double expect_cum_weibull_tvc_single_cpp(double t,
             NumericVector lin_pred,
             NumericVector t_breaks,
             double k,
             double b){
  NumericVector t_breaks_use = t_breaks[Range(0, sum(t_breaks < t))];
  int t_breaks_use_length  = t_breaks_use.size();
  NumericVector lin_pred_use = lin_pred[Range(0, t_breaks_use_length - 1)];
  t_breaks_use[t_breaks_use_length - 1] = t;
  NumericVector t_breaks_base = expect_cum_weibull_cpp(t_breaks_use, k, b);
  t_breaks_base.push_front(0);
  return sum(diff(t_breaks_base) * exp(lin_pred_use)) ;

}

// // [[Rcpp::export]]
// double mysum(LogicalVector x) { return(sum(x)); }

// [[Rcpp::export]]
NumericVector expect_cum_weibull_tvc_cpp(Nullable<NumericVector> t_,
                                         NumericVector lin_pred = 0,
                                         NumericVector t_breaks = NumericVector::create(),
                                         double k = 1,
                                         double b = 1){
  if(t_.isNull()){
    return  NumericVector::create();
  }

  NumericVector t;
  t = t_;
  if(t.size() == 0){
    return t;
  }

  if(is_true(any(t < 0))){
    stop("some t is less than 0");
  }

  if (t_breaks.size() == 0)
  {
    t_breaks = R_PosInf;
  }

  if(is_true(any(t > tail(t_breaks, 1)[0]))){
    stop("some t is larger than the max covariate node");
  }

  if(lin_pred.size() != t_breaks.size()){
    stop("length of linear effects values is not equal to that of nodes");
  }

  NumericVector m = sapply(t,
                           [&](double tt){return expect_cum_weibull_tvc_single_cpp(tt,lin_pred, t_breaks, k, b);});

  return m;

}



double expect_cum_weibull_tvc_inverse_single_cpp(double m,
                                          NumericVector lin_pred,
                                          NumericVector t_breaks,
                                          NumericVector m_breaks,
                                          double k,
                                          double b){
  int left_index = sum(m > m_breaks) - 1;
  // double t_left = t_breaks[left_index];
  // double m_left = m_breaks[left_index];
  return pow( ( (m - m_breaks[left_index]) / (b * exp(lin_pred[left_index])) +
              pow(t_breaks[left_index], k) ), 1 / k) ;
}



// [[Rcpp::export]]
NumericVector expect_cum_weibull_tvc_inverse_cpp(Nullable<NumericVector> m_,
                                                 NumericVector lin_pred = 0,
                                                 NumericVector t_breaks = NumericVector::create(),
                                                 double k = 1,
                                                 double b = 1){
  if(m_.isNull()){
    return  NumericVector::create();
  }

  NumericVector m;
  m = m_;
  if(m.size() == 0){
    return m;
  }

  if(is_true(any(m < 0))){
    stop("some m is less than 0");
  }

  if (t_breaks.size() == 0)
  {
    t_breaks = R_PosInf;
  }

  if(lin_pred.size() != t_breaks.size()){
    stop("length of linear effects values is not equal to that of nodes");
  }

  NumericVector m_mid_vals = expect_cum_weibull_cpp(t_breaks, k, b);
  m_mid_vals.push_front(0);
  NumericVector m_breaks = cumsum(diff(m_mid_vals) * exp(lin_pred));

  m_breaks.push_front(0);
  t_breaks.push_front(0);

  if(is_true(any(m > tail(m_breaks, 1)[0]))){
    stop("some m is larger than the max covariate node corresponding value");
  }




  return sapply(m,
                [&](double mm){return expect_cum_weibull_tvc_inverse_single_cpp(mm, lin_pred, t_breaks, m_breaks, k,b);});

}
// [[Rcpp::export]]
std::vector<double> generate_events_poisson1(double t_start, double t_end){

  double t;
  std::vector<double> z;
  z.reserve(100);
  t = t_start;
  // z.clear();

  while(t <= t_end){
    t = t - log(R::runif(0, 1));
    z.push_back(t);
  }

  z.pop_back();
  return z;
}

// [[Rcpp::export]]
std::vector<double> generate_events_no0_poisson1(double t_start, double t_end){

  double t;
  std::vector<double> z;
  z.reserve(100);
  t = t_start;
  // z.clear();
  t = t - log(1 - (R::runif(0, 1)) * (1 - exp(t_start - t_end)));
  z.push_back(t);

  while(t <= t_end){
    t = t - log(R::runif(0, 1));
    z.push_back(t);
  }

  z.pop_back();
  return z;
}


// [[Rcpp::export]]
std::vector<double> ltrgamma(int n, NumericVector shapes, NumericVector scales, double truncate){
  std::vector<double> values(n);
  double value;
  // values.reserve(2*n);
  int tries = 0;

  for(int i = 0; i < n;i++){
    tries = 0;
    do{
      tries++;
      value = R::rgamma(shapes[i],scales[i]);
      if(tries > 2E7){
        stop("Maximum number of attempts reached for ltrgamma");
      }

      // Rcout << value  << truncate <<"\n";

    } while (value < truncate);

    values[i] = value;

  }

  return values;
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# expect_cum_weibull_cpp(c(),1,1)
# expect_cum_weibull_inverse_cpp(c(1,2),2,4)
# # expect_cum_weibull_tvc_cpp(0.5, t_breaks = c(1,2))
# # mysum(rep(T,3))
# # expect_cum_weibull_tvc_single_cpp(2, c(0.5,1.5,2), c(1, 3,5), k = 2, b = 3)
# expect_cum_weibull_tvc(c(2,3), c(0.5,1.5,2), c(1, 3,5), k = 2, b = 3)
# expect_cum_weibull_tvc_cpp(c(2,3), c(0.5,1.5,2), c(1, 3,5), k = 2, b = 3)
# expect_cum_weibull_tvc_cpp(c(), c(0.5,1.5,2), c(1, 3,5), k = 2, b = 3)
# expect_cum_weibull_tvc(c(), c(0.5,1.5,2), c(1, 3,5), k = 2, b = 3)
# # expect_cum_weibull_tvc_inverse_single_cpp(25, c(0), Inf, Inf ,k = 2, b = 3)
# # expect_cum_weibull_tvc_inverse(25, c(0), Inf, Inf, k = 2, b = 3)
# expect_cum_weibull_tvc_inverse_cpp(c(25,36), c(0.5,1.5,2), c(1, 3,5), k = 2, b = 3)
# expect_cum_weibull_tvc_inverse(c(25,36), c(0.5,1.5,2), c(1, 3,5), k = 2, b = 3)
# expect_cum_weibull_tvc_inverse_cpp(c(), c(0.5,1.5,2), c(1, 3,5), k = 2, b = 3)
# generate_events_no0_poisson1(0,0.0001)
# generate_events_poisson1(0, 0.0000001)
*/
