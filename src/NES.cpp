#include <Rcpp.h>
using namespace Rcpp;

//' Function to compute NES
//'
//' @param x ES_p vector
//' @param mp mean value of ES_p>=0
//' @param mn mean value of ES_p<0
//' @export
//' @rdname NES
// [[Rcpp::export]]
NumericVector NES(NumericVector x, double mp, double mn) {
  NumericVector  NES_k(x.size());
  int i = 0;
  for(NumericVector::iterator it = x.begin(); it != x.end(); ++it){
    if(*it>0){
      NES_k(i)=*it/mp;
    } else{
      NES_k(i)=*it/mn;
    }
    i++;
  }
  return NES_k;
}
