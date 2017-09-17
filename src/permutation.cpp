#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace Rcpp;

//' permutation
//'
//' @param x NumericVector
//' @param n int
//' @export
//' @rdname permutation
// [[Rcpp::export]]
NumericMatrix permutation(NumericVector x, int n) {
  int s = x.size();
  NumericMatrix y(n,s);
  for(int i = 0; i < n; i++){
    std::random_shuffle( x.begin(), x.end() );
    y(i,_) = x;
  }
 return y;

}

