#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace Rcpp;

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

