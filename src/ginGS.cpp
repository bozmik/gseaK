#include <Rcpp.h>
using namespace Rcpp;

//' ginGS
//'
//' @param g CharacterVector
//' @param gs CharacterVector
//' @param stat NumericVector
//' @export
//' @rdname ginGS
// [[Rcpp::export]]
List ginGS(CharacterVector g, CharacterVector gs, NumericVector stat) {
  NumericVector statGS;
  CharacterVector namesGS;
  List res;
  for(int k = 0; k < gs.size(); k++){
    for(CharacterVector::iterator it = g.begin(); it != g.end(); ++it){
      if(*it==gs(k)){
      statGS.push_back(stat(it-g.begin()));
      namesGS.push_back(*it);
    }
  }
  }
  res["names"]=namesGS;
  res["stat"] = statGS;
 return res;
}

