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
  NumericVector pos;
  for(int k = 0; k < gs.size(); k++){
    for(CharacterVector::iterator it = g.begin(); it != g.end(); ++it){
      if(*it==gs(k)){
        pos.push_back(it-g.begin());
        statGS.push_back(stat(it-g.begin()));
        namesGS.push_back(*it);
      }
    }
  }

  return List::create(
    _["names"] = namesGS,
    _["stat"]  = statGS,
    _["position"]  = pos);
}

