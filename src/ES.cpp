#include <Rcpp.h>
#include <cstdlib>
using namespace Rcpp;

//' Functions that compute Enrichmnet Score
//'
//' @param stat NumericVector
//' @param geneOrd CharacterVector
//' @param Pmiss double
//' @param Nr double
//' @param poz NumericVector
//' @export
//' @rdname ES
// [[Rcpp::export]]
NumericMatrix ES(NumericVector stat, CharacterVector geneOrd, double Pmiss, double Nr, NumericVector poz) {
  int N = geneOrd.size();
  NumericVector P_hit_vec(N);
  NumericVector P_miss_vec(N);
  NumericVector is_in_GS(N);
  NumericVector ES(N);
  double sum_t_stat;
  double Phit;
  double cumSumP = 0;
  double cumSumM = 0;
  NumericVector::iterator itP = &P_hit_vec(0);
  NumericVector::iterator itM = &P_miss_vec(0);
  NumericMatrix ESm(N,3);

  for(int k = 0; k < N; k++){
    int it2 = 0;
    for(NumericVector::iterator it = poz.begin(); it != poz.end(); ++it){
      if(k==*it){
        it2++;
      }
    }
    if(it2==0){
      is_in_GS(k) = 0;
      *itP = cumSumP;
      cumSumM = cumSumM + Pmiss;
      *itM = cumSumM;
    }else{
      is_in_GS(k) = 1;
      sum_t_stat = std::abs(stat(k));
      Phit = sum_t_stat/Nr;
      cumSumP = cumSumP + Phit;
      *itP = cumSumP;
      *itM = cumSumM;
      }
    ES(k) = *itP - *itM;
    ++itP;
    ++itM;
  }

  ESm.column(0) = P_hit_vec ;
  ESm.column(1) = P_miss_vec ;
  ESm.column(2) = ES;

 return ESm;
}
