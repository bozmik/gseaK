#include <Rcpp.h>
#include <cstdlib>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix ES(NumericVector stat, CharacterVector geneOrd, double Pmiss, double Nr, CharacterVector gs) {
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
    for(CharacterVector::iterator it = gs.begin(); it != gs.end(); ++it){
      if(*it==geneOrd(k)){
        it2++;
      }
    }
    if(it2==0){
     is_in_GS(k) = 0;
    // P_hit_vec(k) = cumSumP;
     *itP = cumSumP;
     cumSumM = cumSumM + Pmiss;
     *itM = cumSumM;
     //P_miss_vec(k) = cumSumM;
    }else{
     is_in_GS(k) = 1;
     sum_t_stat = std::abs(stat(k));
     std::cout <<sum_t_stat <<"\t";
     Phit = sum_t_stat/Nr;
     cumSumP = cumSumP + Phit;
     *itP = cumSumP;
     *itM = cumSumM;
     // P_hit_vec(k) = cumSumP;
     // P_miss_vec(k) = cumSumM;
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
