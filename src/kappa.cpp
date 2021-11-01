#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double odds_no_sharing(const double kappa, const NumericVector p, const int ndis) {
  int nsnps=p.size()-1;
  int nd=ndis - 1;
  int maxsnps=50;
  if(maxsnps > nsnps)
    maxsnps = nsnps;
  // Rcout << "nsnps: " << nsnps << "\n";
  double pn=0.0;
  for(int i=0; i<=maxsnps; i++) {
    double tmp = 0.0;
    for(int j=0; j<=maxsnps; j++) {
      double denom = 1 + kappa * ( Rf_choose(nsnps,j) / Rf_choose(nsnps-i,j) - 1);
      tmp = tmp + p(j) / denom;
      // Rcout << i << ' ' << j << ' ' << p(i) * p(j) / denom << "\n";
      // pn = pn + p(i) * p(j) / denom;
    }
    pn = pn + pow(tmp,nd) * p(i);
  }
  return( log(pn) - log(1-pn) );
  // return(pn);
}
