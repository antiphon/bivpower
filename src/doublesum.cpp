#include <Rcpp.h>
using namespace Rcpp;


// For K12 function

// [[Rcpp::export]]
NumericVector doublesum(NumericMatrix xy, NumericVector mark, NumericVector r) {
  int n = xy.nrow();
  int dim = xy.ncol();
  int nr = r.size();
  NumericVector out(nr);
  int i,j,k;
  double d,d2;
  double rmax = max(r);
  double rmax2 = rmax*rmax;
  
  for(i = 0; i < n-1; i++)
    for(j = i + 1 ; j < n; j++) {
      if(mark(i)!=mark(j)){
        d2 = 0;
        for(k = 0; k < dim; k++) d2 += pow(xy(i,k)-xy(j,k),2);
        if(d2 < rmax2){
          d = sqrt(d2);
          for(k = nr-1; k > -1; k--) {
            if(d < r(k)) out(k) += 1.0;
            else break;
          }
        }
      }
    }
  return out;
}



// For K11 function
// [[Rcpp::export]]
NumericVector doublesum22(NumericMatrix xy, NumericVector r) {
  int n = xy.nrow();
  int dim = xy.ncol();
  int nr = r.size();
  NumericVector out(nr);
  int i,j,k;
  double d,d2;
  double rmax = max(r);
  double rmax2 = rmax*rmax;
  
  for(i = 0; i < n-1; i++)
    for(j = i + 1 ; j < n; j++) {
      d2 = 0;
      for(k = 0; k < dim; k++) d2 += pow(xy(i,k)-xy(j,k),2);
      if(d2 < rmax2){
        d = sqrt(d2);
        for(k = nr-1; k > -1; k--) {
          if(d < r(k)) out(k) += 1.0;
          else break;
        }
      }
    }
  return out;
}


/*  For g12 function
 *  Epanechnikov kernel, one bandwidth
 * 
 * 
 */ 

double epa(double d) {
  if(d < 0) d *= -1.0;
  if(d > 1) return 0;
  return 0.75 * (1 - d*d);
}

double box(double d) {
  if(d < 0) d *= -1.0;
  if(d > 1) return 0;
  return 0.5;
}


// [[Rcpp::export]]
NumericVector g_doublesum(NumericMatrix xy, NumericVector mark, NumericVector r, double bw, int kern = 1) {
  int n = xy.nrow();
  int dim = xy.ncol();
  int nr = r.size();
  NumericVector out(nr);
  int i,j,k;
  double d,d2;
  double rmax = max(r) + bw;
  double rmax2 = rmax*rmax;
  
  double (*kfun)(double) = &epa;
  if(kern == 1) kfun = &box;
  
  for(i = 0; i < n-1; i++)
    for(j = i + 1 ; j < n; j++) {
      if(mark(i)!=mark(j)){
        d2 = 0;
        for(k = 0; k < dim; k++) d2 += pow(xy(i,k)-xy(j,k),2);
        if(d2 < rmax2){
          d = sqrt(d2);
          for(k = 0; k < nr; k++) {
            out(k) += kfun((d-r(k))/bw)/bw;
          }
        }
      }
    }
    return out;
}


// [[Rcpp::export]]
NumericMatrix g_doublesum_asym(NumericMatrix xy, NumericVector mark, NumericVector r, double bw12, double bw21, int kern = 1) {
  int n = xy.nrow();
  int dim = xy.ncol();
  int nr = r.size();
  NumericMatrix out(nr,2);
  int i,j,k;
  double d,d2;
  double rmax = max(r) + fmax(bw12,bw21);
  double rmax2 = rmax*rmax;
  
  double (*kfun)(double) = &epa;
  if(kern == 1) kfun = &box;

  for(i = 0; i < n-1; i++)
    for(j = i + 1 ; j < n; j++) {
      if(mark(i) != mark(j)){
        d2 = 0;
        for(k = 0; k < dim; k++) d2 += pow(xy(i,k)-xy(j,k),2);
        if(d2 < rmax2){
          d = sqrt(d2);
          for(k = 0; k < nr; k++) {
            out(k, 0) += kfun((d-r(k))/bw12)/bw12;
            out(k, 1) += kfun((d-r(k))/bw21)/bw21;
          }
        }
      }
    }
    return out;
}


