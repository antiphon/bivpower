#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mc_triple_integral_ball(NumericMatrix bbox, NumericVector r, int n) {

  int i,j,k;
  int nr = r.size();
  int dim = bbox.ncol();
  NumericVector out(nr);

  NumericMatrix x(n, 2), y(n,2), z(n, 2);
  double Vol = 1;
  for(i = 0; i < dim; i++) {
    x(_,i) = runif(n, bbox(0,i), bbox(1,i));
    y(_,i) = runif(n, bbox(0,i), bbox(1,i));
    z(_,i) = runif(n, bbox(0,i), bbox(1,i));
    Vol *= bbox(1,i) - bbox(0,i);
  }
  double dy, dz;
  NumericVector sy(nr), sz(nr);
  int oo=0;

  for(i = 0; i < n; i++) {
    checkUserInterrupt();
    for(k =0; k < nr; k++){
      sz(k) = 0;
      sy(k) = 0;
    }
    for(j = 0; j < n; j++) {
      dy = 0;
      dz = 0;
      for(k=0;k<dim;k++){
        dy += pow(x(i,k) - y(j,k),2);
        dz += pow(x(i,k) - z(j,k),2);
      }
      dy = sqrt(dy);
      dz = sqrt(dz);
      oo = 0;
      for(k = nr - 1; k  > -1; k-- ) {
        if(dy < r(k)) sy(k)++;
        else oo = 1;
        if(dz < r(k)) sz(k)++;
        else if(oo) break;
      }
    }
    for(k=0; k < nr; k++) out(k) += sy(k)*sz(k);
  }
  double nn = 1.0/(1.0*n*n*n);
  for(k=0; k < nr; k++) out(k) *= pow(Vol,3) * nn;
  return out;
}

// [[Rcpp::export]]
NumericMatrix mc_triple_integral_ball_matrix(NumericMatrix bbox, NumericVector r, int n) {

  int i,j,k,l;
  int nr = r.size();
  int dim = bbox.ncol();
  NumericMatrix out(nr, nr);

  NumericMatrix x(n, 2), y(n,2), z(n, 2);
  double Vol = 1;
  for(i = 0; i < dim; i++) {
    x(_,i) = runif(n, bbox(0,i), bbox(1,i));
    y(_,i) = runif(n, bbox(0,i), bbox(1,i));
    z(_,i) = runif(n, bbox(0,i), bbox(1,i));
    Vol *= bbox(1,i) - bbox(0,i);
  }
  double dy, dz;
  NumericVector sy(nr), sz(nr);
  int oo=0;

  for(i = 0; i < n; i++) {
    checkUserInterrupt();
    for(k =0; k < nr; k++){
      sz(k) = 0;
      sy(k) = 0;
    }
    for(j = 0; j < n; j++) {
      dy = 0;
      dz = 0;
      for(k=0;k<dim;k++){
        dy += pow(x(i,k) - y(j,k),2);
        dz += pow(x(i,k) - z(j,k),2);
      }
      dy = sqrt(dy);
      dz = sqrt(dz);
      oo = 0;
      for(k = nr - 1; k  > -1; k-- ) {
        if(dy < r(k)) sy(k)++;
        else oo = 1;
        if(dz < r(k)) sz(k)++;
        else if(oo) break;
      }
    }
    for(k=0; k < nr; k++) {
      for(l = k; l < nr; l++) {
        out(k,l) += sy(k)*sz(l);
      }
    }
  }
  double nn = 1.0/(1.0*n*n*n);
  double Vol3 = pow(Vol,3);

  for(k=0; k < nr; k++){
    for(l = k; l < nr; l++) {
      out(k,l) *= Vol3 * nn;
      if(k!=l) out(l,k) = out(k,l);
    }
  }
  return out;
}




// [[Rcpp::export]]
NumericVector mc_triple_integral_annuli(NumericMatrix bbox, NumericVector r, double h, int n) {

  int i,j,k;
  int nr = r.size();
  int dim = bbox.ncol();
  NumericVector out(nr);

  NumericMatrix x(n, 2), y(n,2), z(n, 2);
  double Vol = 1;
  for(i = 0; i < dim; i++) {
    x(_,i) = runif(n, bbox(0,i), bbox(1,i));
    y(_,i) = runif(n, bbox(0,i), bbox(1,i));
    z(_,i) = runif(n, bbox(0,i), bbox(1,i));
    Vol *= bbox(1,i) - bbox(0,i);
  }
  double dy, dz;
  NumericVector sy(nr), sz(nr);

  for(i = 0; i < n; i++) {
    checkUserInterrupt();
    for(k = 0; k < nr; k++){
      sz(k) = 0;
      sy(k) = 0;
    }
    for(j = 0; j < n; j++) {
      dy = 0;
      dz = 0;
      for(k=0;k<dim;k++){
        dy += pow(x(i,k) - y(j,k),2);
        dz += pow(x(i,k) - z(j,k),2);
      }
      dy = sqrt(dy);
      dz = sqrt(dz);
      for(k = 0; k  < nr; k++ ) {
        if(dy < r(k)+h && dy > r(k)-h) sy(k)++;
        if(dz < r(k)+h && dz > r(k)-h) sz(k)++;
      }
    }
    for(k=0; k < nr; k++) out(k) += sy(k)*sz(k);
  }
  double nn = 1.0/(1.0*n*n*n);
  double Vol3 = pow(Vol,3);
  for(k=0; k < nr; k++) out(k) *= Vol3 * nn;
  return out;
}


// [[Rcpp::export]]
NumericMatrix mc_triple_integral_annuli_matrix(NumericMatrix bbox, NumericVector r,
                                               NumericVector h, int n) {

  int i,j,k,l;
  int nr = r.size();
  int dim = bbox.ncol();
  NumericMatrix out(nr, nr);

  NumericMatrix x(n, 2), y(n,2), z(n, 2);
  double Vol = 1;
  for(i = 0; i < dim; i++) {
    x(_,i) = runif(n, bbox(0,i), bbox(1,i));
    y(_,i) = runif(n, bbox(0,i), bbox(1,i));
    z(_,i) = runif(n, bbox(0,i), bbox(1,i));
    Vol *= bbox(1,i) - bbox(0,i);
  }
  double dy, dz;
  NumericVector sy(nr), sz(nr);
  int oo=0;

  for(i = 0; i < n; i++) {
    checkUserInterrupt();
    for(k =0; k < nr; k++){
      sz(k) = 0;
      sy(k) = 0;
    }
    for(j = 0; j < n; j++) {
      dy = 0;
      dz = 0;
      for(k=0;k<dim;k++){
        dy += pow(x(i,k) - y(j,k),2);
        dz += pow(x(i,k) - z(j,k),2);
      }
      dy = sqrt(dy);
      dz = sqrt(dz);
      oo = 0;
      for(k = 0; k  < nr; k++ ) {
        if(dy < r(k)+h(0) && dy > r(k)-h(0)) sy(k)++;
        if(dz < r(k)+h(1) && dz > r(k)-h(1)) sz(k)++;
      }
    }
    for(k=0; k < nr; k++){
      for(l = k; l < nr; l++) {
        out(k,l) += sy(k)*sz(l);
      }
    }
  }

  double nn = 1.0/(1.0*n*n*n);
  double Vol3 = pow(Vol,3);

  for(k=0; k < nr; k++){
    for(l = k; l < nr; l++) {
      out(k,l) *= Vol3 * nn;
      if(l!=k) out(l,k) = out(k,l);
    }
  }
  return out;
}


