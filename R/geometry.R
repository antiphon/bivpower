#' Isotropised Set Covariance
#'
#' @param W owin
#' @param R range
#'
#' @details
#'
#' Let gamma_W(x) = |W cap W_x|.
#'
#' Compute E gamma_W(Ru) over u distributed as unif(unit sphere)
#'
#' @export
iso_set_cov <- function(W, R) {
  if(length(R) > 1) return(sapply(R, iso_set_cov, W=W))
  aa <- diff(W$xrange)
  bb <- diff(W$yrange)
  a <- min(aa,bb)
  b <- max(aa,bb)
  A <- a*b # area
  d <- b/a #ratio
  x <- R/sqrt(A/d)
  V <- A/pi
  if(x <= 1) {
    V * ( pi - 2*x - 2*x/d + x^2/d )
  }
  else if(x <= d){
    u <- sqrt(x^2-1)
    V * ( 2*asin(1/x) - 1/d - 2*(x-u) )
  }
  else if( x < sqrt(d^2+1)){
    u <- sqrt(x^2-1)
    v <- sqrt(x^2-d^2)
    V * ( 2*asin((d-u*v)/x^2) + 2*u + 2*v/d - d - (1+x^2)/d )
  }
  else 0
}


#' Double integral over ball inclusion
#'
#' 2D rectangular window.
#'
#' @param W rectangular owin window object
#' @param R radius of the b(o,R)
#' @param r optional lower bound
#'
#' @details
#' Evaluate int_W int_W 1_B(x-y)dxdy where B=b(o,r)
#'
#' @export
double_integral_ball <- function(W, R, r = 0) {
  if(length(R) > 1) return(sapply(R, double_integral_ball, W=W, r = r))
  ff <- function(r) r * iso_set_cov(W, r)
  2 * pi * integrate(ff, r, R)$val
}

#' Triple integral over two ball inclusion of differences
#'
#' Use MC sampling
#'
#' @param W owin
#' @param R radius of ball b(o,R)
#' @param n sample size
#' @param pairs see details
#' @details
#' Let B=b(o,r). Evaluate int_W int_W int_W 1_B(x-y)1_B(x-z)dxdydz
#'
#' If R is a vector and pairs = TRUE, returns a matrix where the integration is with respect to all pairs Bi and Bj of R[i] and R[j].
#' @export

triple_integral_ball_ball_mc <- function(W, R, n = 5000, pairs = FALSE) {
  bbox <- cbind(W$xrange, W$yrange)
  if(!pairs)mc_triple_integral_ball(bbox, R, n)
  else mc_triple_integral_ball_matrix(bbox, R, n)
}



#' Triple integral over two ball inclusion of differences
#'
#' Use MC sampling
#'
#' @param W owin
#' @param R radius of ball b(o,R)
#' @param h annulus halfwidth
#' @param n sample size
#' @param pairs see details
#' @details
#' Let B=b(o,r+h)-b(o,r-h). Evaluate int_W int_W int_W 1_B(x-y)1_B(x-z)dxdydz
#'
#' If R is a vector and pairs = TRUE, returns a matrix where the integration is with respect to all pairs Bi and Bj of R[i] and R[j].
#'
#' @import Rcpp
#' @useDynLib bivpower
#' @export

triple_integral_annulus_annulus_mc <- function(W, R, h, n = 5000, pairs = FALSE) {
  bbox <- cbind(W$xrange, W$yrange)
  if(length(h)<2) h <- c(h,h)
  if(!pairs) mc_triple_integral_annuli(bbox, R, h, n)
  else mc_triple_integral_annuli_matrix(bbox, R, h, n)
}




