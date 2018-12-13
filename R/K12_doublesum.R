#' Doublesum in K12 estimation
#' The double sum of indicators in K12 estimator. No edge correction or anything fancy.
#'
#' @param xy coordinates
#' @param mark 1-base integer vector giving tpyes
#' @param r ranges
#'
#' @details Use double_integral_ball to compute the global edge correction factor.
#'
#' @export

K12_doublesum <- function(xy, mark, r) {
  doublesum(xy, mark, r)
}

#' Doublesum in g12 estimation
#' The double sum of kernels in g12 estimator. No edge correction or anything fancy.
#' @param xy x,y coordineate matrix
#' @param mark 1/2 mark vector
#' @param r ranges
#' @param bw bandwidth
#' @param kernel kernel, 0:epa 1:box
#' @details
#' unnormalised. Scale with iso_set_cov and pi*r^2 and lambda^2
#' @export

g12_doublesum <- function(xy, mark, r, bw, kernel = 0) {
  g_doublesum(xy, mark, r, bw, kernel)
}
