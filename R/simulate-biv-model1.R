#' Simulate Simple Bivariate Interaction Model
#'
#' Species 1 from Poisson, species 2 from product-shot-noise given species 1.
#'
#' @param n vector of length 2, target point counts
#' @param theta interaction parameter, >= -1
#' @param sigma scale parameter of interaction kernel
#' @param W observation window, spatstat owin-object
#' @param free if TRUE, sample point counts as they were Poisson
#' @param ... passed on to rcox-function
#' @param x1 If a cbind(x,y)-matrix is given, use these as species 1.
#' @param asppp If FALSE (default) return a matrix, else a ppp object
#' @return matrix cbind(x,y,m) with columns for x-coordinate, y-coordinate and 1/2 type mark. if asppp=TRUE, returns a spatstat ppp-object.
#' @import spatstat
#' @export
sim_biv_model <- function(n, theta, sigma, W, free = FALSE, ..., x1 = NULL, asppp=FALSE){
  if(free){
    # redraw n from Poisson
    n <- rpois(2, n)
  }
  # window
  bbox <- with(W, cbind(xrange, yrange))
  # type 1
  x <- if(is.null(x1)) apply(bbox, 2, function(ab) runif(n[1], ab[1], ab[2])) else x1
  # field
  la <- n/prod(apply(bbox,2,diff))
  k0 <- 1/(2*pi*sigma^2)
  a1 <- log(la[2]) - theta * la[1]/k0
  lam <- lambda(x, kernel = "gauss", type = "product", alpha = c(a1, theta), sigma = sigma)
  #e <- coxintensity2matrix(lam)
  #plot(e)
  y <- rcox(lam, bbox = bbox, n = n[2], ...)
  out <- rbind(cbind(x,1), cbind(y$x,2))
  if(asppp) out <- xym2ppp(out, W)
  out
}
