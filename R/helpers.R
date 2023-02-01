
#' Id tag generator, pad 000
#'
#' @export
makeid <- function(k, pref = "")  paste0("id_", pref, "_", substr(10000000+k, 2, 9))

#' xym matrix to ppp
#'
#' @importFrom spatstat.geom ppp
#' @export
xym2ppp <- function(xym, W) {
  ppp(xym[,1], xym[,2], window = W, marks = factor(xym[,3]))
}

#' Simplify a ppp object to x,y,m matrix
#'
#' @export
ppp2xym <- function(pp) {
  with(pp, cbind(x, y, as.integer(marks)) )
}

#' List of pps to just x,y,m matrix, for storage
#'
#' @export
ppl2xy <- function(ppl){
  lapply(ppl, function(x) ppp2xym(x))
}
