#' deviation test
#'
#' @export
test_all <- function(e, rok = NULL) {
  f <- e$est_data
  f0 <- e$est_M0
  # use f0 to estimate mean
  s <- apply(f0, 1, sd)
  if(is.null(rok)) rok <- rep(TRUE, length(s))
  # skip s=0
  ok <- s>0 & is.finite(s)
  rok <- rok[ok]
  s <- s[ok]
  f <- f[ok]
  f0 <- f0[ok,]
  fm <- rowMeans(f0)
  m <- ncol(f0)
  # quantiles
  dq <- apply(f0, 1, quantile, prob = c(0.025, 0.975))

  # diffs
  absdelta <- abs(f - fm)
  absdelta0 <- abs(f0 - fm)
  l2delta <- (f - fm)^2
  l2delta0 <- (f0 - fm)^2
  #
  # pointwise
  m1 <- rowSums(l2delta0 >= l2delta)
  pv1 <- (1+m1)/(1+m)
  # add NA to s=0
  pv <- 0*s+NA
  pv[ok] <- pv1

  # Deviances
  l2delta <- l2delta[rok]
  l2delta0 <- l2delta0[rok,]
  # L2 deviance test
  deltaS <- sum(l2delta)
  delta0S <- colSums(l2delta0)
  mS <-  (1+sum(delta0S >= deltaS))/(1+m)
  #
  # studentised
  sdeltaS <- sum(l2delta/s^2)
  sdelta0S <- colSums(l2delta0/s^2)
  smS <-  (1+sum(sdelta0S >= sdeltaS))/(1+m)
  #
  # # directional quantiles
  # lo <- abs( dq[1,]-m)
  # hi <- abs(dq[2,]-m)
  # dql2delta <- l2delta * ( (l2delta < 0)/lo + (l2delta>=0)/hi)
  # dql2delta0 <- l2delta0 * ( (l2delta0 < 0)/lo + (l2delta0>=0)/hi)
  #
  # dqdeltaS <- sum(dql2delta)
  # dqdeltaS0 <- colSums(dql2delta0)
  # dqmS <-  (1+sum(dqdeltaS0 >= dqdeltaS))/(1+m)
  #
  # #
  # # mad
  # do <- max(absdelta)
  # do0 <- apply(absdelta0, 2, max)
  # mad <-  (1+sum(do0 >= do))/(1+m)
  # #
  # # mad std
  # do <- max(absdelta/s)
  # do0 <- apply(absdelta0/s, 2, max)
  # smad <-  (1+sum(do0 >= do))/(1+m)
  # #

  # thats enough.
  list(pointwise = pv, dev2 = mS, dev2_st = smS) #, mad = mad, mad_st = smad, dev2_dq = dqmS)
}

#' Power per list of p-values from test-all
#' @export
testlist_power <- function(pl, alpha = 0.05) {
  nam <- names(pl[[1]])
  pol <- list()
  for(n in nam) {
    pe <- rbind( sapply(pl, getElement, n) )
    po <- rowMeans(pe <= alpha, na.rm = TRUE)
    pol[[n]] <- po
  }
  pol
}

