# Test simulator

devtools::load_all(".")

x <- sim_biv_model(c(50,50), -.2, .1, W=spatstat.geom::square(), asppp=TRUE)


plot(x, cols=1:3)
