source("code/pBIL.R")
## get Noam Ross's profiling code
source("https://raw.githubusercontent.com/noamross/noamtools/refs/heads/master/R/proftable.R")

kds2df <- function(obs, pred, h = 50000) {
  names(obs) <- tolower(names(obs))
  names(pred) <- tolower(names(pred))
  nobs <- length(obs$x)
  npred <- length(pred$x)
  
  #h <- with(obs, c(bandwidth.nrd(x), bandwidth.nrd(y))/4)
  
  ax <- outer(obs$x, pred$x, "-")/h
  ay <- outer(obs$y, pred$y, "-")/h

  n <- length(ax)
  k <- exp(-(ax^2+ay^2)/2)/(2*pi)
  az <- colSums(k)/(nobs*h^2)
  
  result <- pred
  result$z <- (az/sum(az))*npred
  return(result)
}

nobs <- 10000
npred <- 5000
## want to profile ksd2d()
mkdat <- function(n) {
  runif(n) |> matrix(ncol = 2, dimnames = list(NULL, c("x", "y"))) |>
    as.data.frame()
}
set.seed(101)
obs <- mkdat(nobs)
pred <- mkdat(npred)
kk0 <- kds2d(obs, pred)

Rprof("kds2d_prof1.Rout", line.profiling = TRUE)
invisible(kds2df(obs, pred))
Rprof(NULL)
proftable("kds2d_prof1.Rout", lines = 20)

system.time(
  kk1 <- kds2df(obs, pred)
)
stopifnot(all.equal(kk0, kk1))
