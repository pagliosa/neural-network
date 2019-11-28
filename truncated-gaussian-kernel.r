# ===============================
# About: Truncated Guassian Kernel
# Dependences: Canon pkg, dwnn.r
# Author: Daniel Cestari, Lucas Pagliosa
# Data of criation: 05/09/2018
# Last revision: 05/09/18
# ==============================

require("Canon")

source("dwnn.r")
loadPackages("kernlab")

inner <- function(a, b)
{
  return(a %*% b)
}

truncatedGaussian <- function(a, b, sigma, order)
{
  gamma = 1/(2 * sigma^2)
  g = 0;
  for (i in 0:order)
  {
    for (k in 0:i)
    {
      for (p in 0:k)
      {
        e1 = i - k
        e2 = k - p
        t1 = ((-gamma)^i * (-2)^p) / (fatorial(e1) * fatorial(e2) * fatorial(p))
        t2 = inner(a, a)^e1 * inner(b, b)^e2 * inner(a, b)^p
        g = g + t1 * t2
      }
    }
  }
  return(g)
}

test <- function(nop = 100, order = 5, C = 1e10)
{
  x0 = runif(1, 0, 1)
  data = createLogistic(nop, onlyts = F, m = 3, tau = 1, x0 = x0)
  # data = createLorenz(nop, onlyts = F, m = 3, tau = 8, start = runif(3, -50, 50))
  sigma = 1 / std(data$ts)
  emb = data$emb
  cc = ncol(emb) #class column
  y = emb[, cc]
  emb = emb[, -cc]
  n = nrow(emb)
  tx = init(n, n, 0)
  d = as.matrix(dist(emb))
  rx = gaussianKernel(d, sigma)
  for (i in 1:n)
  {
    printf("%d\n", i)
    xi = emb[i,]
    for (j in 1:n)
    {
      tx[i, j] = truncatedGaussian(xi, emb[j,], sigma, order)
    }
  }
  rk = ksvm(tx, y, kernel = "matrix", C = C, epsilon = 1e-2, maxiter = -1)
  tk = ksvm(rx, y, kernel = "matrix", C = C, epsilon = 1e-2, maxiter = -1)
  return(list("real" = rk, "truncated" = tk))
}

