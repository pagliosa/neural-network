# ===============================
# About: Empirical Evidence of Phase States Independence
# Dependences: Canon pkg, lasso.r
# Author: Lucas Pagliosa
# Last revision: 02/07/18
# ==============================

require("Canon")

sourceFiles("~/Canon/R/timeSeries.r")
sourceFiles("lasso.r")

testIndependence <- function(ts, m = 2, tau = 1, g = 2e1, l = 1e2)
{
  emb = embedd(ts, m, tau)
  N = nrow(emb)
  half = floor(N / 2)
  
  ret = list()
  ret$sameSet = zeros(m, m)
  ret$interSet = zeros(m, m)
  
  for (i in 1:g)
  {
    set1 = sample(N, half)
    # set1 = 1:half
    set2 = setdiff(1:N, set1)[1:half]
    emb1 = emb[set1,]
    emb2 = emb[set2,]
    data = t(cbind(emb1, emb2))
    mi = makemim(data)
    mi = mi + diag(2 * m)

    ret$sameSet = ret$sameSet + mi[1:m, 1:m]
    ret$interSet = ret$interSet + mi[(m + 1):(2 * m), 1:m]
  }
  ret$sameSet = ret$sameSet / g
  ret$interSet = ret$interSet / g
  
  sameSetT = ret$sameSet[upper.tri(ret$sameSet, FALSE)] 
  
  # ret$sameSetN = mean(sameSetT)
  # ret$interSetN = mean(ret$interSet)
  ret$sameSetN = max(sameSetT)
  ret$interSetN = max(ret$interSet)
  
  ret$dif = ret$sameSetN - ret$interSetN
  
  betas = zeros(m + 1, 2)
  
  for (i in 1:l)
  {
    nod = m + 1
    ps = embedd(ts, nod, tau)
    x = ps[,1:m]
    y = ps[,nod]
    lasso = glmnetLasso(x, y, plot = F)
    betas = betas + cbind(lasso$min$beta, lasso$fes$beta)
  }
  betas = betas / l

  colnames(betas) = c("min", "fes")  
  ret$betas = betas
  
  return(ret)
}

testIndepLoop <- function(ts, maxm = 10, maxd = 10, threshold = 1e-2)
{
  same = NULL
  inter = NULL
  
  for (m in 2:maxm)
  {
    sameVec = c()
    interVec = c()
    for (tau in 1:maxd)
    {
      printf("%d/%d\n", m, tau)
      a = testIndependence(ts, m, tau)
      printMat(a$sameSet)
      printMat(a$interSet)
      sameVec = c(sameVec, a$sameSetN)
      interVec = c(interVec, a$interSetN)
    }
    same = rbind(same, sameVec)
    inter = rbind(inter, interVec)
  }
  ret = list()
  rownames(same) = rownames(inter) = 2:maxm
  colnames(same) = colnames(inter) = 1:maxd
  ret$same = same
  ret$inter = inter
  return(ret)
}
