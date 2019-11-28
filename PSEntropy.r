# ===============================
# About: Phase space entropy related functions
# Dependences: none
# Author: Lucas Pagliosa
# Last revision 19/10/15
# ==============================

require("Canon")

cofMatrix <- function(mat)
{
  cofMat = NULL
  for (i in 1:nrow(mat))
  {
    cofVec = c()
    for (j in 1:ncol(mat))
      cofVec = c(cofVec, cofactor(mat, i, j))
    cofMat = rbind(cofMat, cofVec)
  }
  return(cofMat)
}

# Maximize this one
openBallsEntropy <- function(phaseSpace, epsilon = 1e-3)
{
  if (missing(phaseSpace))
    stop("Phase space required")
  
  nattr = ncol(phaseSpace) - 1
  nstates = nrow(phaseSpace)
  knownDimensions = phaseSpace[, 1:nattr]
  euclidean = as.matrix(dist(knownDimensions))
  adjacency = c()
  
  for (i in 1:nstates)
  {
    # Probability
    ids = which(euclidean[i,] < epsilon)
    
    # Updating the probabilities
    adjacency = c(adjacency, length(ids))
  }
  return(shannonEntropy(adjacency))
}

# Minimize this one -> not a good one, though
futureStateEntropy <- function(phaseSpace, epsilon = 1e-3)
{
  if (missing(phaseSpace))
    stop("Phase space required")
  
  nattr = ncol(phaseSpace)
  nstates = nrow(phaseSpace)
  futureState = phaseSpace[, 1]
  adjacency = c()
  
  for (i in 1:nstates)
  {
    # Probability
    dist = abs(futureState - phaseSpace[i, 2])
    ids = which(dist <= epsilon)
    
    # Updating the probabilities
    adjacency = c(adjacency, length(ids))
  }
  return(shannonEntropy(adjacency))
}

vonNeumannEntropy <- function(phaseSpace, k = 1)
{
  values = getEingValues(phaseSpace, normalize = F)$values
  # cat("eig:", values, "\n")
  sh = shannonEntropy(values[1:k])
  return(sh) # / length(values))
}

# Nao tem sentido, mas testando igual
covarianceEntropy <- function(phaseSpace)
{
  if (missing(phaseSpace))
    stop("Phase space required")
  
  nattr = ncol(phaseSpace)
  nstates = nrow(phaseSpace)
  
  covMatrix = cov(phaseSpace)
  return(vonNeumannEntropy(covMatrix))
}

numberOfIterations <- function(maxIter = 2, iter)
{
  return(ifelse(iter <= maxIter, TRUE, FALSE))
}

numberOfInstances <- function(vec, n = 100, test = FALSE)
{
  noe = length(vec)
  
  if (test)
    printf("noe: %d\n", noe)
  if (noe <= n || n == 0)
    return(FALSE)
  return(TRUE)
}

subdivedePhaseSpaceAxis <- function(vec, values, eta, criterionFunc, iter, range)
{
  if (!criterionFunc(values, iter, range))
    return(vec)

  half = sum(range) / 2
  leftRange = c(range[1], half)
  rightRange = c(half, range[2])
  left = values[values < half]
  right = values[values >= half]

  if (FALSE)
  {
    cat(values, "\n")
    cat(vec, "\n")
    cat(half, "\n")
    cat(left, "\n")
    cat(right, "\n")
  }
  
  vec = c(vec, half)
  vec = subdivedePhaseSpaceAxis(vec, left, eta, criterionFunc, 
    iter + 1, leftRange)
  vec = subdivedePhaseSpaceAxis(vec, right, eta, criterionFunc, 
    iter + 1, rightRange)
  return(vec)
}

subdvidePhaseSpace <- function(phaseSpace, criterionFunc = numberOfInstances)
{
  numberOfDimensions = ncol(phaseSpace)
  numberOfElements = nrow(phaseSpace)
  eta = 1e-6
  list = NULL

  for (i in 1:numberOfDimensions)
  {
    values = phaseSpace[,i]
    range = range(values)
    vec = c(range[1] - eta, range[2] + eta)
    if (range[1] == range[2])
      ret = min(values)
    else
      ret = subdivedePhaseSpaceAxis(vec, values, eta, criterionFunc, 1, range)
    list = addList(list, sort(ret))
  }
  return(list)
}

validateCriterion <- function(dimVec, values, criterionFunc = numberOfInstances)
{
  printf("Testing criterion:\n")
  for (j in 2:length(dimVec))
  {
    b = dimVec[j - 1]
    e = dimVec[j]
    criterionFunc(values[values > b & values <= e], test = TRUE)
  }
}

getCritFunc <- function(phaseSpace, p, maxIter)
{
  critFunc = NULL
  
  if (missing(p) || is.null(p))
  {
    if (maxIter < 1)
      maxIter = 1
    critFunc <- function(values, iter, range)
    {
      numberOfIterations(maxIter, iter)
    }
    printf("Max number of iteration: %d\n",  maxIter)
  }
  else
  {
    noe = nrow(phaseSpace)
    noi = round(noe * p)
    critFunc <- function(values, iter, range)
    {
      (numberOfInstances(values, n = noi, F) && range[1] != range[2]
        && min(values) != max(values))
    }
    printf("Number of instances: %d (%g%% of %d)\n", noi, p, noe)
  }
  return(critFunc)
}

testSubdivision <- function(phaseSpace = rand(10, 2), p, maxIter = 10)
{
  phaseSpace = as.matrix(phaseSpace)
  sortPhaseSpace = apply(phaseSpace, 2, sort)
  
  par(pty = "s")
  if(ncol(phaseSpace) == 1)
    plot(cbind(phaseSpace, phaseSpace))
  else
    plot3d(phaseSpace)
  
  criterionFunc = getCritFunc(phaseSpace, p, maxIter)
  list = subdvidePhaseSpace(sortPhaseSpace, criterionFunc)
  numberOfDims = length(list)
  numberOfDraws = ifelse(numberOfDims > 3, 3, numberOfDims)
  for (i in 1:numberOfDraws)
  {
    eval(parse(text = paste("abline(",ifelse(i > 1, "h", "v"),"=",list[[i]],")")))
    validateCriterion(list[[i]], phaseSpace[,i])
  }
  return(list)
}

subdiv1 <- function(p = 0.3)
{
  testSubdivision(matrix(c(1:10, rep(5, 5), rep(10, 5)), ncol = 2), p = p)
}

subdiv2 <- function(p = 0.3)
{
  testSubdivision(rand(50, 2), p = p)
}

subdiv3 <- function(maxIter = 2)
{
  testSubdivision(rand(50, 2), maxIter = maxIter)
}

subdiv4 <- function(N = 50, maxIter = 2)
{
  testSubdivision(createLogistic(N = N, onlyts = F, tau = 1)$emb, maxIter = maxIter)
  testSubdivision(createLogistic(N = N, onlyts = F, tau = 5)$emb, maxIter = maxIter)
  # testSubdivision(createLogistic(n, onlyts = F, tau = 10)$emb, maxIter = maxIter)
  # testSubdivision(createLogistic(n, onlyts = F, m = 3, tau = 10)$emb, 
    # maxIter = maxIter)  
}

subdiv5 <- function(n = 500, p = 0.2)
{
  testSubdivision(createLogistic(n, onlyts = F, tau = 1)$emb, p = p)
  testSubdivision(createLogistic(n, onlyts = F, tau = 4)$emb, p = p)
}

subdiv6 <- function(p = 0.3)
{
  n = 1000
  testSubdivision(createLorenz(n, onlyts = F, tau = 8)$emb, p = p)
  testSubdivision(createLorenz(n, onlyts = F, tau = 20)$emb, p = p)
}

subdiv7 <- function(p = 0.3)
{
  rand = rand(100, 2)
  testSubdivision(embedd(rand, 2, 3), p = p)
  testSubdivision(embedd(rand, 2, 8), p = p)
}

# ==============================
# All entropies
# ==============================

entropyOptimizationWalk <- function(serie, m, d)
{
  phaseSpace = embedd(serie, m, d)
  plot(phaseSpace)
  cat("The lower this value, the best the estimation:", 
    futureStateEntropy(phaseSpace), "\n")
  cat("the higher this value, the best the estimation:", 
    openBallsEntropy(phaseSpace), "\n")
  #return(phaseSpace)
}
