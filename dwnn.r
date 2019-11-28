# ===============================
# About: Radial Basis Function Distance-Weighted Nearest Neighbors
# implementation
# Dependences: utils.r
# Author: Rodrigo Mello, Lucas Pagliosa
# Last revision 16/10/15
# ==============================

loadPackages("pdist", "RTisean", "tseriesChaos")

sqtDist <- function(a, b)
{
  return(sum((a - b)^2))
}

euclidean <- function(a, b)
{
	return(sqrt(sqtDist(a, b)))
}

gaussianKernel <- function(dist, sigma)
{
  return(exp(-dist^2 / (2 * sigma^2)))
}

computeError <- function(realClass, estimatedClass)
{
  meanSquaredError = sqtDist(realClass, estimatedClass) / length(realClass)
  return(list(meanError = sqrt(meanSquaredError), 
    meanSquaredError = meanSquaredError))
}

# ===============================
# trainSet: training set
# testSet: test set
# sigma: radial kernel radius
# threshold: minimun influence needed to be classified as relevant neighbor
# getnon: should this function also return the number of relevant neighbors?
# ===============================
dwnn <- function(trainSet, testSet, sigma = 1, threshold = 1e-5, getnon = F)
{
  classAttrId = ncol(trainSet)
  nAttrs = classAttrId - 1

  neighbors = c()
  comp = c()
  
	for (i in 1:nrow(testSet))
	{
	  estimatedClass = 0
	  totalWeight = 0
	  weights = c()
	  testInstance = testSet[i, 1:nAttrs]
	  
	  for (j in 1:nrow(trainSet))
    {
  		dist = euclidean(trainInstance, trainSet[j, 1:nAttrs])
  		weight = gaussianKernel(dist, sigma)
  		weights = c(weights, weight)
  		estimatedClass = estimatedClass + weight * trainSet[i, classAttrId]
  		totalWeight = totalWeight + weight
	  }
	  realClass = testSet[j, classAttrId]
	  comp = rbind(comp, cbind(realClass, estimatedClass / totalWeight))
	  neighbors = c(neighbors, sum(weights > threshold))
	}
  error = computeError(comp[,1], comp[,2])$meanError
	if (getnon)
	  return(list(error = error, numberOfRelevantNeighbors = neighbors))
	return(error)
}

dwnnOptimum2 <- function(base, testSet, sigma = 1) {
  
  nAttrs = ncol(base)-1
  class = ncol(base)
  
  base = as.matrix(base)
  testSet = as.matrix(testSet)
  
  euclidean = matrix(pdist(matrix(testSet[,1:nAttrs], ncol=nAttrs), 
    matrix(base[,1:nAttrs], ncol=nAttrs))@dist, 
    nrow=nrow(testSet), byrow=T)
  weights = gaussianKernel(euclidean, sigma)
  
  num = apply(weights, 1, function(x) { x %*% base[,class] })
  
  #obtained = c()
  #for (i in 1:nrow(weights)) {
  #	obtained = c(obtained, (weights[i,] %*% base[,class]) / sum(weights[i,]))
  #}
  
  den = rowSums(weights)
  obtained = num / den
  
  ret = list()
  ret$obtained = obtained
  ret$absError = testSet[,class] - obtained
  ret$absError = abs(testSet[,class] - obtained) # generalization
  ret$expected = testSet[,class]
  return(ret)
}

# Schefler, Statistics...
dwnnOptimum <- function(base, testSet, sigma = 1, relevantDistance = NULL)
{
  nAttrs = ncol(base) - 1
  class = ncol(base)
  
  base = matrix(base, ncol = class)
  testSet = matrix(testSet, ncol = class)
  
  euclidean = matrix(pdist(matrix(testSet[,1:nAttrs], ncol=nAttrs), 
    matrix(base[,1:nAttrs], ncol=nAttrs))@dist, 
    nrow=nrow(base), ncol=nrow(testSet))
  
  base = as.matrix(base[,class])
  
  # Relevant Neighbors
  if (!is.null(relevantDistance))
  {
    rnIds = euclidean <= relevantDistance
    euclidean = as.matrix(euclidean[rnIds])
    base = as.matrix(base[rnIds])
  }
  
  nrn = dim(euclidean)[1]
  printf("Number of relevant neighbors: %d\n", nrn)

  weights = gaussianKernel(euclidean, sigma)
  
  num = apply(weights, 2, function(x) { sum(x * base) })
  den = colSums(weights)
  
  obtained = num / den
  
  ret = list()
  ret$obtained = obtained
  ret$error = testSet[,class] - obtained
  ret$absError = abs(testSet[,class] - obtained) # generalization
  ret$expected = testSet[,class]
  ret$meanError = computeError(obtained, testSet[,class])$meanError
  
  return(ret)
}

timeSeriesDWNN <- function(data, sigma = 1, trainSetSize = 0.9, plot = F, random = T) {

  data = as.matrix(data)
  trainIds = NULL
  testIds = NULL
  
  if (random)
  {
	  trainIds = sample(1:nrow(data), size = round(trainSetSize * nrow(data)))
	  testIds = setdiff(1:nrow(data), trainIds)
  }
  else
  {
    len = round(trainSetSize * nrow(data))
    trainIds = 1:len
    testIds = (len + 1):nrow(data)
  }

	trainSet = as.matrix(data[trainIds,])
	testSet = as.matrix(data[testIds,])

	cat("Training set size = ", nrow(trainSet), "\n")
	cat("Test set size = ", nrow(testSet), "\n")

  res = dwnnOptimum(trainSet, testSet, sigma)
  # res2 = dwnnOptimum2(trainSet, testSet, sigma)
  
  cat("Mean error:", res$meanError, "\n")
  # cat("Mean error:", computeError(res2$obtained, res2$expected)$meanError, "\n")
  
  if (plot)
  {
    par(mfrow = c(1, 2))
    plot(res$obtained)
    plot(res$expected, col = "blue")
    par(mfrow = c(1, 1))
  }
  return(res)
}
