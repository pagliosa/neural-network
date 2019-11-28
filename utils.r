# ===============================
# About: Utilities functions
# Dependences: none
# Author: Lucas Pagliosa
# Last revision: 29/01/17
# ==============================

# Evaluate string
seval <- function(string)
{
  return(eval(parse(text = string)))
}

allTrue <- function(vec)
{
  return(length(which(vec)) == length(vec))
}

canonPath <- function()
{
  return("~/Canon/R")
}

reset <- function(dir = "~/Canon/R", env = globalenv())
{
  clearPlots()
  rm(list = ls(all = TRUE, envir = env), envir = env)
  source(sprintf("%s/utils.r", dir))
  loadGraphics()
  loadPackages("FNN")
  sourceDir()
  sourceDir(dir)
  cls()
}

cls <- function()
{
  cat("\014")
}

allFuncs <- function()
{
  return(ls(all = TRUE, envir = globalenv()))
}

# ------------------
# implicit function
# myPlot(function(){...})
myPlot <- function(plotFunc, save = TRUE, filename = "R-image", dir = "./", 
  newPlot = F)
{
  if (newPlot)
    dev.new()
  par(mar = c(4.4, 5, 3, 0.4) + 0.1, cex.lab = 2.5, cex.axis = 1.8,
    cex.main = 1.8)
  plotFunc()
  if (save)
    savePDF(plotFunc, filename, dir)
}

plotInterval <- function(listOfSpaces, margin = 0.2)
{
  nod = ncol(listOfSpaces[[1]])
  ranges = repmat(c(Inf, -Inf), nod, 1)
  margin = c(-margin, margin)

  for (space in listOfSpaces)
  {
    for (d in 1:nod)
    {
      r = range(space[,d])
      if (ranges[d, 1] > r[1])
        ranges[d, 1] = r[1]
      if (ranges[d, 2] < r[2])
        ranges[d, 2] = r[2]
    }
  }
  plot(1, type = "n", xlim = ranges[1,] + margin, ylim = ranges[2,] + margin)
}

savePDF <- function(plotFunc, filename = "R-image",
  pathDir = "C:/Users/pagliosa/Desktop", width = 7, height = 5)
{
  path = sprintf("%s/%s.pdf", pathDir, filename)
  pdf(file = path, width = width, height = height)
  par(mar = c(4.4, 5, 1, 0.4) + 0.1, cex.lab = 2.5, cex.axis = 2.5,
    cex.main = 2.8)
  plotFunc()
  dev.off()
}

getFiles <- function(path = "./", fileFormat = "[.][Rr]$", recursive = T)
{
  fileData = list()
  dirs = list.dirs(path = path, full.names = TRUE, recursive = recursive)
  count = 1

  for (dir in dirs)
  {
    printf("Searching dir: %s\n", dir)

    files = NULL
    allFiles = list.files(dir, pattern = fileFormat, all.files = FALSE)
    nof = length(allFiles)

    for (rFiles in allFiles)
    {
      files = c(files, file.path(dir, rFiles))
      printf("> File found: %s\n", rFiles)
    }

    if (nof > 0)
    {
      fileData[[count]] = list(files = allFiles, name = getPathFile(dir),
        paths = files, nof = length(allFiles))
      count = count + 1
    }
  }

  return(fileData)
}

sourceFiles <- function(...)
{
  files = c(...);

  for (i in 1:length(files))
    suppressWarnings(suppressMessages(source(files[i])))
}

sourceDir <- function(path = "./")
{
  for (rFiles in list.files(path, pattern = "[.][Rr]$"))
  {
    file = file.path(path, rFiles)
    cat("Loading ", file, "\n")
    suppressWarnings(suppressMessages(source(file)))
  }
}

clearPlots <- function()
{
  graphics.off()
  par(mfrow = c(1, 1))
}

loadPackages <- function(...)
{
  packages = c(...)

  for (p in packages)
  {
    if (p %in% rownames(installed.packages()) == FALSE)
      smwp(install.packages(p))
    else
      smwp(require(p, character.only = TRUE))
  }
}

loadGraphics <- function()
{
  loadPackages("ggplot2", "plotrix", "scatterplot3d", "gplots", "lattice", "rgl")
}

error <- function(message)
{
  message(sprintf("ERROR: %s", message))
  return(0)
}

printfCounter = 0

openHeader <- function(...)
{
  printf(...)
  printfCounter <<- printfCounter + 1
}

closeHeader <- function(...)
{
  pfc = printfCounter - 1
  if (pfc >= 0)
    printfCounter <<- pfc
  printf(...)
}

printfInline <- function(...)
{
  invisible(cat(sprintf(...)))
}

printf <- function(...)
{
  for (i in seq(1, 1, length.out = printfCounter))
    cat(">")
  if (printfCounter > 0)
    cat(" ")
  printfInline(...)
}

highlightPrintf <- function(..., borderSpace = 0)
{
  for (i in seq(1, 1, length.out = borderSpace))
    printf("\n")
  printf("================\n")
  printf(...)
  printf("\n")
  printf("================\n")
  for (i in seq(1, 1, length.out = borderSpace))
    printf("\n")
}

printMat <- function(mat, message, onlySize = F)
{
  if (is.vector(mat))
    mat = matrix(mat, nrow = 1)

  if (missing(message))
    message = "mat"

  printf("%s (%dx%d)\n", message, nrow(mat), ncol(mat))

  if (onlySize)
    return()

  for (i in 1:nrow(mat))
  {
    printf("")
    for (j in 1:ncol(mat))
      printfInline("%g ", mat[i, j])
    printfInline("\n")
  }
}

# Create dir paths if dont exist
createDirPath <- function(dirPath)
{
  dirs = strsplit(dirPath, "/")[[1]]
  dir = "./"
  for (i in 1:length(dirs))
  {
    dir = sprintf("%s%s/", dir, dirs[i])
    if (!dir.exists(dir))
    {
      printf("Creating directory %s\n", dir)
      dir.create(dir)
    }
  }
  return(dirPath)
}

# Create file
createPath <- function(dirPath, file, ext = "txt")
{
  createDirPath(dirPath)
  return(sprintf("%s%s.%s", dirPath, file, ext))
}

getPathWithoutExtension <- function(path)
{
  return(strsplit(path, ".", fixed = TRUE)[[1]][1])
}

addList <- function(list, newItem)
{
  list[[length(list) + 1]] = newItem
  return(list)
}

# Suppres messages, warnings and prints
smwp <- function(func, debug = F)
{
  sw <- function(func) suppressMessages(suppressWarnings(func))

  if (debug)
    res = sw(func)
  else
    sw(capture.output({
      res = func;
    }))

  return(res)
}

getOnlyResultsThatWorked <- function()
{
  set.seed(123)
  x <- stats::rnorm(50)
  doit <- function(x)
  {
    x <- sample(x, replace = TRUE)
    if(length(unique(x)) > 30) mean(x)
    else stop("too few unique points")
  }
  ## alternative 1
  res = lapply(1:100, function(i)
  {
    cat(i)
    try(doit(x), TRUE)
  })
  return(res)
}

numDigits <- function(data, number = 4)
{
  data = format(round(data, number), nsmall = number)
  return(apply(data, 2, as.numeric))
}

bnorm <- function(p, mu = 0, sd = 1)
{
  pn = pnorm(p, mu, sd)

  return(2 * pn - 1)
}

normalize <- function(min, max, VALUE, MIN, MAX)
{
  # value - min       VALUE - MIN
  # _____________ =  _____________
  #  max - min         MAX - MIN
  return((VALUE - MIN) * (max - min) / (MAX - MIN) + min)
}

normalizeVec <- function(vec, min = 0, max = 1)
{
  MIN = min(vec)
  MAX = max(vec)

  norm = c()
  for (i in 1:length(vec))
    norm = c(norm, normalize(min, max, vec[i], MIN, MAX))
  return(norm)
}

normalizeSum <- function(vec) {
  s = sum(vec)
  if (s == 0)
    s = 1
  return(vec / s)
}

sqtDist <- function(a, b)
{
  return(sum((a - b)^2))
}

euclidean <- function(a, b, ...)
{
  return(sqrt(sqtDist(a, b)))
}

init <- function(n, m, value, min = 0, max = 1, byrow = T)
{
  if (missing(n))
    n = 1
  if (missing(m))
    m = 1

  # Fill matrix with increasing order
  if (missing(value))
    return(matrix(1:(n * m), n, m, byrow = byrow))

  # Fill matrix with random numbers
  if (is.logical(value) && !value)
    return(matrix(runif(n * m, min, max), n, m))

  # Fill matrix with repeated value
  return(repmat(value, n, m))
}

contains <- function(v1, v2)
{
  for (i in 1:length(v1))
    if (!any(v2 == v1[i]))
      return(FALSE)
  return(TRUE)
}

shannonEntropy <- function(vec, k = length(vec))
{
  sum = sum(vec)
  if (sum == 0)
    return(0)

  probs = vec / sum
  probs[probs == 0] = 1e-20
  probs = probs[1:k]
  shannon = -probs %*% log2(probs)
  return(shannon)
}

computeStatistics <- function(vec)
{
  n = length(vec)

  # Statistical moments
  mean = mean(vec)
  variance = var(vec)

  # Shannon
  shannon = shannonEntropy(vec)

  # Energy
  energy = sum(abs(vec))/n

  return(list("mean" = mean, "var" = variance, "shannon" = shannon,
    "energy" = energy))
}

zeros <- function(n, m = 1)
{
  ret = matrix(0, n, m)
  if (missing(m))
    ret = as.vector(ret)
  return(ret)
}

lastOccurence <- function(string, char)
{
  indices = gregexpr(char, string)[[1]]
  lastOccurence = indices[length(indices)]
  return(lastOccurence)
}

getPathFile <- function(path)
{
  lastSlide = lastOccurence(path, "/") + 1
  file = substr(path, lastSlide, nchar(path))
  return(file)
}

globalClock = proc.time()

startClock <- function()
{
  # <<- sets global variables
  clock = proc.time()
  globalClock <<- clock
  return(clock)
}

endClock <- function(clock)
{
  if (missing(clock))
    clock = globalClock
  diff = proc.time() - clock
  return(diff)
}

getOS <- function()
{
  return(as.vector(Sys.info()["sysname"]))
}

fatorial <- function(arg)
{
  if (arg == 0)
    return(1)

  fat = 1
  for (i in 1:arg)
    fat = fat * i
  return(fat)
}

centralizeEmbddings <- function(emb1, emb2)
{
  for (i in 1:ncol(emb1))
    emb1[,i] = emb1[,i] + mean(emb2[,i]) - mean(emb1[,i])
  return(emb1)
}

computeSlope <- function(xvec, yvec, plot = F)
{
  lm = lm(yvec ~ xvec)
  y = as.numeric(lm$coefficients[1])
  slope = as.numeric(lm$coefficients[2])

  if (plot)
    abline(y, slope)

  return(list("y" = y, "slope" = slope))
}

computeKnnData <- function(space, k)
{
  knn = get.knn(space, k, algorithm = "kd_tree")
  return(list("id" = knn$nn.index, "dist" = knn$nn.dist))
}

differentiate <- function(vec, h = 1e-3)
{
  len = length(vec)
  derivative = (vec[3:len] - vec[1:(len - 2)])/(2 * h)
  return(derivative)
}

lnorm <- function(vec, n = 2)
{
  abs = abs(vec)
  infNorm = max(abs)
  norm = sum(abs^n)^(1/n)

  if (n == Inf || norm == Inf)
    return(infNorm)
  return(norm)
}

getEingValues <- function(X, normalize = T, percentage = 1)
{
  # The singular values of a MxN matrix X are the square roots of the eigenvalues
  # of the NxN matrix X%*%t(X) (where t(X) stands for the transpose-conjugate
  # matrix if it has complex coefficients, or the transpose if it has real
  # coefficients).
  # Thus, if X is NxN real symmetric matrix with non-negative eigenvalues, then
  # eigenvalues and singular values coincide, but it is not generally the case!

  if (ncol(X) == nrow(X))
    S = X
  else {
    S = t(X) %*% X

    # svd = svd(X)$d^2
    # fro = norm(X, "F")^2
    # eig = eig(S)
    # sum(eig) = sum(svd) = fro
  }

  eig = eigen(S)
  values = eig$values
  vectors = eig$vectors

  normalizedValues = values / sum(values)
  p = cumsum(normalizedValues)
  index = 1

  for (i in 1:length(values))
  {
    if (isGreaterEqual(p[i], percentage))
    {
      index = i
      break
    }
  }
  if (normalize)
    ret = normalizedValues[1:index]
  else
    ret = values[1:index]
  vectors = vectors[1:index,]
  return(list("values" = ret, "vectors" = as.matrix(vectors)))
}

isZero <- function(a, eps = 1e-8)
{
  return(abs(a) - eps <= 0)
}

isEqual <- function(a, b)
{
  return(isZero(a - b))
}

isGreater <- function(a, b, eps = 1e-8)
{
  return((a - b) > eps)
}

isGreaterEqual <- function(a, b, eps = 1e-8)
{
  return((a - b) >= -eps)
}

isLesser <- function(a, b, eps = 1e-8)
{
  return(a - b < -eps);
}

isLesserEqual <- function(a, b, eps = 1e-8)
{
  return(a - b < eps);
}


heaviside <- function(v, threshold = 0.5)
{
  v[v > threshold] = 1
  v[v <= threshold] = 0
  return(v)
}

minor <- function(A, i, j)
{
  return(det(A[-i,-j]))
}

cofactor <- function(A, i, j)
{
  return((-1)^(i+j) * minor(A,i,j))
}

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

plotIris <- function()
{
  pairs(iris[1:4], main = "Edgar Anderson's Iris Data", pch = 21,
    bg = c("red", "green3", "blue")[unclass(iris$Species)])
}

standarlize <- function(X = matrix(1:10, ncol = 2))
{
  mean = apply(X, 2, mean)
  std = apply(X, 2, std)
  tmp = rbind(1:ncol(X), X)
  s = apply(tmp, 2, function(c)
    {
      id = c[1]
      c = (c - mean[id]) / std[id]
    })
  s1 = s[-1,]

  n = nrow(X)
  one = rep(1, n)
  meanx = drop(one %*% X) / n
  s2 = scale(X, mean, F)
  stdx = sqrt(drop(one %*% (s2^2)))
  s2 = scale(s2, F, std)

  ret = list("hand" = s1, "step" = s2, "R" = scale(X, T, T))
  return(ret$R)
}

bibtex <- function(packageName = "gplots")
{
  toBibtex(citation(packageName))
}

findMaximum <- function(vec, k = 1)
{
  v = vec[1]
  count = 0

  for (i in 2:length(vec))
  {
    if (v > vec[i])
      count = count + 1
    if (count == k)
      return(i - 1)
    v = vec[i]
  }
  return(1)
}

parseAnscombe <- function(plot = T)
{
  data = with(anscombe, data.frame(x = c(x1, x2, x3, x4), y = c(y1, y2, y3, y4),
    group = gl(4, nrow(anscombe))))

  printf("Mean:\n")
  print(as.matrix(aggregate(.~group, data = data, mean)), quote = F)

  printf("Variance:\n")
  print(as.matrix(aggregate(.~group, data = data, var)), quote = F)

  if (plot)
  {
    coefs = lm(y ~ x, data = data)$coefficients
    pf = ggplot(data, aes(x = x, y = y, color = "orange")) + geom_point() +
        facet_wrap(~group) + theme_bw() + theme(legend.position = "none") +
        geom_abline(aes(intercept = coefs[1], slope = coefs[2]))
    myPlot(function() {print(pf)}, T, "Anscombe's quartet")
  }

  return(data)
}

seeKernel <- function(n = 100, order = 2)
{
  source("R/timeSeries.r")
  inputSpace = cbind(createSine(N = n, onlyts = F, tau = 10)$emb)
  inputSpace = rbind(inputSpace, cbind(rnorm(n, 0, 0.1), rnorm(n, 0, 0.1)))
  featureSpace = (inputSpace %*% t(inputSpace))^order
  classes = rbind(c(rep(1, n), rep(2, n)))

  plot(inputSpace, col = classes)
  plot3d(featureSpace, col = classes)
}

MImatrix <- function(A, b, nbins = 4)
{
  if (missing(A))
  {
    A = matrix(rnorm(6 * 3), ncol = 2)
    A = cbind(A, A[,1])
  }
  n = ncol(A)
  mim = eye(n)
  c = nbins - 1
  bins = NULL

  D = cbind(A, b)
  d = n + 1
  for (i in 1:d)
  {
    v = D[,i]
    min = min(v)
	  max = max(v)
    seq = seq(min, max, (max - min) / c)
    bins = cbind(bins, findInterval(v, seq, left.open = T,
      all.inside = T))
  }
  for (i in 1:n)
    for (j in (i:n))
      mim[i, j] = mim[j, i] = mutinformation(bins[,i], bins[,j])
  class = zeros(n)
  for (i in 1:n)
    class[i] = mutinformation(bins[,i], bins[,d])
  return(list(x = mim, y = class))
}

myNewtonRaphson <- function(n = 10, p = 3)
{
  X = matrix(rnorm(n * p), ncol = p)
  y = X[,1] + 0.566 + X[,3] + 0.432

  lambda = 2
  beta = as.matrix(rep(1, p))

  f = 1
  iter = 0

  while (!all(isZero(f)) & iter < 10)
  {
    cat(f, "\n")
    iter = iter + 1
    f = t(X) %*% (y - X %*% beta) - sign(beta) * lambda
    J = t(X) %*% X
    beta = beta + solve(J, -f)
    printf("Iter: %d\n", iter)
  }
}

printPascal <- function(n = 5)
{
  for (line in 1:n)
  {
    C = 1;  # used to represent C(line, i)
    for (i in 1:line)
    {
      printf("%d ", C)  # The first value in a line is always 1
      C = C * (line - i) / i;
    }
    printf("\n");
  }
}

weight <- function(dist, sigma) {
	return(exp(-dist^2 / (2 * sigma^2)))
}

dwnn <- function(dataset, query, sigma = 0.5) {

	classId = ncol(dataset)

	dist =  apply(dataset, 1, function(row) {
	  sqrt(sum((row[1:(classId - 1)] - query)^2))
  })

	w = weight(dist, sigma)
	Y = dataset[,classId]

	y = sum(w * Y) / sum(w)

	return(y)
}

dwnnManager <- function(dataset, timeDelay, trainFrom, trainTo, sigmas, window,
  plot = T, save = F, name = "DNWW", recursive = TRUE)
{
  # Put here to nice cleaner
  if (missing(trainFrom))
    trainFrom = 1
  if (missing(trainTo))
    trainTo = floor(nrow(dataset) * 0.75)
  if (missing(sigmas))
    sigmas = computeSigmas(dataset, window)

	originalDataset = dataset
	ncol = ncol(dataset)

	pclass = c()

	testRange = (trainTo + 1):nrow(dataset)
	for (to in testRange)
	{
		query = as.numeric(dataset[(to - timeDelay), 2:ncol])
		y = dwnn(dataset[trainFrom:(to - 1),], query, sigmas[to])

		if (recursive)
		  dataset[to,] = c(query, y)
		pclass = c(pclass, y)
	}

	oclass = originalDataset[testRange, ncol]

	ret = list()
	ret$originalDataset = originalDataset
  ret$predictedDataset = dataset
	ret$originalClass = oclass
	ret$predictedClass = pclass
  ret$error = euclidean(oclass, pclass) / length(pclass)

  if (plot)
  {
    myPlot(function() {
      plot(oclass, type = "l", ylab = "Value", xlab = "Time",
      	      lty = 2, main = sprintf("Error: %g", ret$error))
      lines(pclass, col = "blue")
      par(mfrow = c(2, 1))
      if (ncol == 3)
      {
        scatterplot3d(originalDataset)
        scatterplot3d(dataset)
      }
      else
      {
        range = (ncol - 1):ncol
        plot(originalDataset[,range])
        plot(dataset[,range])
      }
      par(mfrow = c(1, 1))
    }, save, name)
  }

	return(ret)
}

computeSigmas <- function(data, windowLength = nrow(data) / 3)
{
  k = floor(log(windowLength)) + 1

  # Computing distance from each state to all states
  distMatrix = as.matrix(dist(data, upper = T, diag = T))

  # Sorting each state distances
  distMatrixSortedByRow = t(apply(distMatrix, 1, sort))

  # Getting the distances for each kth closest neighbor
  cut = distMatrixSortedByRow[,k]

  return(cut)
}
