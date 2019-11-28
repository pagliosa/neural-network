# ===============================
# About: Multidimensional Grid
# Dependences: none
# Author: Rodrigo Mello, Lucas Pagliosa
# Last revision: 03/11/17
# ==============================

require("Canon")

grid <- 0

mdquantize <- function(data, nbins, plot = F)
{
  ndim = ncol(data) # Number of dimensions
  nofs = nrow(data) # Number of states
	cells = NULL
  c = nbins - 1
	
  if (plot)
  {
    if (ndim > 2)
      plot3d(data[,1:3])
    else
      plot(data)
  }			
	for (i in 1:ndim)
	{
	  min = min(data[,i])
	  max = max(data[,i])
	  seq = seq(min, max, (max - min) / c)
	  cells = cbind(cells, findInterval(data[,i], seq, left.open = T, all.inside = T))
	  if (plot && ndim < 3)
      eval(parse(text = paste("abline(",ifelse(i > 1, "h", "v"),"=",seq,")")))
  }

	grid <<- list()
	classBinCount = list()
	dataBinCount = list()
	
  for (i in 1:nofs)
  {
    # Number of points
		key = strcat(as.character(cells[i,]), collapse = "#")
		if (is.null(grid$cells[[key]]))
			grid$cells[[key]] <<- list("nop" = 1)
		else
			grid$cells[[key]]$nop <<- grid$cells[[key]]$nop + 1
		
		# Class bin count
    binKey = strcat(as.character(cells[i, 1:(ndim - 1)]), collapse = "#")
	  grid$cells[[key]]$classBin <<- binKey
	  if (is.null(classBinCount[[binKey]]))
			classBinCount[[binKey]] = 1
		else
			classBinCount[[binKey]] = classBinCount[[binKey]] + 1
	  
	  # Data bin count
    binKey = strcat(as.character(cells[i, ndim]), collapse = "#")
	  grid$cells[[key]]$dataBin <<- binKey
	  if (is.null(dataBinCount[[binKey]]))
			dataBinCount[[binKey]] = 1
		else
			dataBinCount[[binKey]] = dataBinCount[[binKey]] + 1
  }

	grid$classBinCount <<- classBinCount
	grid$dataBinCount <<- dataBinCount
	grid$total <<- nofs

	return (grid)
}

f <- function(grid, coordinate)
{
	key = strcat(as.character(coordinate), collapse = "#")
	if (is.null(grid$cells[[key]]))
	  return (0)
	return (grid$cells[[key]]$nop)
}

getCount <- function(count, bin)
{
	bin = strcat(as.character(bin), collapse = "#")
	if (is.null(count[[bin]]))
	  return (0)
	return (count[[bin]])
}

meanProb <- function(grid, bin = "class", debug = F)
{
  sum = 0
  
  for (i in 1:length(grid$cells))
  {
    cell = grid$cells[[i]]
    if (bin == "class")
    {
      prob = cell$classProb
      count = grid$classBinCount[[cell$classBin]]
      if (debug)
        printf("(C) ")
    }
    else
    {
      prob = cell$dataProb
      count = grid$dataBinCount[[cell$dataBin]]
      if (debug)
        printf("(D) ")
    }
    sum = sum + prob
    if (debug)
      printf("Cell %s: %d / %d (%f)\n", names(grid$cells[i]),
      cell$nop, count, prob)
  }
  return (sum / length(grid$cells))
}

bin <- 0

traverseGrid <- function(coords, sum, level, nbins, rec = F)
{
  if (!rec)
  {
  	for (i in 1:length(grid$cells))
  	{
  	  cell = grid$cells[[i]]
  	  nop = cell$nop
      classProb = nop / grid$classBinCount[[cell$classBin]]
      dataProb = nop / grid$dataBinCount[[cell$dataBin]]
      grid$cells[[i]]$classProb <<- classProb
      grid$cells[[i]]$dataProb <<- dataProb
  	}
    return (0)
  }
  
  ret = 0
  len = length(coords)
  if (level > len)
  {
    nop = f(grid, coords)
    classBinCount = getCount(grid$classBinCount, coords[1:(len - 1)])
    dataBinCount = getCount(grid$dataBinCount, coords[len])
    # cat(coords, nop, classBinCount, "\n")
    key = strcat(as.character(coords), collapse = "#")
    if (nop != 0)
    {
      classProb = nop /classBinCount
      dataProb = nop / dataBinCount
      grid$cells[[key]]$classProb <<- classProb
      grid$cells[[key]]$dataProb <<- dataProb
    }
    return (nop)
  }
  else
  {
    for (i in 1:nbins)
    {
      ret = ret + traverseGrid(coords, sum, level + 1, nbins, rec)
      coords[level] = coords[level] + 1
    }
    # if (level == 2)
    # {
    #   printf("Bin: %d %d\n", bin, ret)
    #   bin <<- bin + 1
    # }
  }
  return (ret)
}

plotResults <- function(ret, name = "Result", ylab = "Probs")
{
  min = min(ret)
  max = max(ret)
  
  plot(1, main = name, type = "n", xlim = c(1, ncol(ret)), 
      ylim = c(min, max), xlab = "Delay", ylab = ylab)
  color = 1
  for (i in seq(1, nrow(ret), 2))
  {
    lines(ret[i,], type = "b", col = color, lty = 2, pch = 20) # circle
    lines(ret[i + 1,], type = "b", col = color, lty = 1, pch = 17) # triangle
    color = color + 1
  }
}

computeJoinProbability <- function(series, N = 200, maxm = 6,
  maxd = 60, nbins = 10, name = "Result", probType = "data", save = T)
{
  nop = length(series)
  ret = zeros(2 * (maxm - 1), maxd)
  count = 1
  
	for (m in 2:maxm)
  {
	  probs = c()
	  vns = c()
		for (d in 1:maxd)
	  {
		  if (N > 1)
      nop = N + (m - 1) * d
			dataset = embedd(series[1:nop], m, d)
			grid <<- mdquantize(dataset, nbins, F)
			bin <<- 0
			
			traverseGrid(zeros(m), 0, 1, nbins + 1, F)
			
			printf("m = %d, d = %d\n", m, d)
			# printf("Number of phase states: %d\n", nrow(dataset))
			# printf("number of non-empty cells: %d\n", length(grid$cells))
			
			prob = meanProb(grid, probType)
			printf("mean probability: %f\n", prob)
			probs = c(probs, prob) 
			vns = c(vns, vonNeumannEntropy(dataset, 1))
		}
 		ret[count,] = probs
 		ret[count + 1,] = vns
 		count = count + 2
	}
  myPlot(function() {plotResults(ret, name = name)}, savePlot = save, 
    filename = sprintf("%sProb%s", probType, name))
  return(probs)
}

runPipeline <- function()
{
  computeJoinProbability(createLogistic(), name = "Logistic", probType = "class")
  computeJoinProbability(createLorenz(), name = "Lorenz", probType = "class")
  computeJoinProbability(createRossler(), name = "Rossler", probType = "class")  
  computeJoinProbability(createLogistic(), name = "Logistic")
  computeJoinProbability(createLorenz(), name = "Lorenz")
  computeJoinProbability(createRossler(), name = "Rossler")
}
