# ===============================
# About: Multilayer Perceptron and Estimation of Embedding Parameters
# Dependences: utils.r and timeSeries.r
# Author: Rodrigo Mello, Lucas Pagliosa
# Date of creation: 19/09/18
# Last revision: 19/09/18
# ==============================
#
# This nomenclature below goes in accordance with Mello's class about neural
# network
#
# Output layer
# ============
#
# w_kj(t+1) = w_kj(t) + stepSize * (y_pk - o_pk) * o_pk * (1-o_pk) * i_pj
# thstepSize_k(t+1) = thstepSize_k(t) + stepSize * (y_pk - o_pk) * o_pk * (1-o_pk) * 1
#
# Hidden layer
# ============
#
# w_ji(t+1) = w_ji(t) + stepSize *
#	 i_pj * (1-i_pj) * x_pi * sum_k (y_pk - o_pk) * o_pk * (1-o_pk) * w^o_kj
#
# thstepSize_j(t+1) = thstepSize_j(t) + stepSize *
#	 i_pj * (1-i_pj) * 1 * sum_k (y_pk - o_pk) * o_pk * (1-o_pk) * w^o_kj
#
# ============
#
# Knowing that:
#
# x_p is the pth data input
# i_pj is the pth exit from the jth hidden neuron
# o_pk is the pth exit from the kth output neuron
#
# Also:
#
# N is the number of input neurons
# L is the number of hidden neurons
# M is the number of output neurons

# Put below the dir path of where you stored those files
source("utils.r")
source("timeSeries.r")

library(ggplot2)

test <- function(nop = 10)
{
  a = matrix(runif(nop), ncol = 2)

  b = impactFactor(a, F)

  # Plot the barplot, save the x-axis position of each bar in n
  n <- barplot(b, ylim = c(0, max(b) + 10))

  # Plot the boxplot using add = TRUE  and the "at" argument for the x-axis
  # positions which are saved in n
  bp = boxplot(a, plot = F)
  bp$stats = bp$stats + t(replicate(5, b))
  bxp(bp, at = n, add = TRUE, axes = FALSE)
}

# Activation function: sigmoid
fnet <- function(net)
{
  return (1 / (1 + exp(-net)))
}

# Partial of fnet: sigmoid
partialNet <- function(net)
{
  return (net * (1 - net))
}

# Net architecture
architecture <- function(N, L, M)
{
  # Each neuron has its connections and a bias, therefore it is added one in the
  # sizes below
  hidden = matrix(runif(min = -0.1, max = 0.1, n = (N + 1) * L), nrow = L)
  output = matrix(runif(min = -0.1, max = 0.1, n = (L + 1) * M), nrow = M)

  # nh = (N + 1) * L
  # t = ones(nh, 1)
  # t[1:(nh / 2)] = 0
  # hidden = matrix(t, nrow = L)
  #
  # no = (L + 1) * M
  # t = ones(no, 1)
  # t[1:(no / 2)] = 0
  # output = matrix(t, nrow = M)

  ret = list()
	ret$N = N
	ret$L = L
	ret$M = M
	ret$hidden = hidden
	ret$output = output
  ret$lastHiddenVariation = zeros(L, N + 1)
  ret$lastOutputVariation = zeros(M, L + 1)

	return (ret)
}

# Forward execution
forward <- function(arch, x_p)
{
	hidden_p = arch$hidden %*% c(x_p, 1)
	h_p = fnet(hidden_p)

	output_p = arch$output %*% c(h_p, 1)
	o_p = fnet(output_p)

	ret = list()
	ret$h_p = h_p
	ret$o_p = o_p

	return (ret)
}

# Forward all instances
forwardAll <- function(model, x)
{
  o = c()
	for (p in 1:nrow(x)) {
		x_p = x[p,]

		run = forward(model, x_p)
		o = c(o, run$o_p)
	}
	o[o < 1e-3] = 0
	o[o > 0.99] = 1
	return (o)
}

# Training network
# arch: neural network architecture
# dataset: is the matrix of instances in the form of XY, that is, instances and
# classes
# stepSize: step size in gradient descent method, a.k.a. learning rate
# momentum: percentage betweein 0 and 1. To avoid local minima, this heuristic
# includes a factor of previous variation in the current iteration. A.k.a.
# momentum
# tol: error tolerance threshold
# epochs: maximum number of iterations allowed
# lwf: learning with forgetting
# forgetFactor: forgetting factor. If equals to 0, do nothing
# hucFactor: hidden unit clarification factor. If equals to 0, do nothing
# selectiveFactor: selective learning threshold. If equals to 1, do notinhg
train <- function(arch, dataset, stepSize = 0.1, momentum = 0.2, tol = 1e-3,
  epochs = 1e5, forgetFactor = 0, hucFactor = 0, selectiveThreshold = Inf)
{
  # Making sure data is a matrix
  dataset = as.matrix(dataset)

  # Make it easier
  N = arch$N
  L = arch$L

  # Predefining ranges
	inputRange = 1:N
	classRange = (N + 1):ncol(dataset)

	# Number of elements
	numberOfInstances = nrow(dataset)

	# Mean squared error
  mse = 2 * tol

  # Number of epochs counter
  epoch = 0

	while (mse > tol && epoch < epochs)
  {
    mse = 0
    epoch = epoch + 1

    # Simulating some data picking randomness
    order = sample(1:numberOfInstances)

    for (r in 1:numberOfInstances)
    {
	    p = order[r]
      x_p = dataset[p, inputRange]
		  y_p = dataset[p, classRange]

		  run = forward(arch, x_p)
		  error = y_p - run$o_p
		  mse = mse + sum(error^2)

  		# Computing deltas
  		# ========

		  # Delta output: delta^o_pk = - (y_pk - o_pk) * (o_pk)'
		  d_o = - error * partialNet(run$o_p)

		  # Delta hidden: delta^h_pj = (i_pj)' * sum_k delta_pk * w^o_kj
		  sum_k = c()
		  for (j in 1:L)
		    sum_k = c(sum_k, sum(d_o * arch$output[,j]))
		  d_h = partialNet(run$h_p) * sum_k

		  # Current output (hidden) direction: partial and momentum
		  od = d_o %*% c(run$h_p, 1) + momentum * arch$lastOutputVariation
		  hd = d_h %*% c(x_p, 1) + momentum * arch$lastHiddenVariation

		  # Learning with Forgetting
		  olwf = forgetFactor * sign(od)
		  hlwf = forgetFactor * sign(hd)

		  # Selective Forgetting
	    olwf = olwf[abs(od) > selectiveThreshold] = 0
	    hlwf = hlwf[abs(hd) > selectiveThreshold] = 0

  	  # Hidden Unit Clarification
	    s = rep(1, L)
	    s[run$h_p > 0.5] = -1

	    # Adding no HUC bias
	    inc = cbind(matrix(rep(x_p, L), ncol = N, byrow = T), rep(0, L))
	    huc = hucFactor * s * inc

	    # Selective Learning with Forgetting
	    od = od + olwf
	    hd = hd + hlwf + huc

		  # Variations: step size times direction
		  ov = stepSize * od
		  hv = stepSize * hd

		  # Updating weights and bias according to the gradient descent
		  # ========
		  arch$output = arch$output - ov
      arch$hidden = arch$hidden - hv

      # Updating variation
      arch$lastOutputVariation = ov
      arch$lastHiddenVariation = hv
    }

    mse = mse / numberOfInstances

    cat(sprintf("(%d/%d) MSE: %g\n", epoch, epochs, mse))
	}

	return (arch)
}

# Structural Learning with Forgetting
slf <- function(arch, dataset, stepSize = 0.1, momentum = 0.2, tol = 1e-3,
  epochs = 1e5, forgetFactor = 1e-5, hucFactor = 1e-3, selectiveThreshold = 0.1,
  debug = FALSE)
{
  sampleSize = 0.5
  id = sample(1:nrow(dataset), size = nrow(dataset) * sampleSize)
	set1 = dataset[id,]
	set2 = dataset[-id,]

  lwf = train(arch, set1, stepSize, momentum, tol, epochs, forgetFactor, 0, Inf)
  huc = train(lwf, set2, stepSize, momentum, tol, epochs, forgetFactor,
    hucFactor, selectiveThreshold)

  huc$impact = impactFactor(huc$hidden)

  # Debugguing
  if (debug)
  {
    nod = ncol(dataset)
    x = dataset[, 1:(nod - 1)]
    y = dataset[, nod]
    o1 = forwardAll(lwf, x)
    o2 = forwardAll(huc, x)

    d1 = sum(abs(y - o1))
    d2 = sum(abs(y - o2))

    printMat(cbind(y, o1, o2), "Comparison")
    cat(c(d1, d2), "\n")
  }

  return (huc)
}

# Comput impact factor of the hidden layer
impactFactor <- function(weights, hasBias = TRUE)
{
  if (!is.matrix(weights))
    weights = as.matrix(weights)
  w = abs(weights)
  if (nrow(weights) > 1)
    w = colSums(w)
  if (hasBias)
    w = w[1:(length(w) - 1)]
  return (w)
}

# Visualize adapted hidden layer
plotImpact <- function(impact, name = "Weights", ut = 0.8, lt = 0.1, ylim)
{
  if (missing(ylim))
    ylim = c(0, max(impact))
  n <- barplot(impact, main = name, ylim = ylim, names.arg = 1:length(impact),
    xlab = "Input Neurons", "ylab" = "Input Importance")

  # Stats
  max = max(impact)

  # Plot quantile
  abline(h = quantile(impact, ut), lty = 2, col = "blue")

  # Plot max-based
  abline(h = max * lt, lty = 1, col = "red")

  return (n)
}

booleanTest <- function(stepSize = 0.1, momentum = 0.2, tol = 1e-3,
  epochs = 5e2, forgetFactor = 1e-4, hucFactor = 1e-3, selectiveThreshold = 0.1,
  debug = FALSE)
{
  f <- function(a, b, c, d, e, f, g, h)
  {
    return (as.numeric((a || b) && (c || e)))
  }

  fapply <- function(x)
  {
    return (f(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]))
  }

  x = as.matrix(expand.grid(replicate(8, 0:1, simplify = FALSE)))
  y = apply(x, 1, fapply)
  dataset = cbind(x, y)
  inputNames = c("a", "b", "c", "d", "e", "f", "g", "h")
  colnames(dataset) = c(inputNames, "o")

	# Creating and training neural network
	arch = architecture(ncol(dataset) - 1, 4, 1)
	model = slf(arch, dataset, stepSize, momentum, tol, epochs, forgetFactor,
	  hucFactor, selectiveThreshold, debug)

	# Return list
	ret = list()
  ret$inputNames = inputNames
  ret$ret = cbind(y, forwardAll(model, x))

  return (ret)
}

xorTest <- function(stepSize = 0.1, tol = 1e-2)
{
	arch = architecture(2, 2, 1)
	model = train(arch, read.table("./datasets/xor.dat"), stepSize, tol = tol)

	cat("Hidden layer...\n")
	print(model$hidden)

	x = seq(0, 1, length = 100)

	w_0 = model$hidden[1,] # neuron 0
	w_1 = model$hidden[2,] # neuron 1

	net_h_0 = outer(x, x, function(x_0, x_1) { cbind(x_0, x_1, 1) %*% w_0 } )
	net_h_1 = outer(x, x, function(x_0, x_1) { cbind(x_0, x_1, 1) %*% w_1 } )

	id = which(net_h_0 > 0)
	net_h_0[id] = 1
	net_h_0[-id] = 0

	id = which(net_h_1 > 0)
	net_h_1[id] = 1
	net_h_1[-id] = 0

	result = net_h_0 + net_h_1
	filled.contour(x, x, result, nlevels = 5)
}

irisTest <- function(sampleSize = 0.5, hiddenSize = 3, stepSize = 0.2,
  momentum = 0.2, tol = 1e-2)
{
  # Reading and adjusting dataset
	dataset = read.table("./datasets/iris.dat")
	dataset = as.matrix(dataset)

	# Creating training and test sets
	id = sample(1:nrow(dataset), size = nrow(dataset) * sampleSize)

	trainingSample = dataset[id,]
	testSample = dataset[-id,]

	# Creating and training neural network
	arch = architecture(4, hiddenSize, 3)
	model = train(arch, trainingSample, stepSize, momentum, tol)

	# Range variables
	inputRange = 1:arch$N
	classRange = (arch$N + 1):ncol(dataset)

	# Training...
	counter = 0
	for (p in 1:nrow(trainingSample)) {
		x_p = trainingSample[p, inputRange]
		y_p = trainingSample[p, classRange]

		run = forward(model, x_p)

		cat(y_p, "", run$o_p, "\n")

		resp = run$o_p
		resp[run$o_p > 0.5] = 1
		resp[run$o_p <= 0.5] = 0

		if (sum((y_p - resp)^2) == 0)
			counter = counter + 1
	}
	accuracy = counter / nrow(trainingSample)
	R_emp = 1 - accuracy

	# Testing...
	counter = 0
	for (p in 1:nrow(testSample))
  {
		x_p = testSample[p, inputRange]
		y_p = testSample[p, classRange]

		run = forward(model, x_p)

		cat(y_p, "", run$o_p, "\n")

		resp = run$o_p
		resp[run$o_p > 0.5] = 1
		resp[run$o_p <= 0.5] = 0

		if (sum((y_p - resp)^2) == 0)
			counter = counter + 1
	}
	accuracy = counter / nrow(testSample)
	R_test = 1 - accuracy

	cat("R_emp: ", R_emp, " R_test: ", R_test, "\n")
}

# Provisorial Embedding Vector format
pevFormat <- function(series, m, tau)
{
  c = (m - 1) * tau
  n = length(series) - c
  embed = zeros(n, c + 1)
  for (i in 1:n)
    embed[i,] = series[i:(i + c)]
  return (embed)
}

# Time series  neural network
tsTest <- function(ts = createHenon(), maxM = 6, maxTau = 3, sampleSize = 0.75,
  K = 5, stepSize = 0.1, momentum = 0.2, tol = 1e-3, epochs = 5e2,
  forgetFactor = 1e-3, hucFactor = 0, selectiveThreshold = Inf, uppert = 0.8,
  lowert = 0.1, test = TRUE, plot = TRUE, dir = "Henon", save = TRUE,
  recursive = FALSE, random = T)
{
  # Check for list
  if (is.list(ts))
    ts = ts$ts

  # Putting extra information on directory name
  dir = sprintf("Results/%s (%d,%d,%g)", dir, maxM, maxTau, uppert)

  # Create dir if does not exist
  if (save)
    createDirPath(dir)

  # Normalize series
  ts = normalize(0, 1, ts, min(ts), max(ts))

  # Creating provisional phase space
  dataset = pevFormat(ts, maxM, maxTau)

  # Number of dimensions
  nod = ncol(dataset)

  # Cross Validation Impact
  cvImpact = c()

  # Number of training samples
  nos = nrow(dataset)
  trainTo = floor(nos * sampleSize)
  testSize = nos - trainTo - 1
  
  # Cross validation parameters
  for (k in 1:K)
  {
    startTime = Sys.time()

    foldDataset = dataset
    if (recursive)
  	  trainingSample = dataset[1:trainTo,]
    else
    {
      if (random)
      {
        rids = sample(1:nos, trainTo)
        trainingSample = dataset[rids,]
        testSample = dataset[-rids,]
      }
      else
      {
        trainingSample = dataset[1:trainTo,]
        testSample = dataset[(trainTo + 1):nos,]
      }
    }
    
  	# The last column will be the predict value (class)
  	N = nod - 1

  	# The number of hidden layers is set as log(NumberOfImputs) + 1. This is a
  	# based on the SLT to avoid over and underfitting
    L = ceiling(log(N)) + 1

    # Create architecture
  	arch = architecture(N, L, 1)

  	# Training
  	model = slf(arch, trainingSample, stepSize, momentum, tol, epochs,
  	  forgetFactor, hucFactor, selectiveThreshold)

  	endTime = Sys.time()

  	cat(sprintf("Processing time: %g\n", endTime - startTime))

    if (plot)
    {
      name = sprintf("Fold %d (upper = %g, m_max = %d, t_max = %d)", k, uppert,
        maxM, maxTau)
      myPlot(function() {plotImpact(model$impact, name, uppert, lowert)},
        save, name, dir)
    }

  	# Testing recursivelly
  	if (test)
  	{
    	vec.y_p = c()
    	vec.o_p = c()
    	for (p in 1:testSize)
      {
    		# Important: our model does not yield good results due to the butterfly
    		# effect. For a recusive forecasting, use another method on the phase
    		# space reconstructed using our estimation
    		if (recursive)
  		  {
      	  oid = trainTo + 1 + p
      	  
    		  x_p = foldDataset[(oid - maxTau), 1:arch$N]
    		  y_p = dataset[oid, arch$N + 1]
    		  
      		run = forward(model, x_p)
    		  foldDataset[oid,] = c(x_p, run$o_p)
    		}
    	  else
    	  {
    	    x_p = testSample[p, 1:arch$N]
    	    y_p = testSample[p, arch$N + 1]
      		run = forward(model, x_p)
    	  }

    		vec.y_p = c(vec.y_p, y_p)
    		vec.o_p = c(vec.o_p, run$o_p)
    	}

    	error = euclidean(vec.y_p, vec.o_p) / length(vec.y_p)

  	  # Plot forecasted results
    	if (plot)
    	{
    	  myPlot(function() {
          par(mfrow = c(2, 1))
          range = (nod - 1):nod
          plot(dataset[,range])
          plot(foldDataset[,range])
          par(mfrow = c(1, 1))
    	    plot(vec.y_p, type = "l", ylab = "Value", xlab = "Observation",
    	      lty = 2, main = sprintf("Error: %g", error))
	        lines(vec.o_p, col = "blue")
    	  }, save, sprintf("Forecasting fold %d (upper = %g)", k, uppert), dir)
    	}
  	}

  	# Store impact
  	cvImpact = rbind(cvImpact, model$impact)
  }

  # Impact average
	finalImpact = impactFactor(cvImpact, FALSE) / nrow(cvImpact)

	# Plot impact
	box = boxplot(cvImpact, plot = F)
	name = sprintf("Aggregation (upper = %g, m_max = %d, t_max = %d)", uppert,
        maxM, maxTau)
  myPlot(function()
  {
    n = plotImpact(finalImpact, name, uppert, lowert,
      c(0, max(box$stats[5,]) * 1.2))
    bp = boxplot(cvImpact, at = n, add = TRUE, axes = FALSE, col = "lightgrey")
  }, save, name, dir)

  # Return list
  ret = list()

  ret$cvImpact = cvImpact
  ret$finalImpact = finalImpact
  ret$box = box

  return (ret)
}

meb <- function(maxM = 5, maxTau = 5)
{
  return(list(maxM = maxM, maxTau = maxTau))
}

nnPipeline <- function(epochs = 500)
{
  tsVec = list(createRandomSeries, createRandomWalk, createLogistic,
    createHenon, createLorenz, createRossler, createSunspot)
  mebVec = list(meb(6, 6), meb(6, 8), meb(5, 5), meb(4, 3), meb(5, 3),
    meb(8, 8), meb(3, 5))
  for (i in 1:length(tsVec))
  {
    tdata = tsVec[[i]]()
    mdata = mebVec[[i]]
    tsTest(ts = tdata$ts, dir = tdata$name, K = 5, mdata$maxM, mdata$maxTau,
      epochs = epochs, forgetFactor = 1e-3, uppert = 0.8, lowert = 0.1)
  }
}
