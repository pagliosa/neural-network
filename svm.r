# ===============================
# About: Support Vector Machines
# Dependences: Canon pkg, LowRankQP
# Author: Rodrigo Mello, Lucas Pagliosa
# Last revision: 28/02/18
# ==============================

require("Canon")

loadPackages("LowRankQP", "e1071")

genData <- function(max = 10, n = 30)
{
  y = normalizeVec(cumsum(abs(rnorm(30, sd = 10))), max, 0)
  data = cbind(1:n, y)
  
  plot(data)
  
  colnames(data) = c("X", "Y")
  write.table(numDigits(data, 2), file = "SVM.csv", sep=",", col.names = T, 
    row.names = F)
}

readData <- function()
{
  data = read.csv("SVM.csv", header = T)
  head(data)

  #Scatter Plot
  plot(data, main ="Scatter Plot", type = "b")
  return(data)
}

linearFit <- function(data)
{
  #Fit linear model using OLS
  model = lm(Y ~ X, data)
  
  #Predict Y using Linear Model
  predY <- predict(model, data)
  
  #Overlay best-fit line on scatter plot
  points(data$X, predY, col = "blue", pch = 16, type = "b")
}

svmFit <- function(data)
{
  #Regression with SVM
  modelsvm = svm(Y ~ X, data)

  #Predict using SVM regression
  predYsvm = predict(modelsvm, data)
  
  points(data$X, predYsvm, col = "red", pch = 16, type = "b")
}

pipe <- function()
{
  genData()
  data = readData()  
  # linearFit(data)
  svmFit(data)
}

mySvm <- function(X, Y, C = Inf, order = 1, threshold = 1e-8)
{
	# LowRankQP
	# =========
 	#	min        d^T alpha + 1/2 alpha^T Q alpha
  #	such that  A alpha = b
  #          	 0 <= alpha <= u
	# SVM / Dual
	# ==========
	#	max	e^T alpha - 1/2 alpha^T Q alpha
	#	such that  y^T alpha = 0
	#		         0 <= alpha <= C
	# Adapting the objective function:
	#	max	e^T alpha - 1/2 alpha^T Q alpha
  #                ||
	#	min	-e^T alpha + 1/2 alpha^T Q alpha
	#

	Qmat <- (Y %*% t(Y)) * (1+(X %*% t(X)))^order
	dvec <- rep(-1, nrow(X))
	Amat <- t(Y)
	bvec <- 0
	uvec <- rep(C, nrow(X)) # C is set 1000, should be theoretically inf

	res <- LowRankQP(Qmat, dvec, Amat, bvec, uvec, method = "CHOL")
	alphas <- res$alpha

	# Retrieving the support vectors and alphas
	supportVectors <- which(alphas < C - threshold & alphas > threshold)
	supportAlphas <- alphas[supportVectors]

	# computing the b value
	b <- Y[supportVectors] - t(supportAlphas * Y[supportVectors]) %*% (1+(X[supportVectors,] %*% t(X[supportVectors,])))^order

	return (list(X = X, Y = Y, order = order, supportVectors = supportVectors,
	  supportAlphas = supportAlphas, b = mean(b)))
}
