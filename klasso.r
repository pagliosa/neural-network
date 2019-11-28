# ===============================
# About: Linear Model Selection and Regularization
# Dependences: Canon pkg
# Author: Rodrigo Mello, Lucas Pagliosa
# Last revision: 10/01/18
# ==============================

require("Canon")

loadPackages("limSolve")

printLog <- function(lambda, beta0, beta)
{
	openHeader("Parameters\n")
  printf("Lambda: %g\n", lambda)
	printMat(as.matrix(c(beta0, beta)), "Betas:")
  closeHeader("Done\n")
}

verifySolution <- function(g, A, lambda, beta0, beta, cVec, oneVec, n)
{
  cond1 = all(isEqual(beta0, t(cVec) %*% oneVec / n))
  if (length(A) == 0)
    return(cond1)
  cond2 = all(isEqual(g[A], sign(beta[A] * lambda)))
  return(cond1 & cond2)
}

klasso <- function(X, y, lambda = 5, eta = 0.01, ikv = 0, order = 1, lambdaThreshold = 1e-3, 
  errorThreshold = 1e-3, debug = F)
{
  # ikv stands for irregular kernel value
	K = (X%*%t(X) + ikv)^order

	# System dimensios
	n = nrow(K)
	p = ncol(K)
	d = 2 * p
  	
	# First and second half
	fh = 1:p
	sh = (p + 1):d
	
	# Standardizing data
	K = scale(K, T, T)

	# Centering classes
	y = scale(y, T, F)
	
  # Hermitian  kernel matrix
  H = t(K) %*% K
	  
	# Vector of ones
	oneVec = as.matrix(rep(1, n))
	
	# System initial coefficients
	beta = as.matrix(zeros(n))
	beta0 = mean(y)

	# Initial active set
	A = NULL
	
	# Initial gradient
	g = K %*% y - beta0
	
	# Initial cvec (wrong value to enter while)
	cVec = beta0 + oneVec
	
	# At each iteration, lambda is reduced by the following amount: 
	lambdaReductionFactor = lambda # * 0.01
	
	# Initial debug
	if (debug)
	  printLog(lambda, beta0, beta)

	# Solution paths
	lambdas = c()
	betas0 = c()
	betas = c()
	norms = c()	
	errors = c()
	
	# We will solve the system using the Newton-Raphson method: 
	# x_{n + 1} = x_n - F(x_n) / J(x_n)
	while (lambda > lambdaThreshold)
  {
	  iter = 0
	  
	  while (!verifySolution(g, A, lambda, beta0, beta, cVec, oneVec, n) & iter < 1e3)
	  {
  		# Creating Jacobian J
  		J = zeros(d, d)
  		J[fh, fh] = diag(p)
  		J[fh, sh] = diag(t(K) %*% oneVec / n)
  		J[sh, fh] = J[fh, sh] * n
  		J[sh, sh] = H
  		
  		# Common vec
  		cVec = y - K %*% beta
  		
  		# Defining the function F
  		f = zeros(d, 1)
  		f[fh] = (-beta0 + t(cVec) %*% oneVec / n) / p
  		f[sh] = (K %*% (cVec - beta0 * oneVec)) - sign(beta) * lambda
  
      # Adding small constants to overcome singularity
  		J = J + diag(rep(1e-7, d))
  		
  		# Applying solver
  		delta = solve(J, f)
  
  		# x_{n + 1} = x_n + delta
  		beta0 = beta0 + eta * sum(delta[fh])
  		beta = beta + eta * delta[sh]
    	
  		# Updating cVec
  		cVec = y - K %*% beta
  		  		
    	# Computing active set
    	A = which(!isZero(beta))
    		
    	# Gradient
    	g = K %*% (y - K %*% beta - beta0 * oneVec)

    	# Counting iterations
    	iter = iter + 1
    	
    	printf("Iter %d) %e\n", iter, lnorm(cVec))
	  }
    
	  # Infos
	  bnorm = lnorm(beta, 1)
	  error = lnorm(y - K %*% beta - beta0)
        
	  # Return data
	  lambdas = rbind(lambdas, lambda)
    betas0 = rbind(betas0, beta0)
    betas = rbind(betas, beta)
    norms = rbind(norms, bnorm)
    errors = rbind(errors, error)	  

  	# Debugging solutions
  	if (debug)
  	  printLog(lambda, beta0, beta)
  
  	# Decreasing lambda		
  	lambda = lambda - lambdaReductionFactor
	}
	return(list(lambdas = lambdas, betas0 = betas0, betas = betas, norms = norms, 
	  erros = errors, K = K, y = y))
}

klassoTest <- function(mu = 5, n = 4)
{
  n = 100
  p = 10
  X = matrix(rnorm(n * p), ncol = p)
  y = X[,4] + 0.566 + X[,8] + 0.432

	klasso(X, y)
}
