# ===============================
# About: Phase space related functions
# Dependences: Canon package
# Author: Rodrigo Mello, LCP
# Last revision 06/08/17
# ==============================

require("Canon")
sourceDir(canonPath())

loadPackages("hydroTSM", "tseriesChaos")

# ==============================
# Von Neumann entropy
# ==============================

maxEntropy <- function(series, maxm = 4, maxd = 10)
{
	ret = NULL

	for (m in 2:maxm)
  {
	  dvec = c()
	  
		for (d in 1:maxd)
		{
			data = embedd(series, m = m, d = d)
			plot(data)
			# H = futureStateEntropy(data)
			# H = vonNeumannEntropy(data, normalize = TRUE)
			# H = openBallsEntropy(data)
			# H = covarianceEntropy(data)
			H = vonNeumannEntropy(data)
			cat("m=",m, "d=",d, H, "\n")
    	dvec = c(dvec, H)
		}
	  ret = rbind(ret, dvec)
	}
	colnames(ret) = 1:ncol(ret)
	rownames(ret) = 1:nrow(ret) + 1
	return(ret)
}

Lambda <- function(Is_vector, data, lambda) {
	#Is_vector = heaviside(Is_vector)
	#cat("Is_vector = ", Is_vector, " lambda = ", lambda, "\n")
	term1 = von.Neumann.Entropy(data%*%diag(Is_vector)) 
	term2 = lambda * norm1(diag(Is_vector))
	#cat("term1 = ", term1, " term2 = ", term2, "\n")
	#-term1/term2
	return (term1 - term2)
}

shrinkage <- function(dHvn_dI, epsilon) {

	r = svd(dHvn_dI)
	dHvn_dI = r$d

	id = which(dHvn_dI >= -epsilon &&  dHvn_dI <= epsilon)
	if (length(id) > 0) dHvn_dI[id] = rep(0, length(id))

	id = which(dHvn_dI > epsilon)
	if (length(id) > 0) dHvn_dI[id] = dHvn_dI[id] - epsilon

	id = which(dHvn_dI < -epsilon)
	if (length(id) > 0) dHvn_dI[id] = dHvn_dI[id] + epsilon

	return (r$u %*% diag(dHvn_dI) %*% t(r$v))
}

my.optim <- function(series, m=100, lambda=0.1, alpha=0.1, iter=10) {
	M = embedd(series, m=m, d=1)
	hat_I = rnorm(mean=0, sd=0.01, n=m)

	for (i in 1:iter) {
		#cat("hat_I = ")
		#print(hat_I)

		cat("Objective = ", Lambda(hat_I, M, lambda), "\n")


		# Maximizing Von Neumann Entropy
		A = M %*% diag(hat_I)
		R = svd(A)
		Sigma = R$d/sum(R$d)

		#cat("Sigma\n")
		#print(Sigma)

		# Derivative Von Neumann Entropy in the direction of Lambda
		d_von.Neumann.Entropy__d_Lambda = rep(0, length(Sigma))
		d_von.Neumann.Entropy__d_Lambda[Sigma > 0] = -(1+log(Sigma[Sigma > 0]))
#		d_von.Neumann.Entropy__d_Lambda = d_von.Neumann.Entropy__d_Lambda[Sigma > 0]
#		r = length(d_von.Neumann.Entropy__d_Lambda)
#		if (r > 0) {
#			M_r = as.matrix(R$u[,1:r]) %*% diag(as.vector(R$d[1:r])) %*% t(as.matrix(R$v[,1:r]))
#		}
#		else {
#			M_r = as.matrix(1)
#		}
		
		#print(diag(d_von.Neumann.Entropy__d_Lambda))

		# Derivative Lambda in the direction of hat_I
		# 	Applying a numerical perturbation on hat_I
		dHvn_dI = diag(d_von.Neumann.Entropy__d_Lambda) %*% (t(M) %*% R$u %*% t(R$v))
		dHvn_dI = shrinkage(-dHvn_dI, lambda)
		#lambda = lambda * 0.9

		hat_I = hat_I + alpha * colSums(dHvn_dI)
		hat_I[hat_I > 1] = 1
		hat_I[hat_I < -1] = -1
	}

	hat_I
}
