# ===============================
# About: Linear Model Selection and Regularization
# Dependences: Canon pkg, svm.r
# Author: Lucas Pagliosa
# Last revision: 10/01/18
# ==============================

require("Canon")

sourceFiles("svm.r")
loadPackages("lars", "genlasso", "glmnet", "entropy", "infotheo", "c3net")

OLS <- function(data)
{
  if (missing(data))
  {
    seq = 1:9
    cs = cumsum(seq)
    data = cbind(seq, cs)
  }
  n = nrow(data)
  x = data[,1]
  X = matrix(c(x, x^2), ncol = 2)
  y = data[,2]
  Xt = t(X)
  inv = solve(Xt %*% X)
  B = inv %*% Xt %*% y
  plot(x, y)
  lines(x, X %*% B)
  return(B)
}

testRandomLasso <- function(n = 100, p = 5)
{
  X = matrix(rnorm(n * p), ncol = p)
  y = X[,1] + rnorm(n) + X[,2] + rnorm(n)
  object <- lars(X, y, type="lasso")
  plot(object)
  coef4.1 <- coef(object, s = 4.1, mode = "norm")
  return(list("obj" = object, "beta" = coef4.1))
}

testLassoIris <- function()
{
  X = as.matrix(iris[,1:4])
  y = unclass(iris[,5])
  y = y[1:length(y)]
  object <- lars(X, y, type = "lasso")
  plot(object)
  coef4.1 <- coef(object, s = 4.1, mode = "norm")
  return(list("obj" = object, "beta" = coef4.1))
}

testGenLassoRandom <- function(n = 100, p = 10)
{
  X = matrix(rnorm(n * p), ncol = p)
  y = X[,1] + 0.566 + X[,2] + 0.432
  D = diag(1, p)
  out = genlasso(y = y, X = X, D = D)
  plot(out)
  beta = coef(out, lambda = sqrt(n * log(p)))$beta

  return(list("X" = X, "y" = y, "beta" = beta))
}

testGenLassoIris <- function()
{
  X = as.matrix(iris[,1:4])
  y = drop(unclass(iris[,5]))
  n = nrow(X)
  p = ncol(X)
  D = diag(1, p)
  out = genlasso(y = y, X = X, D = D)
  plot(out)
  beta = coef(out, lambda = sqrt(n * log(p)))$beta

  return(list("X" = X, "y" = y, "beta" = beta))
}

errorTest <- function(data, beta = NULL)
{
  fesb = data$fes$beta
  printf("Lasso: %f\n", mean((data$ytest - data$xtest %*% 
    fesb[2:length(fesb)])^2))
  if (!is.null(beta))
    printf("Guess: %f\n", mean((data$ytest - data$xtest %*% beta)^2))
}

glmnetLasso <- function(x, y, percentage = 0.8, alpha = 1,
  df = -1, plot = T)
{
  if (missing(x) | missing(y))
  {
    n = 100
    p = 10
    x = matrix(rnorm(n * p), ncol = p)
    y = x[,4] + 0.566 + x[,8] + 0.432 + x[,3]
  }

  n = nrow(x) # Number of observations
  p = ncol(x) # Number of predictors

  # Defining train and test features
  trainSize = floor(n * percentage)
  trainSample = sample(1:n, trainSize)
  xtrain = x[trainSample,]
  xtest = x[-trainSample,]

  # Defining train and test responses
  ytrain = y[trainSample]
  ytest = y[-trainSample]

  # Lasso
  fit = glmnet(xtrain, ytrain, alpha = alpha, intercept = TRUE)

  # Cross Validation Lasso
  cvfit = cv.glmnet(xtrain, ytrain, alpha = alpha, grouped = TRUE)

  # Plot
  if (plot)
  {
    par(mfrow = c(2, 1))
    plot(fit, xvar = "norm", label = TRUE)
    plot(fit, xvar = "lambda", label = TRUE)
    par(mfrow = c(1, 1))
    plot(cvfit)
  }

  # Lambda
  minLambda = cvfit$lambda.min
  fesLambda = cvfit$lambda.1se

  # Beta
  minBeta = as.matrix(coef(fit, s = minLambda))
  fesBeta = as.matrix(coef(fit, s = fesLambda))

  # Prediction for min lambda and one-standard-error lambda
  # lambda.1se gives the most regularized model such that
  # error is within one standard error of the minimum
  minPred = predict(cvfit, newx = xtest, s = "lambda.min")
  fesPred = predict(cvfit, newx = xtest, s = "lambda.1se")

  # Residual sum of squares (RSS)
  minRSS = sum((ytest - minPred)^2)
  fesRSS = sum((ytest - fesPred)^2)

  min = list("lambda" = minLambda,
    "beta" = minBeta,
    "pred" = minPred,
    "rss" = minRSS)
  fes = list("lambda" = fesLambda,
    "beta" = fesBeta,
    "pred" = fesPred,
    "rss" = fesRSS)

  if (df >= 0 & df <= n)
  {
    id = which(fit$df >= df)[1]
    dfLambda = fit$lambda[id]
    dfPred = predict(cvfit, newx = xtest, s = dfLambda)
    df = list("lambda" = dfLambda,
    "beta" = as.matrix(coef(fit, s = dfLambda)),
    "pred" = dfPred,
    "rss" = sum((ytest - dfPred)^2))
    if (plot)
      abline(v = log(dfLambda), lty = 3, col = 4)
  }
  else
    df = NULL

  ret = list("xtrain" = xtrain, "xtest" = xtest, "ytrain" = ytrain,
    "ytest" = ytest, "min" = min, "fes" = fes, "df" = df)
  return(ret)
}

irisLasso <- function()
{
  x = as.matrix(iris[,1:4])
  y = drop(unclass(iris[,5]))
  ret = glmnetLasso(x, y)
}

tsLasso <- function(ts = "Log", m = 5, tau = 1, alpha = 1)
{
  family = c("Logistic", "Ikeda", "Henon", "Lorenz", "Rossler", "Sunspot", 
    "Sine")
  ts = family[pmatch(ts, family)]
  nod = m + 1
  stringFunc = sprintf("create%s(m = %d, tau = %d, onlyts = F)$emb", ts, nod,
    tau)
  ps = eval(parse(text = stringFunc))
  x = ps[,1:m]
  y = ps[,nod]
  ret = glmnetLasso(x, y, alpha = alpha)
  betas = cbind(ret$min$beta, ret$fes$beta)
  colnames(betas) = c("min", "fes")
  return (betas)
}

kernelLasso <- function(n = 10, p = 20)
{
  x = matrix(rnorm(n * p), ncol = p)
  y = x[,4] + 0.566 + x[,8] + 0.432
    
  # Computing mutual information
  data = cbind(x, y)
  mim = makemim(t(data))
  # x = mim[1:p, 1:p]
  # x = x + diag(p)
  # y = as.matrix(mim[p+1, 1:p])
  model = glmnetLasso(x, y)

  return (model)
}
