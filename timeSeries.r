# ===============================
# About: Personal time series functions
# Dependences: utils.r
# Author: Lucas Pagliosa
# Last revision: 15/10/15
# ==============================

source("utils.r")
loadGraphics()
loadPackages("tseriesChaos", "nonlinearTseries", "deSolve")

createMackeyGlass <- function(nop = 1000, noise = 0, m = 5, tau = 17,
  onlyts = T, plot = F, save = F)
{
  filename = sprintf("%s/datasets/%s", canonPath(), "mackey-glass.txt")
  ts = read.table(filename)[,2]
  return(createTimeSeries(list(x = ts[1:nop]), m, tau, "Mackey-Glass", onlyts, plot,
    subSampleTS = nop, save = save))
}

createRandomSeries <- function(nop = 1000, min = -1, max = 1, name = "Guassian")
{
  return(list(ts = runif(nop, min, max), name = name))
}

createRandomWalk <- function(nop = 1000, m = 2, noise = 0, name = "RandomWalk",
  plot = T, save = F)
{
  walk = apply(matrix(runif(nop * m, min = -1, max = 1), ncol = m), 2, cumsum)

  if (plot)
    myPlot(function() {plot(walk, type = "s", pch = 15)}, save, "Randow Walk")
  return(list(ts = walk, name = name))
}

createLogistic <- function(nop = 1000, noise = 0, m = 2, tau = 1, r = 3.8,
  x0 = 0.5, onlyts = T, plot = F, save = F, N = -1)
{
  if (N != -1)
    nop = N + (m - 1) * tau
  if (nop < 10)
    error("Enter at least a ten length time series...")

  states = list(x = logisticMap(r = r, start = x0, n.sample = nop,
    n.transient = nop / 10, do.plot = F) + rnorm(nop, 0, noise))

  return(createTimeSeries(states, m, tau, "Logistic", onlyts, plot,
    subSampleTS = nop, save = save))
}

createRossler <- function(nop = 1000, m = 3, tau = 8, alpha = 0.2, beta = 0.2,
  gamma = 5.7, onlyts = T, plot = F, save = F, N = -1)
{
  if (N != -1)
    nop = N + (m - 1) * tau
  timeInc = 0.1
  states = rossler(a = alpha, b = beta, w = gamma, start = c(-2, -10, 0.2),
    time = seq(0, timeInc * nop, by = timeInc), do.plot = F)
  return(createTimeSeries(states, m, tau, "Rossler", onlyts, plot, save = save))
}

plotLorenzAttraction <- function(n = 4, color = F, save = T)
{
  range = c(-100, 100)
  plot(createLorenz(start = runif(3, range[1], range[2]), onlyts = F)$emb,
    xlim = 1.3 * range, ylim = 1.3 * range, type = "l")

  if (color)
    for (i in 2:n)
      lines(createLorenz(start = runif(3, range[1], range[2]), onlyts = F)$emb,
        col = i)
  else
    for (i in 2:n)
      lines(createLorenz(start = runif(3, range[1], range[2]), onlyts = F)$emb)
        #, lty = i)
  rc = recordPlot()
  myPlot(function(){replayPlot(rc)}, save, filename = "Lorenz-attraction-test")
}

createLorenzHand <- function(nop = 1000, m = 3, tau = 8, plot = T)
{
  ts = createLorenz(nop)$ts
  max = length(ts) - (m - 1) * tau
  emb = zeros(max, 3)

  for (i in 1:max)
    for (j in 1:m)
      emb[i, j] = ts[i + (j - 1) * tau]

  if (plot)
    plot3d(emb, type = "l")

  return(emb)
}

createLorenz <- function(nop = 1000, m = 3, tau = 8, sigma = 10, beta = 8/3,
  rho = 28, start = c(-13, -14, 47), timeInc = 0.01, onlyts = T, plot = F,
  save = F, pure = F, N = -1, oneframe = F)
{
  if (N != -1)
    nop = N + (m - 1) * tau
  states = lorenz(sigma = sigma, beta = beta, rho = rho, start = start,
    time = seq(0, timeInc * (nop - 1), by = timeInc), do.plot = F)
  if (pure)
    return(states)
  return(createTimeSeries(states, m, tau, "Lorenz", onlyts, plot, save = save,
    oneframe = oneframe))
}

createLorenzUsingDerivatives <- function(nop = 1000, tau = 8)
{
  ts = createLorenz(nop)

  min = tau + 1
  max = nop - tau

  emb = NULL

  for (i in min:max)
  {
    p = ts[i]
    nextp = ts[i + tau]
    prevp = ts[i - tau]

    fd = (nextp + p) / tau
    sd = (nextp + 2*p - prevp) / tau^2

    emb = rbind(emb, c(p, fd, sd))
  }
  return(emb)
}

compareLorenzMethods <- function(nop = 5e3, save = T)
{
  params = list(xlab = "x", ylab = "y", zlab = "z", type = "l")
  myPlot(function() {
    states = createLorenz(nop, pure = T)
    do.call(scatterplot3d, c(list(states$x, states$y, states$z), params))
  }, save, "Pure")

  myPlot(function() {
    do.call(plot, c(list(createLorenz(nop, onlyts = F)$emb), params))
  }, save, "Delay")

  myPlot (function() {
    do.call(plot, c(list(createLorenzUsingDerivatives(nop)), params))
  }, save, "Derivatives")
}

createHenon <- function(nop = 1000, m = 2, tau = 1, a = 1.4, b = 0.3,
  onlyts = T, plot = F, save = F)
{
  states = henon(n.sample = nop, n.transient = nop / 10, a = a, b = b,
    do.plot = F, start = c(-0.006423277,-0.473545134))
  return(createTimeSeries(states, m, tau, "Henon", onlyts, plot,
    subSampleTS = 100, save = save))
}

createIkeda <-function(nop = 1000, m = 2, tau = 1, onlyts = T, plot = F,
  save = F, N = -1)
{
  if (N != -1)
    nop = N + (m - 1) * tau
  states = ikedaMap(a = 0.85, b = 0.9, cc = 7.7, k = 0.4, start = runif(2),
    n.sample = nop, n.transient = nop / 10, do.plot = F)
  return(createTimeSeries(states, m, tau, "Ikeda", onlyts, plot,
    subSampleTS = 100, save = save))
}

createClifford <- function(nop = 1000, m = 2, tau = 1, onlyts = T, plot = F,
  save = F, N = -1)
{
  if (N != -1)
    nop = N + (m - 1) * tau
  a=-1.24458046630025
  b=-1.25191834103316
  c=-1.81590817030519
  d=-1.90866735205054
  states = cliffordMap(a = -1.4, b = 1.6, cc = 1, d = 0.7, start = c(0,0),
    n.sample = nop, n.transient = nop / 10, do.plot = F)
  return(createTimeSeries(states, m, tau, "Clifford", onlyts, plot, save = save,
    embTypePlot = "p"))
}

createSine <- function(nop = 1000, noise = 0, m = 2, tau = 1, A = 1,
  start = 0, end = 9, onlyts = T, plot = F, save = F, name = "Sine",
  embt = "l", N = -1)
{
  if (N != -1)
    nop = N + (m - 1) * tau
  states = list(x = A * sin(2*pi*seq(start, end, len = nop)) +
      rnorm(nop, 0, noise))
  return(createTimeSeries(states, m, tau, name, onlyts, plot,
    embTypePlot = embt, save = save))
}

createCombinedSine <- function(nop = 1e3, n = 10, A = NULL)
{
  if (is.null(A))
    A = seq(1e-2, 1, len = n)

  ts = NULL
  sum = 0

  for (i in n:1)
  {
    #sum = sum + createSine(nop, A = A[i])
    # ts = c(ts, sum)
    ts = c(ts, createSine(nop, A = A[i]))
  }

  emb = embedd(ts, 2, 1)
  plot(ts, type = "l")

  return(ts)
}

createSunspot <- function(nop = 1000, m = 2, tau = 1, nos = 100, onlyts = T,
  plot = F, save = F, N = -1)
{
  if (N != -1)
    nop = N + (m - 1) * tau
  ts = as.matrix(sunspot.month)[1:nop]
  states = list(x = spline(1:length(ts), ts, n = nos)$y)
  return(createTimeSeries(states, m, tau, "Sunspot", onlyts, plot,
    length(states$x), "l",  save))
}

createPendulum <- function(g = 9.81, L = 1, min = 0, max = 10, t_s = 0.01,
  plot = T, type = "p")
{
  angleSeq = seq(pi/2, pi/2, pi/5)
  angularVelocitySeq = seq(100, 100, 10)

  parameters = c(w = g/L)
  times = seq(min, max, by = t_s)
  pendulum = function(time, state, parameters)
  {
    with(as.list(c(state, parameters)),
    {
      dO = angularVelocity
      dV = w * -sin(angle)
      list(c(dO, dV))
    })
  }

  count = 1

  for (angle in angleSeq)
  {
    for (angularVelocity in angularVelocitySeq)
    {
      state = c(angle = angle, angularVelocity = angularVelocity)
      ret = ode(y = state, times = times, func = pendulum, parms =
          parameters)[,2:3]

      if (plot)
      {
        if (count == 1)
        {
          plot(ret[,2], ret[,1], type = type, pch = 19)
          #, xlim = range(angleSeq), ylim = range(angularVelocitySeq))
        }
        else
        {
          points(ret[,2], ret[,1], type = type, pch = 19)
        }
      }
      count = count + 1
    }
  }
  return(ret)
}

createTimeSeries <- function(states, m, tau, name, onlyts, plot,
  subSampleTS = length(states$x), embTypePlot = "p", save = F, oneframe = F)
{
  ts = as.vector(states$x)
  retList = list(ts = ts, name = name, m = m, tau = tau)

  cat("Creating ", sprintf("%s... ", name))
  if (onlyts)
  {
    if (plot)
      myPlot(function(){plot(ts[1:subSampleTS], main = name, type = "l",
        ylab = "Observation", xlab = "Time")}, save, sprintf("./%s.pdf", name))
    cat("Done\n")
    return(retList)
  }

  emb = embedd(ts, m = m, d = tau)
  cnames = c("x(t)")

  for (i in 1:(m - 1))
    cnames = c(cnames, sprintf("x(t + %d)", i * tau))
  colnames(emb) = cnames

  if (plot)
    visualizePhaseSpace(states, emb, name, subSampleTS, embTypePlot, save,
      oneframe)
  cat("Done\n")
  printf("Number of embedding states: %d\n", nrow(emb))

  retList$emb = emb
  return(retList)
}

visualizePhaseSpace <- function(ts, emb, name, subSampleTS, embTypePlot = "p",
  save = T, oneframe = F)
{
  pf = plot
  if (ncol(emb) == 3)
    pf = plot3d
  tf = function(){plot(ts$x[1:subSampleTS], type = "l",
      ylab = "Observation x(t)", xlab = "Time t")}
  ef = function(){pf(emb, type = embTypePlot)}

  if (oneframe)
  {
    myPlot(function()
    {
      par(mfrow=c(1, 2))
      tf()
      ef()
      par(mfrow=c(1, 1))
    }, save, sprintf("%s-ts-embed", name), "./")
  }
  else
  {
    myPlot(tf, save, sprintf("%s-timeSeries", name), "./")
    myPlot(ef, save, sprintf("%s-embed", name), "./")
  }
}

saveAllPaperIIDPlots <- function()
{
  ret = createLogistic(plot = T, onlyts = F, m = 2, tau = 1, save = T)
  ret = createIkeda(plot = T, onlyts = F, m = 2, tau = 1, save = T)
  ret = createHenon(plot = T, onlyts = F, m = 2, tau = 1, save = T)
  ret = createLorenz(plot = T, onlyts = F, m = 3, tau = 8, save = T)
  ret = createRossler(plot = T, onlyts = F, m = 3, tau = 8, save = T)
  ret = createSunspot(plot = T, onlyts = F, m = 2, tau = 1, save = T)
  ret = createSine(plot = T, onlyts = F, m = 2, tau = 1, save = T)
  ret = createSine(plot = T, onlyts = F, noise = 0.05,
    name = "Sine-default-0.05", save = T)
  ret = createSine(plot = T, onlyts = F, noise = 0.05, m = 4, tau = 5,
    name = "Sine-learned-0.05", save = T)
  ret = createSine(plot = T, onlyts = F, noise = 0.1,
    name = "Sine-default-0.1", save = T)
  ret = createSine(plot = T, onlyts = F, noise = 0.1, m = 4, tau = 10,
    name = "Sine-learned-0.1", save = T)
}

# ===============================
# Create array of time series. Each time series has its information, such as:
# ts: time series array
# mMax: max embeed dimension
# mInc: increment used in the search for ideal dimension
# mMax: max time delay
# mInc: increment used in the search for time delay
# tsName: time series name for logging
# ==============================
createTimeSeriesArrayIidPaper <- function()
{
  ts = NULL

  if (F)
  {
    ts = addList(ts, list(ts = createLogistic(), mMax = 6, dMax = 6,
      tsName = "Logistic", relevantNeighbors = 1.96))

    ts = addList(ts, list(ts = createHenon(), mMax = 6, dMax = 6,
      tsName = "Henon", relevantNeighbors = 1.3))

    ts = addList(ts, list(ts = createLorenz(), mMax = 6, dMax = 10,
      tsName = "Lorenz", relevantNeighbors = 1.5))

    ts = addList(ts, list(ts = createRossler(), mMax = 6, dMax = 15,
      tsName = "Rossler", relevantNeighbors = 4))

    ts = addList(ts, list(ts = createIkeda(), mMax = 6, dMax = 6,
      tsName = "Ikeda", relevantNeighbors = 1.5))

    ts = addList(ts, list(ts = createSine(), mMax = 10, dMax = 10,
      tsName = "Sine", relevantNeighbors = 0.5))

    ts = addList(ts, list(ts = createSine(noise = 0.05), mMax = 10, dMax = 10,
      tsName = "Sine-0.05", relevantNeighbors = 0.8))

    ts = addList(ts, list(ts = createSine(noise = 0.1), mMax = 10, dMax = 10,
      tsName = "Sine-0.1", relevantNeighbors = 1.2))
  }

  ts = addList(ts, list(ts = createSunspot(nos= 1000), mMax = 10, dMax = 15,
    tsName = "Sunspot", relevantNeighbors = 0.8))

  return(ts)
}

createTimeSeriesArray <- function()
{
  ts = NULL

  ts = addList(ts, list(data = createLogistic(), name = "Logistic"))

  ts = addList(ts, list(data = createHenon(100), name = "Henon"))

  ts = addList(ts, list(data = createLorenz(), name = "Lorenz"))

  ts = addList(ts, list(data = createRossler(), name = "Rossler"))

  ts = addList(ts, list(data = createIkeda(), name = "Ikeda"))

  ts = addList(ts, list(data = createSine(), name = "Sine"))

  ts = addList(ts, list(data = createSine(noise = 0.05), name = "Sine-0.05"))

  ts = addList(ts, list(data = createSine(noise = 0.1), name = "Sine-0.1"))

  ts = addList(ts, list(data = createSunspot(nos= 1000), name = "Sunspot"))

  return(ts)
}

saveTS <- function(ts, filename = "timeSeries.txt", append = T)
{
  write(ts, filename, append = append, ncolumns = 1)
}

saveAllTimeSeriesFiles <- function(dir = './TimeSeries/')
{
  createDirPath(dir)
  timeSeries = createTimeSeriesArray()

  for (i in 1:length(timeSeries))
  {
    ts = timeSeries[[i]]
    saveTS(ts$data$ts, sprintf("%s%s.txt", dir, ts$name))
  }
}

AMI <- function(ts = createLorenz(), plot = T, save = F)
{
  if (is.list(ts))
    ts= ts$ts

  ami = as.matrix(mutual(ts, lag.max = 30))

  id = 0
  prev = ami[1]

  for (i in 2:length(ami))
  {
    curr = ami[i]
    if (prev < curr)
    {
      id = i - 2
      break
    }
    prev = curr
  }

  if (plot)
    myPlot(function(){plot(ami, type = "h")}, save, "AMI")
  return(list(ami = ami, id = id))
}

FNN <- function(ts, tau = 1, maxm = 10, plot = T, save = F)
{
  if (is.list(ts))
    ts = ts$ts

  fnn = false.nearest(ts, maxm, tau, floor(length(ts) / 18))
  if (plot)
    myPlot(function(){plot(fnn, type = "h")}, save, "FNN")
  return(fnn)
}

extimateEmbeddinParameters <- function(ts, plot = T)
{
  ami = AMI(ts, plot)
  fnn = FNN(ts, plot = plot)
  return(list("ami" = ami, "fnn" = fnn))
}

plotRandomEmbeddingPoints <- function(embed, nop, plot = T, save = F)
{
  points = sample(1:nrow(embed), nop)
  randomPoints = embed[points,]

  mplot <- function()
  {
    plot(embed, cex = 0.5)
    points(randomPoints, pch = 1, cex = 4)
  }

  if (plot)
    myPlot(mplot())
  if (save)
    savePDF(mplot(), filename = "Random Embedding Plot")
}

showBifurcation <- function()
{
  LotVmod <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      dx = x*(alpha - beta*y)
      dy = -y*(gamma - delta*x)
      return(list(c(dx, dy)))
    })
  }

  n <- 100 # number of simulations
  param.name <- "gamma" # choose parameter to perturb
  param.seq <- seq(0,1,length = 50) # choose range of parameters

  Pars <- c(alpha = 1, beta = .001, gamma = 1, delta = .001)
  Time <- seq(0, 10, length = n)
  State <- c(x = .5, y = .9)

  param.index <- which(param.name == names(Pars))
  out <- list()
  for (i in 1:length(param.seq))
    out[[i]] <- matrix(0, n, length(State))

  for (i in 1:length(param.seq)) {
    # set params
    Pars.loop <- Pars
    Pars.loop[param.index] <- param.seq[i]
    # converge
    init <- ode(State, Time, LotVmod, Pars.loop)
    # get converged points
    out[[i]] <- ode(init[n,-1], Time, LotVmod, Pars.loop)[,-1]
  }

  range.lim <- lapply(out, function(x) apply(x, 2, range))
  range.lim <- apply(do.call("rbind", range.lim), 2, range)
  plot.variable <- "x" # choose which variable to show
  plot(0, 0, pch = "", xlab = param.name, ylab = plot.variable,
       xlim = range(param.seq), ylim = range.lim[,plot.variable], type = "l")
  for (i in 1:length(param.seq))
    points(rep(param.seq[i], n), out[[i]][,plot.variable], type = "l")
}

predictTimeSeries <- function(tsdata, plot = T, save = T, name = tsdata$name,
  window = 10, recursive = TRUE)
{
	# Testing recursivelly
	return(dwnnManager(tsdata$emb, tsdata$tau, window = window, plot = plot,
	  save = save, name = name, recursive = recursive))
}
