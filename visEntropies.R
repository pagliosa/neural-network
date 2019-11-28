# ===============================
# About: Visualize entropies from C++ code
# Dependences: Canon pkg
# Author: Lucas Pagliosa
# Last revision: 10/09/17
# ==============================

require("Canon")
sourceDir(canonPath())

readEntropies <- function(path = "../VS Codes/Entropies/", regex = ".txt")
{
  path = sprintf("%s/%s", pwd(), path)
  fileData = getFiles(path, regex, T)
  dirs = list()

  for (dir in fileData)
  {
    dirs = addList(dirs, list(dirName = dir$path, nof = dir$nof))
    nof = dir$nof

    for (i in 1:nof)
      parseEntropies(dir$paths[i], dir$files[i], partial = F)
  }
}

parseEntropies <- function(path, filename = "Plot", partial = T)
{
  if (partial)
    path = sprintf("../VS Codes/Entropies/%s", path)
  
  printf("Reading file: %s...", path)

  file = as.matrix(read.matrix(path, header = F))
  parameters = file[1,]
  file = file[-1,]
  entropies = file[2:nrow(file),]
  
  maxm = parameters[1]
  maxd = parameters[2]
  
  plotAtt <- function(cid, ylab, pch)
  {
    plot(1, main = filename, type = "n", xlim = c(1, maxd), 
      ylim = range(entropies[,cid]), xlab = "Delay", ylab = ylab)

    for (j in 1:(maxm - 1))
    {
      s = (j - 1) * maxd + 1
      e = j * maxd
      lines(file[s:e, cid], type = "b", col = j, lty = 1, pch = pch)
    }
  }
  
  myPlot(plotFunc = function()
  {
    par(mfrow = c(1, 1))
    plotAtt(1, "VNE", 17)  
    plotAtt(2, "EBE", 20)  
    par(mfrow = c(1, 1))
  }, savePlot = F, filename = filename)
  printf("Done\n")
}
