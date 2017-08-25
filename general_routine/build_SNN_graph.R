library(ica)
library(tcltk)
library(FNN)
library(Matrix)
#From Seurat package
CalcSNNSparse <- function(
  cell.names,
  k.param,
  nn.large,
  nn.ranked,
  prune.SNN,
  print.output
) {
  n.cells <- length(cell.names)
  counter <- 1
  idx1 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  idx2 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  edge.weight <- vector(mode = "double", length = n.cells ^ 2 / k.param)
  id <- 1
  # fill out the adjacency matrix w with edge weights only between your target
  # cell and its k.scale*k.param-nearest neighbors
  # speed things up (don't have to calculate all pairwise distances)
  # define the edge weights with Jaccard distance
  if (print.output) {
    print("Constructing SNN")
    pb <- txtProgressBar(min = 0, max = n.cells, style = 3)
  }
  for (i in 1:n.cells) {
    for (j in 1:ncol(x = nn.large)) {
      s <- intersect(x = nn.ranked[i, ], y = nn.ranked[nn.large[i, j], ])
      u <- union(nn.ranked[i, ], nn.ranked[nn.large[i, j], ])
      e <- length(x = s) / length(x = u)
      if (e > prune.SNN) {
        idx1[id] <- i
        idx2[id] <- nn.large[i, j]
        edge.weight[id] <- e
        id <- id + 1
      }
    }
    if (print.output) {
      setTxtProgressBar(pb = pb, value = i)
    }
  }
  if (print.output) {
    close(con = pb)
  }
  idx1 <- idx1[! is.na(x = idx1) & idx1 != 0]
  idx2 <- idx2[! is.na(x = idx2) & idx2 != 0]
  edge.weight <- edge.weight[! is.na(x = edge.weight) & edge.weight != 0]
  w <- sparseMatrix(
    i = idx1,
    j = idx2,
    x = edge.weight,
    dims = c(n.cells, n.cells)
  )
  diag(x = w) <- 1
  rownames(x = w) <- cell.names
  colnames(x = w) <- cell.names
  return(w)
}

args = commandArgs(trailingOnly=TRUE)
file_in   = args[1]
file_out  = args[2]
num_ICs   = as.numeric(args[3])
k.param   = as.numeric(args[4])
k.scale   = as.numeric(args[5])
prune.SNN = as.numeric(args[6])

print("Reading data...")
cell.data = read.csv(file_in,header=FALSE)

print("Performing ICA...")
#icafast wants cells in rows and genes in columns
data.use = icafast(t(cell.data),nc=num_ICs,center=F)$S

n.cells <- nrow(x = data.use)
my.knn <- get.knn(
  data <- as.matrix(x = data.use),
  k = min(k.scale * k.param, n.cells - 1)
)
nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param-1)])
nn.large <- my.knn$nn.index

w = CalcSNNSparse(
  cell.names = colnames(cell.data),
  k.param = k.param,
  nn.large = nn.large,
  nn.ranked = nn.ranked,
  prune.SNN = prune.SNN,
  print.output = T
)

#Force w to be symmetric (equivalent to making an undirected graph)
w = ((w+t(w))/2)

#Save in matrix market format
writeMM(w,file_out)
