library(ica)
args = commandArgs(trailingOnly=TRUE)
file_in  = args[1]
file_out = args[2]
num_ICs   = as.numeric(args[3])
cell.data = read.csv(file_in,header=FALSE)
data.use = icafast(t(cell.data),nc=num_ICs,center=F)$S
write.csv(data.use,file_out)
