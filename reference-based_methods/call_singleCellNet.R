source("run_scc.R")
args <- commandArgs(trailingOnly = TRUE)

ad = read_h5(sprintf("../write/%s.h5",args[1]))
obsRef<-ad$obs
expRef<-ad$rawX
obsRef<-obsRef[,"cell"]# tuple

ad = read_h5(sprintf("../write/%s.h5",args[2]))
obsPred<-ad$obs
expPred<-ad$rawX
obsPred<-rep('UNK',dim(expPred)[2])

x<-run_singleCellNet(obsRef,expRef,expPred)
write.csv(rownames(x)[argmax(x, rows = FALSE)],sprintf("%s_%s_singleCellNet.csv", args[1], args[2]))