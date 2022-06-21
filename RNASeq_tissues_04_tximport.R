library("tximport")
library("readr")
library("rjson")
library(doParallel) 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

setwd("D:/Work/! projects/MTM/RNAseq/07 - tximport and batch effects/MTM_muscle_both batches combined")

dir <- "D:/Work/! projects/MTM/RNAseq/03 - quant_salmon/muscle_combined"

samples <- read.csv("MTM_all_metadata_2019-02-21_muscle_outliers_removed.csv", header=TRUE)
names(samples)
rownames(samples) <- samples$RNAseq_muscle_Lab.Sample
files <- file.path(dir, samples$RNAseq_muscle_Lab.Sample, "quant.sf")
names(files) <- samples$RNAseq_muscle_Lab.Sample

tx2gene <- read.csv("tx2gene_GRCh38.90.csv")

txi <- tximport(files, type="salmon", tx2gene=tx2gene, varReduce=TRUE)

save(samples, files, txi, file = "txi.RData")

sessionInfo = paste("session-info ",Sys.time(),".txt", sep="")
sessionInfo = gsub(":", "-", sessionInfo)
sink(sessionInfo)
sessionInfo()
sink()

