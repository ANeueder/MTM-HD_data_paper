library("tximport")
library("readr")
library("rjson")
library(doParallel) 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

setwd("D:/Work/! projects/MTM/fibroblasts/RNAseq/05 - tximport and batch effects")

dir <- "D:/Work/! projects/MTM/fibroblasts/RNAseq/03 - salmon/salmon"

samples <- read.csv("MTM_all_metadata_2019-09-25_fibro_edit.csv", header=TRUE)
names(samples)
rownames(samples) <- samples$fibro_RNAseqID
files <- file.path(dir, samples$fibro_RNAseqID, "quant.sf")
names(files) <- samples$fibro_RNAseqID

tx2gene <- read.csv("tx2gene_GRCh38.90.csv")

txi <- tximport(files, type="salmon", tx2gene=tx2gene, varReduce=TRUE)

save(samples, files, txi, file = "txi.RData")

sessionInfo = paste("session-info ",Sys.time(),".txt", sep="")
sessionInfo = gsub(":", "-", sessionInfo)
sink(sessionInfo)
sessionInfo()
sink()


