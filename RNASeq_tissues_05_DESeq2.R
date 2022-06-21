library("DESeq2") 
library("vsn")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("BiocParallel")
register(SnowParam(8))
library("tximport")
library("readr")
library("rjson")
library(openxlsx)
library(doParallel) 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)


setwd("D:/Work/! projects/MTM/RNAseq/12 - DESeq2 analysis/MTM_muscle_both batches combined - group")

sessionInfo = paste("session-info ",Sys.time(),".txt", sep="")
sessionInfo = gsub(":", "-", sessionInfo)
sink(sessionInfo)
sessionInfo()
sink()

lnames = load(file = "txi.RData")
lnames
names(samples)

# subset for columns of interest
samples1 = samples[ , c("RNAseq_muscle_batch","Site","Sex","Group","HD","mean_CAG_small", "mean_CAG_large","MTM_DBS","MTM_DBS_zero" ,"MTM_weight", "MTM_height",  "MTM_BMI" , "MTM_age_months", "MTM_tfcscore",  "MTM_motscore", "MTM_fiscore", "MTM_indepscl", "rel_mtDNA")]
head(samples1)

cutBMI = cut(samples1$MTM_BMI, 5)
cutAGE = cut(samples1$MTM_age_months, 5)
factoring = cbind(as.character(cutBMI), as.character(cutAGE))
rownames(factoring) = rownames(samples1)
colnames(factoring) = c("MTM_BMI", "MTM_age_months")
write.csv(factoring, file="factoring_for_continuous_variables.csv")

cutBMI = cut(samples1$MTM_BMI, 5, c(paste("BMI",1:5,sep="")))
cutAGE = cut(samples1$MTM_age_months, 5, c(paste("AGE",1:5,sep="")))

samples2 = cbind(cutBMI, cutAGE, samples1)
head(samples2)

design = ~ RNAseq_muscle_batch + Sex + Site + cutBMI + cutAGE  + Group
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples2, design = design)
ddsTxi

dds <- ddsTxi 
dds

eval.dds <- DESeq(dds, parallel=TRUE)
eval.dds
# coefficients:
resultsNames(eval.dds)

# get normalization and size factors
# normalization factor for each gene
normFactors = normalizationFactors(eval.dds)
# size factor for library used for the gene norm factors
nm <- assays(eval.dds)[["avgTxLength"]]
sf <- estimateSizeFactorsForMatrix(counts(eval.dds) / nm)
sf1 = data.frame(sf)
sf2 = cbind(rownames(sf1),sf1)
colnames(sf2) = c("sample.ID","size.factors")
write.csv(sf2, file="size_factors_used_to_correct_libraries.csv", row.names = F)

normalized.counts <- as.data.frame(counts(eval.dds, normalized=TRUE ))
normalized.counts <- cbind(rownames(normalized.counts), normalized.counts)
colnames(normalized.counts)[1] = "geneid"
write.table(normalized.counts, file = "MTM_muscle_combined_outliers_removed_normalized_counts_DESeq2_from_proper_design_matrix.csv", sep=",", row.names=F)

VTSdat <- varianceStabilizingTransformation(eval.dds, blind=TRUE)
VTSdat1 <- assay(VTSdat)
VTSdat1 <- cbind(rownames(VTSdat1), VTSdat1)
colnames(VTSdat1)[1] = "geneid"
write.table(VTSdat1, file = "MTM_muscle_combined_outliers_removed_normalized_counts_DESeq2_from_proper_design_matrix_VST.csv", sep=",", row.names=F)

# ------------------------------------------
# specify contrast here, c('factorName','numeratorLevel','denominatorLevel')
# need to do for every pair-wise comparison
traitcolumn = "Group"
trait1 = "pre"
trait2 = "control"
contrast = vector()
contrast[1]=paste(traitcolumn)
contrast[2]=paste(trait1)
contrast[3]=paste(trait2)
contrast
# run with shrunken LFCs to be able to cross-compare
res <- results(eval.dds, contrast=contrast, parallel = TRUE)
res
resShrunk <- lfcShrink(eval.dds, contrast=contrast, parallel = TRUE)
resShrunk
resShrunk2 <- lfcShrink(eval.dds, contrast = contrast, type="ashr", parallel = TRUE)
resShrunk2
#resShrunk2 <- lfcShrink(eval.dds, coef = 2, type="apeglm", parallel = TRUE)
resDat <- cbind(rownames(res), res)
colnames(resDat)[1] = "geneid"
write.table(resDat, file = paste(trait1,"_vs_", trait2,"_DESeq2_results.csv", sep=""), sep=",", row.names=F)
resDat2 <- cbind(rownames(resShrunk), resShrunk)
colnames(resDat2)[1] = "geneid"
write.table(resDat2, file = paste(trait1,"_vs_", trait2,"_DESeq2_results_shrunken_LFCs_normal.csv", sep=""), sep=",", row.names=F)
resDat3 <- cbind(rownames(resShrunk2), resShrunk2)
colnames(resDat3)[1] = "geneid"
write.table(resDat3, file = paste(trait1,"_vs_", trait2,"_DESeq2_results_shrunken_LFCs_ashr.csv", sep=""), sep=",", row.names=F)
# LFCs
pdf(file = paste(trait1,"_vs_", trait2,"_DEseq_results_LFCs.pdf", sep=""), w=30, h=10);
maxRes = max(res$log2FoldChange, na.rm=T) + 0.1 
maxResShrunk = max(resShrunk$log2FoldChange, na.rm=T) + 0.1
maxResShrunk2 = max(resShrunk2$log2FoldChange, na.rm=T) + 0.1
minRes = min(res$log2FoldChange, na.rm=T) - 0.1 
minResShrunk = min(resShrunk$log2FoldChange, na.rm=T) - 0.1
minResShrunk2 = min(resShrunk2$log2FoldChange, na.rm=T) - 0.1
par(mfrow=c(1,3))
plotMA(res, main="DESeq2 – unshrunken LFCs", ylim=c(minRes, maxRes))
plotMA(resShrunk, main="DESeq2 – shrunken LFCs (normal)", ylim=c(minResShrunk, maxResShrunk))
plotMA(resShrunk2, main="DESeq2 – shrunken LFCs (ashr)", ylim=c(minResShrunk2, maxResShrunk2))
dev.off()
# ------------------------------------------



# ------------------------------------------
# specify contrast here, c('factorName','numeratorLevel','denominatorLevel')
# need to do for every pair-wise comparison
traitcolumn = "Group"
trait1 = "early"
trait2 = "control"
contrast = vector()
contrast[1]=paste(traitcolumn)
contrast[2]=paste(trait1)
contrast[3]=paste(trait2)
contrast
# run with shrunken LFCs to be able to cross-compare
res <- results(eval.dds, contrast=contrast, parallel = TRUE)
res
resShrunk <- lfcShrink(eval.dds, contrast=contrast, parallel = TRUE)
resShrunk
resShrunk2 <- lfcShrink(eval.dds, contrast = contrast, type="ashr", parallel = TRUE)
resShrunk2
#resShrunk2 <- lfcShrink(eval.dds, coef = 2, type="apeglm", parallel = TRUE)
resDat <- cbind(rownames(res), res)
colnames(resDat)[1] = "geneid"
write.table(resDat, file = paste(trait1,"_vs_", trait2,"_DESeq2_results.csv", sep=""), sep=",", row.names=F)
resDat2 <- cbind(rownames(resShrunk), resShrunk)
colnames(resDat2)[1] = "geneid"
write.table(resDat2, file = paste(trait1,"_vs_", trait2,"_DESeq2_results_shrunken_LFCs_normal.csv", sep=""), sep=",", row.names=F)
resDat3 <- cbind(rownames(resShrunk2), resShrunk2)
colnames(resDat3)[1] = "geneid"
write.table(resDat3, file = paste(trait1,"_vs_", trait2,"_DESeq2_results_shrunken_LFCs_ashr.csv", sep=""), sep=",", row.names=F)
# LFCs
pdf(file = paste(trait1,"_vs_", trait2,"_DEseq_results_LFCs.pdf", sep=""), w=30, h=10);
maxRes = max(res$log2FoldChange, na.rm=T) + 0.1 
maxResShrunk = max(resShrunk$log2FoldChange, na.rm=T) + 0.1
maxResShrunk2 = max(resShrunk2$log2FoldChange, na.rm=T) + 0.1
minRes = min(res$log2FoldChange, na.rm=T) - 0.1 
minResShrunk = min(resShrunk$log2FoldChange, na.rm=T) - 0.1
minResShrunk2 = min(resShrunk2$log2FoldChange, na.rm=T) - 0.1
par(mfrow=c(1,3))
plotMA(res, main="DESeq2 – unshrunken LFCs", ylim=c(minRes, maxRes))
plotMA(resShrunk, main="DESeq2 – shrunken LFCs (normal)", ylim=c(minResShrunk, maxResShrunk))
plotMA(resShrunk2, main="DESeq2 – shrunken LFCs (ashr)", ylim=c(minResShrunk2, maxResShrunk2))
dev.off()
# ------------------------------------------



# ------------------------------------------
# specify contrast here, c('factorName','numeratorLevel','denominatorLevel')
# need to do for every pair-wise comparison
traitcolumn = "Group"
trait1 = "early"
trait2 = "pre"
contrast = vector()
contrast[1]=paste(traitcolumn)
contrast[2]=paste(trait1)
contrast[3]=paste(trait2)
contrast
# run with shrunken LFCs to be able to cross-compare
res <- results(eval.dds, contrast=contrast, parallel = TRUE)
res
resShrunk <- lfcShrink(eval.dds, contrast=contrast, parallel = TRUE)
resShrunk
resShrunk2 <- lfcShrink(eval.dds, contrast = contrast, type="ashr", parallel = TRUE)
resShrunk2
#resShrunk2 <- lfcShrink(eval.dds, coef = 2, type="apeglm", parallel = TRUE)
resDat <- cbind(rownames(res), res)
colnames(resDat)[1] = "geneid"
write.table(resDat, file = paste(trait1,"_vs_", trait2,"_DESeq2_results.csv", sep=""), sep=",", row.names=F)
resDat2 <- cbind(rownames(resShrunk), resShrunk)
colnames(resDat2)[1] = "geneid"
write.table(resDat2, file = paste(trait1,"_vs_", trait2,"_DESeq2_results_shrunken_LFCs_normal.csv", sep=""), sep=",", row.names=F)
resDat3 <- cbind(rownames(resShrunk2), resShrunk2)
colnames(resDat3)[1] = "geneid"
write.table(resDat3, file = paste(trait1,"_vs_", trait2,"_DESeq2_results_shrunken_LFCs_ashr.csv", sep=""), sep=",", row.names=F)
# LFCs
pdf(file = paste(trait1,"_vs_", trait2,"_DEseq_results_LFCs.pdf", sep=""), w=30, h=10);
maxRes = max(res$log2FoldChange, na.rm=T) + 0.1 
maxResShrunk = max(resShrunk$log2FoldChange, na.rm=T) + 0.1
maxResShrunk2 = max(resShrunk2$log2FoldChange, na.rm=T) + 0.1
minRes = min(res$log2FoldChange, na.rm=T) - 0.1 
minResShrunk = min(resShrunk$log2FoldChange, na.rm=T) - 0.1
minResShrunk2 = min(resShrunk2$log2FoldChange, na.rm=T) - 0.1
par(mfrow=c(1,3))
plotMA(res, main="DESeq2 – unshrunken LFCs", ylim=c(minRes, maxRes))
plotMA(resShrunk, main="DESeq2 – shrunken LFCs (normal)", ylim=c(minResShrunk, maxResShrunk))
plotMA(resShrunk2, main="DESeq2 – shrunken LFCs (ashr)", ylim=c(minResShrunk2, maxResShrunk2))
dev.off()
# ------------------------------------------




