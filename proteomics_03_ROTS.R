# early vs control
library(doParallel) 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
library("readr")
library(openxlsx)
library(ROTS)
library(MKmisc)
library(ggrepel)

setwd("D:/Work/! projects/MTM/proteomics/05 - dysregulation/adipose protein")

sessionInfo = paste("session-info ",Sys.time(),".txt", sep="")
sessionInfo = gsub(":", "-", sessionInfo)
sink(sessionInfo)
sessionInfo()
sink()

# load data
set1Data = read.csv(file="40-087_Adipose_phase_4_protein_QS_edit_limma_sex_BMI_age_corrected.csv");
names(set1Data)

# prepare
rownames(set1Data) = set1Data$PID_edit
# remove descriptors, keep only data
data = data.matrix(set1Data[ , -c(1)])

# load trait data
traitData1 = read.csv(file="DecoderRing_MTM_proteomics_combined_adipose.csv");
names(traitData1)
rownames(traitData1) = traitData1$sampleID
# match
order = match(colnames(data), rownames(traitData1))
traitData = traitData1[order , ]
table(colnames(data)==rownames(traitData))

# condition to test for
# gives different results for subset data with only the two groups of interest. But as with RNAseq (DESeq2), probably better to have all data for fitting. Therfore, I didn't subset for more contrasts
names(traitData)
groups = traitData$Status_ctr_early_num
groups

# run ROTS 
# bootstrap (B) should be >1000
# top-ranked features for reproducibility optimization (K) should be > ?? If no value is given, 1/4 of the features are used.
# log A logical (deafult TRUE) indicating whether input data is log2 scaled. This information is only used to calculate log fold change.

results = ROTS(data = data, groups = groups , B = 1000 , seed = 1234567890, log = T)
names(results)
sink(file="summary_early_vs_control.txt")
summary(results, fdr=0.05)
print("reproducibility Z should be larger than 2, according to manual")
sink()


Data = results$data

# calculate fold change
order2 = match(colnames(Data), rownames(traitData1))
FCtrait = traitData1[order2, ]
table(colnames(Data)==rownames(FCtrait))
# class 1 is reference (e.g. control), class 2 is group of interest
groups
class1 = which(FCtrait$Status_ctr_early_num=="1",arr.ind = T)
class2 = which(FCtrait$Status_ctr_early_num=="3",arr.ind = T)
# for log transformed data
logFC = rowMeans(Data[ , class2]) - rowMeans(Data[ , class1])
# for normal data
# logFC = rowMeans(log2(Data[ , class2]+1)) - rowMeans(log2(Data[ , class1]+1))

# export results
resultsDat = as.data.frame(logFC)
resultsDat$optimised.ROTS.statistics.d.value = results$d
resultsDat$p.value = results$pvalue
resultsDat$FDR = results$FDR

annot = read.csv(file="40-087_Adipose_phase_4_protein_QS_edit_annotation.csv")
names(annot)
traitRows = match(rownames(resultsDat), annot$PID_edit);
datAnnot = annot[traitRows, ];

# resultsDat2 = cbind(datAnnot, resultsDat, Data)
# resultsDatWrite = resultsDat2[ order(resultsDat2$FDR) , ]
resultsDatWrite = cbind(datAnnot, resultsDat, Data)
write.xlsx(resultsDatWrite, file="40-087_Adipose_phase_4_protein_QS_edit_limma_sex_BMI_age_corrected_ROTS_results_early_vs_control.xlsx", row.names = F)

# this plots raw p-values instead of FDR, as otherwise a lot of points will fall onto the same spot. Somehow controlled by the p-value cutoff further down: -log10(0.001) = 3
resultsDatWrite$log.pval = -log10(resultsDatWrite$p.value)
x = ggplot(resultsDatWrite, aes(x = logFC, y =log.pval )) + 
    geom_point(size=0.5 )+
    theme_bw(base_size = 16) + # change theme
    xlab(expression("log2 expression")) + # x-axis label
    ylab(expression(" -log10(p-value)")) + # y-axis label
    # geom_vline(xintercept = c(-1,1), colour = "red") + # Add fold change cutoffs
    geom_hline(yintercept = 3, colour = "red") + # Add significance cutoffs
    geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
    scale_colour_gradient(low = "black", high = "black", guide = FALSE)+
   # pvalue and FC cutoff
   # geom_text_repel(data=subset(resultsDatWrite, abs(logFC)>1 & log.pval > 3),
   #                 aes( logFC, log.pval ,label=gene)) # add gene label
  # only pvalue cutoff 
   geom_text_repel(data=subset(resultsDatWrite, log.pval > 3),
                    aes( logFC, log.pval ,label=proteinName)) + # add gene label 
 coord_cartesian(ylim = c(0,(max(resultsDatWrite$log.pval) + 0.5)) )
pdf(file = "ROTS_result_volcano_early_vs_control.pdf", h=8, w=12);
print(x)
dev.off()


