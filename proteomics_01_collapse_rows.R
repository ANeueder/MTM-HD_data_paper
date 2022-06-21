
workingDir = "D:/Work/! projects/MTM/proteomics/02 - collapse rows and batch PCA/adipose peptide";
setwd(workingDir);
library(WGCNA);
options(stringsAsFactors = FALSE);
allowWGCNAThreads()

set1DATA = read.csv(file = "40-087_Adipose_phase_4_peptide_QS_edit.csv", header=TRUE);
n = length(set1DATA[,1])
RowIdentifier=paste("rowIdent", 1:n, sep="")
datset1DATA = data.frame(RowIdentifier, set1DATA)
names(datset1DATA)

datET= datset1DATA[,-c(1:13)]

# identifier, e.g. probeset ID
rowID= datset1DATA[,1]
dimnames(datET)[[1]]=rowID

# the column to collapse to
rowGroup=datset1DATA[,c("Peptide_upper")]

collapse.object = collapseRows(datET, rowGroup, rowID, method="absMaxMean", connectivityBasedCollapsing=FALSE,  selectFewestMissing=TRUE)

dat1Collapsed=data.frame(collapse.object$group2row[,1], collapse.object$datETcollapsed)
colnames(dat1Collapsed)[1] = "Peptide_upper"
write.csv(dat1Collapsed, file = "40-087_Adipose_phase_4_peptide_QS_edit_collapsed.csv", row.names=F)

length(set1DATA[,1])
length(dat1Collapsed[,1])

sessionInfo = paste("session-info ",Sys.time(),".txt", sep="")
sessionInfo = gsub(":", "-", sessionInfo)
sink(sessionInfo)
sessionInfo()
sink()

