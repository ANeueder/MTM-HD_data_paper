# start from variance stabilised transformed data

workingDir = "D:/Work/! projects/MTM/proteomics/03 - batch correction/adipose protein/01 - sex";
setwd(workingDir);
library("limma")
library("DESeq2")
# options(stringsAsFactors = FALSE);

dat1 <- read.csv("40-087_Adipose_phase_4_protein_QS_edit.csv", header=TRUE);
names(dat1)
rownames(dat1) <- dat1$PID_edit
dat1 <- dat1 [,-c(1:5)]
head(dat1)

coldata = read.csv("DecoderRing_MTM_proteomics_combined_adipose.csv", header=TRUE);
head(coldata)

colnames(coldata)[1] = "identifier"
rownames(coldata)<- coldata$identifier
head(coldata)

traitRows = match(colnames(dat1), coldata$identifier);
# column number for allTraits$identifier
datTraits = coldata[traitRows, -1];
rownames(datTraits) = coldata[traitRows, 1];
dim(coldata)
dim(datTraits)
# show that row names agree
table(rownames(datTraits)==colnames(dat1))
datTraits = cbind(rownames(datTraits), datTraits)
colnames(datTraits)[1] = "identifier"

# needs at least two columns including batch to correct for
coldata = datTraits[, c("identifier", "Site","Sex","Group","HD","MTM_BMI","MTM_age_months")];
dim(coldata)
names(coldata)
colData<-data.frame(coldata)
head(colData)

# one factor design:
f <- factor(colData$Group)
design <- model.matrix(~0+f)

colnames(design)
batch <- coldata$Sex
vsdCor <- removeBatchEffect(dat1, batch=batch, design=design)
# for continuous confounders, i.e. numeric variables
# vsdCor <- removeBatchEffect(dat1, covariates=batch, design=design)
head(vsdCor)


vsdCor1=cbind(rownames(vsdCor), vsdCor)
colnames(vsdCor1)[1] = "PID_edit"
write.table(vsdCor1, file = "40-087_Adipose_phase_4_protein_QS_edit_limma_sex_corrected.csv", row.names=F, sep=",")


sessionInfo = paste("session-info ",Sys.time(),".txt", sep="")
sessionInfo = gsub(":", "-", sessionInfo)
sink(sessionInfo)
sessionInfo()
sink()


