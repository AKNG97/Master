library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)

AML_Normal_BM <- readRDS("AML_Normal_BM.rds")
ALL <- readRDS("rnas_raw_TAP2.RDS")
MM_raw <- readRDS("rnas_raw_MM.RDS")
AMLNB <- readRDS("rnas_raw_AMLNB.RDS")

B_ALL <- ALL[ , (ALL$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia" & ALL$sample_type == "Primary Blood Derived Cancer - Bone Marrow") | (ALL$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia" & ALL$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow")]
T_ALL <- ALL[ , ALL$primary_diagnosis == "T lymphoblastic leukemia/lymphoma" & ALL$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
MM <- MM_raw[ , MM_raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow" | MM_raw$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"]
AML_BM <- AML_Normal_BM[, AML_Normal_BM$sample_type == "Primary Blood Derived Cancer - Bone Marrow" | AML_Normal_BM$sample_type =="Recurrent Blood Derived Cancer - Bone Marrow"] 
Normal_BoneMarrow <- AML_Normal_BM[ , AML_Normal_BM$sample_type == "Bone Marrow Normal"]

factorsBALL <- data.frame(Group = "BALL", Sample = colnames(B_ALL))
factorsTALL <- data.frame(Group = "TALL", Sample = colnames(T_ALL))
factorsMM <- data.frame(Group = "MM", Sample = colnames(MM))
factorsAML <- data.frame(Group = "AML", Sample = colnames(AML_BM))
factorsNormalBM <- data.frame(Group = "NormalBM", Sample = colnames(Normal_BoneMarrow))

Metastisic <- rnas1[ , rnas1$sample_type == "Metastatic"]
factorsRichard_1 <- data.frame(Group = "Metastasico", Sample = colnames(Metastatic))

rownames(B_ALL) <- rowData(B_ALL)$external_gene_name
rownames(T_ALL) <- rowData(T_ALL)$external_gene_name
rownames(MM) <- rowData(MM)$external_gene_name
rownames(AML_BM) <- rowData(AML_BM)$external_gene_name
rownames(Normal_BoneMarrow) <- rowData(Normal_BoneMarrow)$external_gene_name

rnas <- cbind(assay(B_ALL), assay(T_ALL), assay(MM), assay(AML_BM), assay(Normal_BoneMarrow))
factors <- rbind(factorsBALL, factorsTALL, factorsMM, factorsAML, factorsNormalBM)
rownames(factors) <- factors$Sample
Ready_factors <- as.data.frame(factors$Group)

dim(rnas)

rnas <- rnas[!duplicated(rownames(rnas)),]
dim(rnas)

rnas_before <- rnas

dataFilt <- TCGAanalyze_Filtering(tabDF = rnas,
                                  method = "quantile",
                                  qnt.cut = 0.25)
threshold <- round(dim(rnas)[2]/2)
ridx <- rowSums(dataFilt == 0) <= threshold
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))
rnas <- rnas[rownames(rnas) %in% rownames(dataFilt), ]
print(dim(rnas))

annot<-read.delim(file="mart_export.txt", sep="\t")
names(annot)<-c("Gene.name", "Chr", "Start", "End", "GC", "Type", "ensembl_gene_id")
annot$Length <- abs(annot$End - annot$Start)
inter <- intersect(rownames(rnas), annot$Gene.name)
dim(inter)
rnas1 <- rnas[rownames(rnas) %in% inter,]
dim(rnas1)
annot <- annot[annot$Gene.name %in% inter,]
dim(annot)
annot <- annot[!duplicated(annot$Gene.name),]
dim(annot)

ln.data <- withinLaneNormalization(rnas1, annot$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annot$GC, which = "full")
Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") 
norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = Ready_factors)
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = TRUE)
rnas_after <- exprs(mydata2corr1)

mydata.before = NOISeq::readData(rnas1, factors = Ready_factors)
myPCA.before = dat(mydata.before, type = "PCA")
png(filename="PCA_before.png")
explo.plot(myPCA.before)

mydata.after = readData(rnas_after, factors = Ready_factors)
myPCA.after = dat(mydata.after, type = "PCA")
png(filename="PCA_after.png")
explo.plot(myPCA.after)

rnas.pca <- prcomp(assay(rnas1), center = TRUE,scale. = TRUE)
summary(rnas.pca)
autoplot(rnas.pca, data = t(assay(rnas1)), colour = as.data.frame(rnas1$grupo))

saveRDS(rnas1, file="rnas_after.RDS")

Everyth <- readRDS("rnas_after_Everything.RDS")

after.pca <- prcomp(t(rnas_after),center = TRUE,scale. = TRUE)
summary(after.pca)
library(devtools)
install_github("vqv/ggbiplot")
install.packages("devtools")
library(ggbiplot)
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=factors$Group)

saveRDS(rnas_after, file="Cancers_vs_Normal.RDS")
saveRDS(factors, file="factors_Cancers_vs_Normal.RDS")

BALL_Norm <- rnas_after[, factors$Group=="BALL"]
TALL_Norm <- rnas_after[, factors$Group=="TALL"]
MM_Norm <- rnas_after[, factors$Group=="MM"]
AML_Norm <- rnas_after[, factors$Group=="AML"]
NormalBM_Norm <- rnas_after[, factors$Group=="NormalBM"]

Aracne_BALL_Norm <- cbind(rownames(BALL_Norm), BALL_Norm)
Aracne_TALL_Norm <- cbind(rownames(TALL_Norm), TALL_Norm)
Aracne_MM_Norm <- cbind(rownames(MM_Norm), MM_Norm)
Aracne_AML_Norm <- cbind(rownames(AML_Norm), AML_Norm)
Aracne_NormalBM_Norm <- cbind(rownames(NormalBM_Norm), NormalBM_Norm)

colnames(Aracne_BALL_Norm)[1] <- "gene"
colnames(Aracne_TALL_Norm)[1] <- "gene"
colnames(Aracne_MM_Norm)[1] <- "gene"
colnames(Aracne_AML_Norm)[1] <- "gene"
colnames(Aracne_NormalBM_Norm)[1] <- "gene"

write.table(Aracne_BALL_Norm, file = "All_BALL_Norm.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(Aracne_TALL_Norm, file = "All_TALL_Norm.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(Aracne_MM_Norm, file = "All_MM_Norm.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(Aracne_AML_Norm, file = "All_AML_Norm.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(Aracne_NormalBM_Norm, file = "All_NormalBM_Norm.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

check <- read.table("All_NormalBM_Norm.tsv")
