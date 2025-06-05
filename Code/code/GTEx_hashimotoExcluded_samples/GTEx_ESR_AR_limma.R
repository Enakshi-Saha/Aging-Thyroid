set.seed(1)
#if (!requireNamespace("BiocManager", quietly = TRUE))   
#  install.packages("BiocManager",repos = "http://cran.us.r-project.org")  
#BiocManager::install("fgsea")
#BiocManager::install("limma")
#BiocManager::install("Biobase")
#install.packages("ggplot2")
#install.packages("igraph")
library(dplyr)
library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)
library(fgsea)
library(readr)
library(biomaRt)
library(edgeR)

# Get disease sample IDs
GTEx_thyroid_diseaseList <- readxl::read_excel("/home/ubuntu/GTEx_thyroid_diseaseList.xlsx")
hashimoto_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("hashimoto", GTEx_thyroid_diseaseList$`Pathology Categories`)]
goiter_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("goiter", GTEx_thyroid_diseaseList$`Pathology Categories`)]
adenoma_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("adenoma", GTEx_thyroid_diseaseList$`Pathology Categories`)]

# Get portal data and use only gtex portal samples
gtex.portal_thyroid <- read.csv("/home/ubuntu/gtex.portal_thyroid.txt", sep="")

##############################
bonobo = data.frame(fread("/home/esaha/Aging_thyroid/data/network/BonoboPanda_GTEx_thyroid_significantGene_subnetwork.txt"))
genes = bonobo$gene
TFs = bonobo$TF
bonobo = bonobo[,match(colnames(gtex.portal_thyroid), colnames(bonobo))]
bonobo_indices =  unlist(lapply(strsplit(colnames(bonobo), split = ".", fixed = T), function(x){paste(paste(x, collapse = "-"), 1, sep = ".")}))
bonobo$TF = TFs

# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)

# Get phenotypic data
phenotypes <- read.csv("~/BONOBO/data/GTEx_thyroid/thyroid_phenotypes.txt", sep="")
GTEx_phenotypes = as.data.frame(cbind(phenotypes$SUBJID, phenotypes$SEX, phenotypes$AGE, phenotypes$RACE, 
                                      phenotypes$BMI, phenotypes$TRISCH, 
                                      phenotypes$gtex.smrin, phenotypes$gtex.smnabtcht,
                                      phenotypes$MHSMKSTS))
rownames(GTEx_phenotypes) = rownames(phenotypes)
colnames(GTEx_phenotypes) = c("subject_ID", "gender", "age", "race", "bmi", "ischemic_time", "rna_degrad", "batch", "smoking")
GTEx_phenotypes =  GTEx_phenotypes[match(bonobo_indices, rownames(GTEx_phenotypes)),]

# Define disease status: normal, hashimoto, goiter
disease = rep("normal", nrow(GTEx_phenotypes))
disease[which(GTEx_phenotypes$subject_ID %in% hashimoto_samples)] = "hashimoto"
disease[which(GTEx_phenotypes$subject_ID %in% goiter_samples)] = "goiter"
disease[which(GTEx_phenotypes$subject_ID %in% adenoma_samples)] = "adenoma"

disease = factor(disease, levels = c("normal", "hashimoto", "goiter", "adenoma"))

bonobo = bonobo[,which(disease == "normal")]

# Extract ESR and AR subnetworks
ESR1_net = bonobo[which(TFs == "ESR1"),]
rownames(ESR1_net) = genes[which(TFs == "ESR1")]
ESR2_net = bonobo[which(TFs == "ESR2"),]
rownames(ESR2_net) = genes[which(TFs == "ESR2")]
AR_net = bonobo[which(TFs == "AR"),]
rownames(AR_net) = genes[which(TFs == "AR")]
PGR_net = bonobo[which(TFs == "PGR"),]
rownames(PGR_net) = genes[which(TFs == "PGR")]


# Define the covariates: gender, age, race etc
gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
gender = factor(gender, levels = c("FEMALE", "MALE"))

age <- as.numeric(as.character(GTEx_phenotypes$age))
age[which(is.na(age))] <- mean(age,na.rm=TRUE)

# Categorize race 1, 4, 99 as single class "others", 2: "black", 3: "white"
race <- as.numeric(as.character(GTEx_phenotypes$race))
race[which(race != 2 & race != 3)] = "others"
race[which(race == 2)] = "black or african american"
race[which(race == 3)] = "white"
race <- as.factor(race)

bmi <- as.numeric(as.character(GTEx_phenotypes$bmi))
bmi[which(is.na(bmi))] <- mean(bmi,na.rm=TRUE)

ischemic_time = GTEx_phenotypes$ischemic_time
hour_minute = strsplit(ischemic_time, split=" ")
ischemic_timeH = sapply(hour_minute, function(x){as.numeric(as.character(x[1]))+as.numeric(as.character(x[1]))/60})

rna_degrad = as.numeric(as.character(GTEx_phenotypes$rna_degrad))
rna_degrad[which(is.na(rna_degrad))] <- mean(rna_degrad,na.rm=TRUE)

batch_effect = as.factor(GTEx_phenotypes$batch)

smoking_status = GTEx_phenotypes$smoking
smoking_status[which(is.na(smoking_status))] = "Unknown"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

###### Fit limma model ######
#### Design Matrix
design = model.matrix(~ age + gender + race + bmi + ischemic_timeH + rna_degrad + batch_effect + smoking_status)
# gender = factor(gender, levels = c("MALE", "FEMALE"))
# design = model.matrix(~ age*gender + race + bmi + ischemic_timeH + rna_degrad + batch_effect + smoking_status)
design = design[which(disease == "normal"),]
fit_ESR1 <- lmFit(ESR1_net, design)
fit_ESR1 <- eBayes(fit_ESR1)
fit_ESR2 <- lmFit(ESR2_net, design)
fit_ESR2 <- eBayes(fit_ESR2)
fit_AR <- lmFit(AR_net, design)
fit_AR <- eBayes(fit_AR)
fit_PGR <- lmFit(PGR_net, design)
fit_PGR <- eBayes(fit_PGR)

# Save table for age effect in female (base gender level)
tb_ESR1 = topTable(fit_ESR1,coef="age",number=Inf)
tb_ESR2 = topTable(fit_ESR2,coef="age",number=Inf)
tb_AR = topTable(fit_AR,coef="age",number=Inf)
tb_PGR = topTable(fit_PGR,coef="age",number=Inf)

gene_id = gsub("\\..*", "", fdata.GTEx$gene_id)
tb_ESR1$gene_name = fdata.GTEx$gene_name[match(rownames(tb_ESR1),gene_id)]
tb_ESR2$gene_name = fdata.GTEx$gene_name[match(rownames(tb_ESR2),gene_id)]
tb_AR$gene_name = fdata.GTEx$gene_name[match(rownames(tb_AR),gene_id)]
tb_PGR$gene_name = fdata.GTEx$gene_name[match(rownames(tb_PGR),gene_id)]

# write.table(tb_ESR1, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_aging_ESR1.txt")
# write.table(tb_ESR2, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_aging_ESR2.txt")
# write.table(tb_AR, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_aging_AR.txt")
# write.table(tb_PGR, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_aging_PGR.txt")

# write.table(tb_ESR1, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_M_aging_ESR1.txt")
# write.table(tb_ESR2, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_M_aging_ESR2.txt")
# write.table(tb_AR, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_M_aging_AR.txt")
# write.table(tb_PGR, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_M_aging_PGR.txt")

# write.table(tb_ESR1, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_F_aging_ESR1.txt")
# write.table(tb_ESR2, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_F_aging_ESR2.txt")
# write.table(tb_AR, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_F_aging_AR.txt")
# write.table(tb_PGR, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/ESR_AR/limma_GTEx_F_aging_PGR.txt")

# Number of significant genes
sum(tb_ESR1$P.Val < 0.05)
sum(tb_ESR2$P.Val < 0.05)
sum(tb_AR$P.Val < 0.05)
sum(tb_PGR$P.Val < 0.05)
unique(c(tb_ESR1$gene_name[tb_ESR1$P.Val < 0.05], 
         tb_ESR2$gene_name[tb_ESR2$P.Val < 0.05], tb_AR$gene_name[tb_AR$P.Val < 0.05]))


# Number of significant genes after FDR correction
sum(tb_ESR1$adj.P.Val < 0.05)
sum(tb_ESR2$adj.P.Val < 0.05)
sum(tb_AR$adj.P.Val < 0.05)
sum(tb_PGR$adj.P.Val < 0.05)
unique(c(tb_ESR1$gene_name[tb_ESR1$adj.P.Val < 0.05], 
         tb_ESR2$gene_name[tb_ESR2$adj.P.Val < 0.05], 
         tb_AR$gene_name[tb_AR$adj.P.Val < 0.05], 
         tb_PGR$gene_name[tb_PGR$adj.P.Val < 0.05]))

