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

# Extract ESR and AR subnetworks
ESR1_net = bonobo[which(TFs == "ESR1"),-ncol(bonobo)]
rownames(ESR1_net) = genes[which(TFs == "ESR1")]
ESR2_net = bonobo[which(TFs == "ESR2"),-ncol(bonobo)]
rownames(ESR2_net) = genes[which(TFs == "ESR2")]
AR_net = bonobo[which(TFs == "AR"),-ncol(bonobo)]
rownames(AR_net) = genes[which(TFs == "AR")]


# Get phenotypic data
phenotypes <- read.csv("~/BONOBO/data/GTEx_thyroid/thyroid_phenotypes.txt", sep="")
GTEx_phenotypes = as.data.frame(cbind(phenotypes$SEX, phenotypes$AGE, phenotypes$RACE, 
                                      phenotypes$BMI, phenotypes$TRISCH, 
                                      phenotypes$gtex.smrin, phenotypes$gtex.smnabtcht,
                                      phenotypes$MHSMKSTS))
rownames(GTEx_phenotypes) = rownames(phenotypes)
colnames(GTEx_phenotypes) = c("gender", "age", "race", "bmi", "ischemic_time", "rna_degrad", "batch", "smoking")
GTEx_phenotypes =  GTEx_phenotypes[match(bonobo_indices, rownames(GTEx_phenotypes)),]

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
design = model.matrix(~ age*gender + gender + race + bmi + ischemic_timeH + rna_degrad + batch_effect +
                        smoking_status)
fit_ESR1 <- lmFit(ESR1_net, design)
fit_ESR1 <- eBayes(fit_ESR1)
fit_ESR2 <- lmFit(ESR2_net, design)
fit_ESR2 <- eBayes(fit_ESR2)
fit_AR <- lmFit(AR_net, design)
fit_AR <- eBayes(fit_AR)

# Save table for age effect in female (base gender level)
tb_ESR1 = topTable(fit_ESR1,coef="age",number=Inf)
tb_ESR2 = topTable(fit_ESR2,coef="age",number=Inf)
tb_AR = topTable(fit_AR,coef="age",number=Inf)

gene_id = gsub("\\..*", "", fdata.GTEx$gene_id)
tb_ESR1$gene_name = fdata.GTEx$gene_name[match(rownames(tb_ESR1),gene_id)]
tb_ESR2$gene_name = fdata.GTEx$gene_name[match(rownames(tb_ESR2),gene_id)]
tb_AR$gene_name = fdata.GTEx$gene_name[match(rownames(tb_AR),gene_id)]

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c5.go.bp.v2022.1.Hs.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c2.cp.reactome.v2022.1.Hs.symbols.gmt")

# Rank genes in limma table: Hashimoto
indegree_rank_ESR1 <- setNames(object=tb_ESR1[,"t"], tb_ESR1$gene_name)
head(indegree_rank_ESR1)
fgseaRes_ESR1 <- fgsea(pathways, indegree_rank_ESR1, minSize=15, maxSize=500)
head(fgseaRes_ESR1)




