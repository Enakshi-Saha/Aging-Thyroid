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

# Compute TF out degree
TF_outdegree = bonobo %>% group_by(TF) %>% 
  summarise(across(everything(), sum),
            .groups = 'drop')  %>%
  as.data.frame()
# apply(bonobo,MARGIN = 2,function(x){newdata = data.frame(x); newdata$TF = TFs; aggregate(newdata$x, by=list(Category=newdata$TF), FUN=sum)$x})
rownames(TF_outdegree) = TF_outdegree$TF
TF_outdegree = TF_outdegree[,-1]

# Get phenotypic data
phenotypes <- read.csv("~/BONOBO/data/GTEx_thyroid/thyroid_phenotypes.txt", sep="")
GTEx_phenotypes = as.data.frame(cbind(phenotypes$SEX, phenotypes$AGE, phenotypes$RACE, 
                                      phenotypes$BMI, phenotypes$TRISCH, 
                                      phenotypes$gtex.smrin, phenotypes$gtex.smnabtcht,
                                      phenotypes$MHSMKSTS))
rownames(GTEx_phenotypes) = rownames(phenotypes)
colnames(GTEx_phenotypes) = c("gender", "age", "race", "bmi", "ischemic_time", "rna_degrad", "batch", "smoking")
GTEx_phenotypes =  GTEx_phenotypes[match(bonobo_indices, rownames(GTEx_phenotypes)),]

# Define disease status: normal, hashimoto, goiter
disease = rep("normal", nrow(GTEx_phenotypes))
disease[which(GTEx_phenotypes$subject_ID %in% hashimoto_samples)] = "hashimoto"
disease[which(GTEx_phenotypes$subject_ID %in% goiter_samples)] = "goiter"
disease[which(GTEx_phenotypes$subject_ID %in% adenoma_samples)] = "adenoma"

disease = factor(disease, levels = c("normal", "hashimoto", "goiter", "adenoma"))

# Define the covariates: gender, age, race etc
gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
gender = factor(gender, levels = c("MALE", "FEMALE"))

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
design = design[which(disease == "normal"),]
fit <- lmFit(TF_outdegree, design)
fit <- eBayes(fit)

# Save table for sex difference: genderMALE
tb_m = topTable(fit,coef="age",number=Inf)
colnames(design) = make.names(colnames(design))
contrast.matrix <- makeContrasts(ageFEMALE = age + age.genderFEMALE, levels=design)
colnames(fit$coefficients) = rownames(contrast.matrix)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tb_f = topTable(fit2,coef="ageFEMALE",number=Inf)

TF0 = unlist(lapply(strsplit(rownames(tb_m), split = "_"), function(x){x[1]}))
tb_m$TF = TF0
tb_m$chr_TF = fdata.GTEx$seqnames[match(TF0,fdata.GTEx$gene_name)]

TF0 = unlist(lapply(strsplit(rownames(tb_f), split = "_"), function(x){x[1]}))
tb_f$TF = TF0
tb_f$chr_TF = fdata.GTEx$seqnames[match(TF0,fdata.GTEx$gene_name)]


head(tb_m)
head(tb_f)

sigTF_f = tb_f$TF[which(tb_f$adj.P.Val<0.05)]

# Table for ESR and AR
tb_f[which(tb_f$TF %in% c("ESR1","ESR2","AR")),]
tb_m[which(tb_m$TF %in% c("ESR1","ESR2","AR")),]


