set.seed(1)
#if (!requireNamespace("BiocManager", quietly = TRUE))   
#  install.packages("BiocManager",repos = "http://cran.us.r-project.org")  
#BiocManager::install("fgsea")
#BiocManager::install("limma")
#BiocManager::install("Biobase")
#install.packages("ggplot2")
#install.packages("igraph")
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
bonobo = data.frame(fread("/home/esaha/Aging_thyroid/data/network/BonoboPanda_GTEx_thyroid_TFoutdegree.txt"))
genes = bonobo$V1
bonobo = bonobo[,match(colnames(gtex.portal_thyroid), colnames(bonobo))]
bonobo_indices =  unlist(lapply(strsplit(colnames(bonobo), split = ".", fixed = T), function(x){paste(paste(x, collapse = "-"), 1, sep = ".")}))

# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)


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
indegree = bonobo
rownames(indegree) = genes
#### Design Matrix
design = model.matrix(~ age*gender + gender + race + bmi + ischemic_timeH + rna_degrad + batch_effect +
                        smoking_status)
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for sex difference: genderMALE
tb = topTable(fit,coef="age",number=Inf)
colnames(design) = make.names(colnames(design))
contrast.matrix <- makeContrasts(ageFEMALE = age + age.genderFEMALE, levels=design)
colnames(fit$coefficients) = rownames(contrast.matrix)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# tb = topTable(fit2,coef="ageFEMALE",number=Inf)

head(tb)

# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_TF.txt")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_TF_male.txt")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_TF_female.txt")
