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

# Load immune cell composition
GTEx_xcell <- read.csv("/home/esaha/Aging_thyroid/data/immuneInfiltration/GTEx_xcell.txt", sep="")
rownames(GTEx_xcell) = GTEx_xcell$cell_type
GTEx_xcell = GTEx_xcell[,-1]

# Get GTEx sample indices
sample_indices = unlist(lapply(strsplit(colnames(GTEx_xcell), split = ".", fixed = T), function(x){paste(paste(x[-length(x)], collapse = "-"), 1, sep = ".")}))

# Get phenotypic data
phenotypes <- read.csv("~/BONOBO/data/GTEx_thyroid/thyroid_phenotypes.txt", sep="")
GTEx_phenotypes = as.data.frame(cbind(phenotypes$SEX, phenotypes$AGE, phenotypes$RACE, 
                                      phenotypes$BMI, phenotypes$TRISCH, 
                                      phenotypes$gtex.smrin, phenotypes$gtex.smnabtcht,
                                      phenotypes$MHSMKSTS))
rownames(GTEx_phenotypes) = rownames(phenotypes)
colnames(GTEx_phenotypes) = c("gender", "age", "race", "bmi", "ischemic_time", "rna_degrad", "batch", "smoking")
GTEx_phenotypes =  GTEx_phenotypes[match(sample_indices, rownames(GTEx_phenotypes)),]

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

# Linear model coefficient of each cell type: age coefficient
# gender = factor(gender, levels = c("FEMALE", "MALE"))
GTEx_lm = apply(GTEx_xcell, MARGIN = 1, function(x){summary(lm(as.numeric(x) ~ age + gender + race + ischemic_timeH + smoking_status + age*gender))$coefficients[2,c(3,4)]})
rownames(GTEx_lm) = c("age_coefficient", "p-value")
GTEx_lm = t(GTEx_lm)
GTEx_lm

# Save immune deconv table
immune_table = GTEx_lm
immune_table = format(immune_table, digits=4)
colnames(immune_table)[1] = "age coefficient"
write.csv(immune_table, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_thyroid_male_immuneInfiltration_aging.csv", row.names = T, quote = F)
# write.csv(immune_table, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_thyroid_female_immuneInfiltration_aging.csv", row.names = T, quote = F)

######################################
# Linear model coefficient of each cell type: sex coefficient
# gender = factor(gender, levels = c("FEMALE", "MALE"))
GTEx_lm = apply(GTEx_xcell, MARGIN = 1, function(x){summary(lm(as.numeric(x) ~ age + gender + race + ischemic_timeH + smoking_status + age*gender))$coefficients[3,c(3,4)]})
rownames(GTEx_lm) = c("sex_coefficientFEMALE", "p-value")
GTEx_lm = t(GTEx_lm)
GTEx_lm

# Save immune deconv table
immune_table = GTEx_lm
immune_table = format(immune_table, digits=4)
colnames(immune_table)[1] = "sex coefficient"
write.csv(immune_table, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_thyroid_immuneInfiltration_sexCoefficient.csv", row.names = T, quote = F)
