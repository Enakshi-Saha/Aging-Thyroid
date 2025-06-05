set.seed(1)
#if (!requireNamespace("BiocManager", quietly = TRUE))   
#  install.packages("BiocManager",repos = "http://cran.us.r-project.org")  

library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)
library(fgsea)
library(readr)
library(biomaRt)
library(edgeR)
library(dplyr)

# Get disease sample IDs
GTEx_thyroid_diseaseList <- readxl::read_excel("/home/ubuntu/GTEx_thyroid_diseaseList.xlsx")
hashimoto_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("hashimoto", GTEx_thyroid_diseaseList$`Pathology Categories`)]
goiter_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("goiter", GTEx_thyroid_diseaseList$`Pathology Categories`)]
adenoma_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("adenoma", GTEx_thyroid_diseaseList$`Pathology Categories`)]

gtex.portal_thyroid <- read.csv("/home/ubuntu/gtex.portal_thyroid.txt", sep="")

##############################
bonobo = data.frame(fread("/home/esaha/Aging_thyroid/data/network/BonoboPanda_GTEx_thyroid_significantGene_subnetwork.txt"))
genes = bonobo$gene
TFs = bonobo$TF
bonobo = bonobo[,match(colnames(gtex.portal_thyroid), colnames(bonobo))]
bonobo_indices =  unlist(lapply(strsplit(colnames(bonobo), split = ".", fixed = T), function(x){paste(paste(x, collapse = "-"), 1, sep = ".")}))

# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)
gene_name = gsub("\\..*", "", fdata.GTEx$gene_id)
genename = fdata.GTEx$gene_name[match(gene_name,genes)]

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

bonobo = bonobo[,which(disease == "normal")]

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

# ESR1
ESR1 = colMeans(bonobo[which(TFs == "ESR1"),])
ESR2 = colMeans(bonobo[which(TFs == "ESR2"),])
AR = colMeans(bonobo[which(TFs == "AR"),])
SP1 = colMeans(bonobo[which(TFs == "SP1"),])
TP53 = colMeans(bonobo[which(TFs == "TP53"),])
FOXE1 = colMeans(bonobo[which(TFs == "FOXE1"),])
PAX8 = colMeans(bonobo[which(TFs == "PAX8"),])
HHEX = colMeans(bonobo[which(TFs == "HHEX"),])
RUNX2 = colMeans(bonobo[which(TFs == "RUNX2"),])
PGR = colMeans(bonobo[which(TFs == "PGR"),])

fulldata = data.frame(age)
fulldata$smoking_status = as.character(smoking_status)
fulldata$sex = gender
fulldata$race = race

fulldata$ESR1 = ESR1
fulldata$ESR2 = ESR2
fulldata$AR = AR
fulldata$SP1 = SP1
fulldata$TP53 = TP53
fulldata$FOXE1 = FOXE1
fulldata$PAX8 = PAX8
fulldata$HHEX = HHEX
fulldata$RUNX2 = RUNX2
fulldata$PGR = PGR

age_cuts = quantile(age, probs = seq(0, 1, 0.05))
age_cuts[1] = age_cuts[1]-1
age_cuts[length(age_cuts)] = age_cuts[length(age_cuts)] + 1
mean_interval = function(x){
  y = rep(0, (length(x)-1))
  for (i in 2:length(x)){
    y[i-1] = (x[i-1] + x[i])/2
  }
  return(y)
}
fulldata$age_group = cut(fulldata$age,
                         breaks=age_cuts,
                         labels=mean_interval(age_cuts))
fulldata$age_group = as.numeric(as.character(fulldata$age_group))
head(fulldata)

############################
# Age trajectory by sex

# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_targeting = mean(ESR1))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_targeting, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: ESR1 targeting of aging genes with age ",
                                                               y="gene score", x = "age")

############################
# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_targeting = mean(ESR2))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_targeting, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: ESR2 targeting of aging genes with age ",
                                                               y="gene score", x = "age")
############################
# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_targeting = mean(AR))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_targeting, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: AR targeting of aging genes with age ",
                                                               y="gene score", x = "age")

############################
# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_targeting = mean(SP1))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_targeting, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: SP1 targeting of aging genes with age ",
                                                               y="gene score", x = "age")

############################
# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_targeting = mean(TP53))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_targeting, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: TP53 targeting of aging genes with age ",
                                                               y="gene score", x = "age")

############################
# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_targeting = mean(RUNX2))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_targeting, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: RUNX2 targeting of aging genes with age ",
                                                               y="gene score", x = "age")
############################
# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_targeting = mean(PGR))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_targeting, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: PGR targeting of aging genes with age ",
                                                               y="gene score", x = "age")

