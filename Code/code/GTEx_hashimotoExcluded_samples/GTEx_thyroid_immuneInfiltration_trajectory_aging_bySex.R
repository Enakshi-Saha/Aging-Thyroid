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

# Define disease status: normal, hashimoto, goiter
disease = rep("normal", nrow(GTEx_phenotypes))
disease[which(GTEx_phenotypes$subject_ID %in% hashimoto_samples)] = "hashimoto"
disease[which(GTEx_phenotypes$subject_ID %in% goiter_samples)] = "goiter"
disease[which(GTEx_phenotypes$subject_ID %in% adenoma_samples)] = "adenoma"

disease = factor(disease, levels = c("normal", "hashimoto", "goiter", "adenoma"))
GTEx_xcell = GTEx_xcell[,which(disease == "normal")]

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

fulldata = data.frame(age)
fulldata$sex = gender

fulldata$CD4 = as.numeric(GTEx_xcell[which(rownames(GTEx_xcell) == "T cell CD4+ Th2"),])
fulldata$TgammaDelta = as.numeric(GTEx_xcell[which(rownames(GTEx_xcell) == "T cell gamma delta"),])
fulldata$CMP = as.numeric(GTEx_xcell[which(rownames(GTEx_xcell) == "Common myeloid progenitor"),])
fulldata$endothelial = as.numeric(GTEx_xcell[which(rownames(GTEx_xcell) == "Endothelial cell"),])

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
  summarise(mean_concentration = mean(CD4))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_concentration, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: Proportion of CD4+ Th2 T cells with age ",
                                                             y="mean cell proportion", x = "age")
# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_concentration = mean(TgammaDelta))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_concentration, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: Proportion of T cells gamma delta with age ",
                                                               y="mean cell proportion", x = "age")
# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_concentration = mean(CMP))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_concentration, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: Proportion of Common Myeloid Progenitor cells with age ",
                                                             y="mean cell proportion", x = "age")
# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(mean_concentration = mean(endothelial))
head(fulldata_aggregated)

fulldata_aggregated$sex = factor(fulldata_aggregated$sex, levels = c("FEMALE", "MALE"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_concentration, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title="GTEx: Proportion of endothelial cells with age ",
                                                             y="mean cell proportion", x = "age")
