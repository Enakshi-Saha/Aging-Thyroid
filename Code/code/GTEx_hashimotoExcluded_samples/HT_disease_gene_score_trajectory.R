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
library(gghalves)

# Get disease sample IDs
GTEx_thyroid_diseaseList <- readxl::read_excel("/home/ubuntu/GTEx_thyroid_diseaseList.xlsx")
hashimoto_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("hashimoto", GTEx_thyroid_diseaseList$`Pathology Categories`)]
goiter_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("goiter", GTEx_thyroid_diseaseList$`Pathology Categories`)]
adenoma_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("adenoma", GTEx_thyroid_diseaseList$`Pathology Categories`)]

##############################
bonobo = data.frame(fread("/home/esaha/BONOBO/network/BonoboPanda_GTEx_thyroid.txt"))
genes = bonobo$V1
bonobo = bonobo[,-1]
bonobo_indices =  unlist(lapply(strsplit(colnames(bonobo), split = ".", fixed = T), function(x){paste(x[1:2], collapse = "-")}))

# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)
gene_id = gsub("\\..*", "", fdata.GTEx$gene_id)
gene_name = fdata.GTEx$gene_name[match(genes,gene_id)]

# Get disease gene list
diseaseEffect_hashimoto1 = data.frame(read.csv("~/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_thyroid_diseaseEffect_hashimoto.txt", sep=""))
diseaseEffect_hashimoto2 = data.frame(read.csv("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE29315/limma_GSE29315_diseaseEffect_hashimoto thyroiditis.txt", sep=""))

disease_gene_increasing1 = diseaseEffect_hashimoto1$gene_name[which(diseaseEffect_hashimoto1$P.Val < 0.05 & diseaseEffect_hashimoto1$t > 0)]
disease_gene_decreasing1 = diseaseEffect_hashimoto1$gene_name[which(diseaseEffect_hashimoto1$P.Val < 0.05 & diseaseEffect_hashimoto1$t < 0)]

disease_gene_increasing2 = rownames(diseaseEffect_hashimoto2)[which(diseaseEffect_hashimoto2$P.Val < 0.05 & diseaseEffect_hashimoto2$t > 0)]
disease_gene_decreasing2 = rownames(diseaseEffect_hashimoto2)[which(diseaseEffect_hashimoto2$P.Val < 0.05 & diseaseEffect_hashimoto2$t < 0)]

disease_gene_increasing = intersect(disease_gene_increasing1, disease_gene_increasing2)
disease_gene_decreasing = intersect(disease_gene_decreasing1, disease_gene_decreasing2)

# Get disease gene indegrees
score_disease_gene_increasing = colMeans(bonobo[which(gene_name %in% disease_gene_increasing),])
score_disease_gene_decreasing = colMeans(bonobo[which(gene_name %in% disease_gene_decreasing),])

# Get phenotypic data
phenotypes <- read.csv("~/BONOBO/data/GTEx_thyroid/thyroid_phenotypes.txt", sep="")
GTEx_phenotypes = as.data.frame(cbind(phenotypes$SUBJID,phenotypes$SEX, phenotypes$AGE, phenotypes$RACE, 
                                      phenotypes$BMI, phenotypes$TRISCH, 
                                      phenotypes$gtex.smrin, phenotypes$gtex.smnabtcht,
                                      phenotypes$MHSMKSTS))
rownames(GTEx_phenotypes) = rownames(phenotypes)
colnames(GTEx_phenotypes) = c("subject_ID", "gender", "age", "race", "bmi", "ischemic_time", "rna_degrad", "batch", "smoking")
GTEx_phenotypes =  GTEx_phenotypes[match(bonobo_indices, GTEx_phenotypes$subject_ID),]

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

# Age trajectory
fulldata = data.frame(age)
fulldata$sex = factor(gender, levels = c("FEMALE", "MALE"))

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
fulldata$increasing_score = score_disease_gene_increasing
fulldata$decreasing_score = score_disease_gene_decreasing

fulldata = fulldata[which(disease == "normal"),]

# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, sex)  %>%
  summarise(across(everything(), list(mean)))
head(fulldata_aggregated)

head(fulldata)

ggplot(fulldata_aggregated, aes(x = age_group , y = increasing_score_1, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title=paste("GTEx: Age-trajectory of genes upregulated in HT"),
                                                               y="gene score", x = "age")
ggplot(fulldata_aggregated, aes(x = age_group , y = decreasing_score_1, color = sex)) +
  geom_point(size = 0.5) + geom_smooth(method = "auto") + labs(title=paste("GTEx: Age-trajectory of genes downregulated in HT"),
                                                               y="gene score", x = "age")

#################################################
# Boxplot by age group
fulldata$age_group = rep("56+", nrow(fulldata))
fulldata$age_group[which(fulldata$age<56)] = "<56"
ggplot(fulldata, aes(x = age_group , y = increasing_score, fill = sex)) +
  scale_fill_manual(values=c("purple1", "lightgreen")) +
  geom_half_boxplot(side = "left") + geom_half_violin(side = "right") + labs(title="Genes upregulated in HT",
                                                                             y="gene score", x = "age group") +
  scale_y_continuous(limits = quantile(fulldata$increasing_score, c(0.05, 0.95))) +
  theme(text=element_text(size=22))
ggplot(fulldata, aes(x = age_group , y = decreasing_score, fill = sex)) +
  scale_fill_manual(values=c("purple1", "lightgreen")) +
  geom_half_boxplot(side = "left") + geom_half_violin(side = "right") + labs(title="Genes downregulated in HT",
                                                                             y="gene score", x = "age group") +
  scale_y_continuous(limits = quantile(fulldata$decreasing_score, c(0.05, 0.95))) +
  theme(text=element_text(size=22))

####### one-sided t-test for each age group ##############
subdata = fulldata[which(fulldata$age_group == "56+"),]
# subdata = fulldata[which(fulldata$age_group == "<56"),]
t.test(subdata$increasing_score[which(subdata$sex == "MALE")], subdata$increasing_score[which(subdata$sex == "FEMALE")], alternative = "less") 
t.test(subdata$decreasing_score[which(subdata$sex == "MALE")], subdata$decreasing_score[which(subdata$sex == "FEMALE")], alternative = "greater") 


