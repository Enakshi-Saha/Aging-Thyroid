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

##############################
bonobo = data.frame(fread("/home/esaha/BONOBO/network/BonoboPanda_GTEx_thyroid.txt"))
genes = bonobo$V1
bonobo = bonobo[,-1]
bonobo_indices =  unlist(lapply(strsplit(colnames(bonobo), split = ".", fixed = T), function(x){paste(paste(x, collapse = "-"), 1, sep = ".")}))

# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)

# Get phenotypic data
phenotypes <- read.csv("~/BONOBO/data/GTEx_thyroid/thyroid_phenotypes.txt", sep="")
GTEx_phenotypes = as.data.frame(cbind(phenotypes$SUBJID,phenotypes$SEX, phenotypes$AGE, phenotypes$RACE, 
                                      phenotypes$BMI, phenotypes$TRISCH, 
                                      phenotypes$gtex.smrin, phenotypes$gtex.smnabtcht,
                                      phenotypes$MHSMKSTS))
rownames(GTEx_phenotypes) = rownames(phenotypes)
colnames(GTEx_phenotypes) = c("subject_ID", "gender", "age", "race", "bmi", "ischemic_time", "rna_degrad", "batch", "smoking")
GTEx_phenotypes =  GTEx_phenotypes[match(bonobo_indices, rownames(GTEx_phenotypes)),]

# Define the covariates: gender, age, race etc
gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
gender = factor(gender, levels = c("MALE", "FEMALE"))

# Select only female samples
bonobo = bonobo[,which(gender == "FEMALE")]
GTEx_phenotypes = GTEx_phenotypes[which(gender == "FEMALE"),]

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

batch_effect = as.factor(GTEx_phenotypes$batch)

smoking_status = GTEx_phenotypes$smoking
smoking_status[which(is.na(smoking_status))] = "Unknown"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

age_group = rep("below50", nrow(GTEx_phenotypes))
age_group[which(age>50)] = "above50"
age_group = factor(age_group, levels = c("below50", "above50"))
# age_group = factor(age_group, levels = c("above50", "below50"))

###### Fit limma model ######
indegree = bonobo
rownames(indegree) = genes
#### Design Matrix
design = model.matrix(~ age*age_group + race + bmi + ischemic_timeH + rna_degrad + batch_effect +
                        smoking_status)
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for disease effect: Hashimoto
tb = topTable(fit,coef="age",number=Inf)

gene_id = gsub("\\..*", "", fdata.GTEx$gene_id)
tb$chr = fdata.GTEx$seqnames[match(rownames(tb),gene_id)]
tb$gene_name = fdata.GTEx$gene_name[match(rownames(tb),gene_id)]

head(tb)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c5.go.bp.v2022.1.Hs.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c2.cp.reactome.v2022.1.Hs.symbols.gmt")

# Rank genes in limma table: Hashimoto
indegree_rank <- setNames(object=tb[,"t"], tb$gene_name)
head(indegree_rank)
fgseaRes <- fgsea(pathways, indegree_rank, minSize=15, maxSize=500)
head(fgseaRes)

save(fgseaRes, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_ageEffect_female_below50_GSEA.RData")

write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_thyroid_ageEffect_female_below50.txt")

# save(fgseaRes, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_ageEffect_female_above50_GSEA.RData")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_thyroid_ageEffect_female_above50.txt")
