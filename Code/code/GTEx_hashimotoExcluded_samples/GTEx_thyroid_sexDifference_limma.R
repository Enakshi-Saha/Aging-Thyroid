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

# Get disease sample IDs
GTEx_thyroid_diseaseList <- readxl::read_excel("/home/ubuntu/GTEx_thyroid_diseaseList.xlsx")
hashimoto_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("hashimoto", GTEx_thyroid_diseaseList$`Pathology Categories`)]
goiter_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("goiter", GTEx_thyroid_diseaseList$`Pathology Categories`)]
adenoma_samples = GTEx_thyroid_diseaseList$`Subject ID`[grep("adenoma", GTEx_thyroid_diseaseList$`Pathology Categories`)]

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

batch_effect = as.factor(GTEx_phenotypes$batch)

smoking_status = GTEx_phenotypes$smoking
smoking_status[which(is.na(smoking_status))] = "Unknown"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

###### Fit limma model ######
indegree = bonobo[,which(disease == "normal")]
rownames(indegree) = genes
#### Design Matrix
design = model.matrix(~ age + gender + race + bmi + ischemic_timeH + rna_degrad + batch_effect +
                        smoking_status)
design = design[which(disease == "normal"),]
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for age effect
# Save table for age effect
tb = topTable(fit,coef="genderFEMALE",number=Inf)
colnames(design) = make.names(colnames(design))

gene_id = gsub("\\..*", "", fdata.GTEx$gene_id)
tb$chr = fdata.GTEx$seqnames[match(rownames(tb),gene_id)]
tb$gene_name = fdata.GTEx$gene_name[match(rownames(tb),gene_id)]

head(tb)


# Rank genes in limma table
indegree_rank <- setNames(object=tb[,"t"], tb$gene_name)
head(indegree_rank)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c5.go.bp.v2022.1.Hs.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c2.cp.reactome.v2022.1.Hs.symbols.gmt")

fgseaRes <- fgsea(pathways, indegree_rank, minSize=15, maxSize=500)
head(fgseaRes)

# save(fgseaRes, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/GTEx_Thyroid_sexDifference_GSEA.RData")

# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/limma_GTEx_sexDifference.txt")
