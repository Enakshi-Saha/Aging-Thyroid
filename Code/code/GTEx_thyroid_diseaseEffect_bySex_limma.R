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
# gender = factor(gender, levels = c("FEMALE", "MALE"))

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
indegree = bonobo
rownames(indegree) = genes
#### Design Matrix
design = model.matrix(~ age + gender + race + bmi + ischemic_timeH + rna_degrad + batch_effect +
                        smoking_status + disease*gender)
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for disease effect: Hashimoto
tb_h = topTable(fit,coef="diseasehashimoto",number=Inf)
tb_g = topTable(fit,coef="diseasegoiter",number=Inf)
tb_a = topTable(fit,coef="diseaseadenoma",number=Inf)

gene_id = gsub("\\..*", "", fdata.GTEx$gene_id)
tb_h$chr = fdata.GTEx$seqnames[match(rownames(tb_h),gene_id)]
tb_h$gene_name = fdata.GTEx$gene_name[match(rownames(tb_h),gene_id)]
tb_g$chr = fdata.GTEx$seqnames[match(rownames(tb_g),gene_id)]
tb_g$gene_name = fdata.GTEx$gene_name[match(rownames(tb_g),gene_id)]
tb_a$chr = fdata.GTEx$seqnames[match(rownames(tb_a),gene_id)]
tb_a$gene_name = fdata.GTEx$gene_name[match(rownames(tb_a),gene_id)]

head(tb_h)
head(tb_g)
head(tb_a)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c5.go.bp.v2022.1.Hs.symbols.gmt")
# pathways <- gmtPathways("/home/ubuntu/c2.cp.reactome.v2022.1.Hs.symbols.gmt")

# Rank genes in limma table: Hashimoto
indegree_rank_h <- setNames(object=tb_h[,"t"], tb_h$gene_name)
head(indegree_rank_h)
fgseaRes_h <- fgsea(pathways, indegree_rank_h, minSize=15, maxSize=500)
head(fgseaRes_h)

# Rank genes in limma table: Goiter
indegree_rank_g <- setNames(object=tb_g[,"t"], tb_g$gene_name)
head(indegree_rank_g)
fgseaRes_g <- fgsea(pathways, indegree_rank_g, minSize=15, maxSize=500)
head(fgseaRes_g)

# Rank genes in limma table: Adenoma
indegree_rank_a <- setNames(object=tb_a[,"t"], tb_h$gene_name)
head(indegree_rank_a)
fgseaRes_a <- fgsea(pathways, indegree_rank_a, minSize=15, maxSize=500)
head(fgseaRes_a)

save(fgseaRes_h, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_hashimoto_male_GSEA.RData")
save(fgseaRes_g, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_goiter_male_GSEA.RData")
save(fgseaRes_a, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_adenoma_male_GSEA.RData")

write.table(tb_h, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_thyroid_diseaseEffect_hashimoto_male.txt")
write.table(tb_g, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_thyroid_diseaseEffect_goiter_male.txt")
write.table(tb_a, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_thyroid_diseaseEffect_adenoma_male.txt")

# save(fgseaRes_h, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_hashimoto_female_GSEA.RData")
# save(fgseaRes_g, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_goiter_female_GSEA.RData")
# save(fgseaRes_a, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_adenoma_female_GSEA.RData")

# write.table(tb_h, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_thyroid_diseaseEffect_hashimoto_female.txt")
# write.table(tb_g, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_thyroid_diseaseEffect_goiter_female.txt")
# write.table(tb_a, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_thyroid_diseaseEffect_adenoma_female.txt")
