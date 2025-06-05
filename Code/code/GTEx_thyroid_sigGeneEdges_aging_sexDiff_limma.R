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
bonobo = data.frame(fread("/home/esaha/Aging_thyroid/data/network/BonoboPanda_GTEx_thyroid_significantGene_subnetwork.txt"))
genes = bonobo$gene
TFs = bonobo$TF
bonobo = bonobo[,match(colnames(gtex.portal_thyroid), colnames(bonobo))]
bonobo_indices =  unlist(lapply(strsplit(colnames(bonobo), split = ".", fixed = T), function(x){paste(paste(x, collapse = "-"), 1, sep = ".")}))

# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)

# Discard TFs on the Y chromosome
### Get names of Y genes 
Y_genes.GTEx = fdata.GTEx$gene_name[which(fdata.GTEx$seqnames == "chrY")]

bonobo = bonobo[-which(TFs %in% Y_genes.GTEx),]
genes = genes[-which(TFs %in% Y_genes.GTEx)]
TFs = TFs[-which(TFs %in% Y_genes.GTEx)]

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
rownames(indegree) = paste0(TFs, "_" , genes)
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

gene_id = gsub("\\..*", "", fdata.GTEx$gene_id)
TF0 = unlist(lapply(strsplit(rownames(tb), split = "_"), function(x){x[1]}))
gene0 = unlist(lapply(strsplit(rownames(tb), split = "_"), function(x){x[2]}))
tb$gene_name = fdata.GTEx$gene_name[match(gene0,gene_id)]
tb$TF = TF0
tb$chr_gene = fdata.GTEx$seqnames[match(gene0,gene_id)]
tb$chr_TF = fdata.GTEx$seqnames[match(TF0,fdata.GTEx$gene_name)]

head(tb)

# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_sigGeneEdges.txt")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_sigGeneEdges_male.txt")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_sigGeneEdges_female.txt")

########################################
library(data.table)
# Get sex-specific significant aging genes
agingGenes_male_increasing <- read.table("~/Aging_thyroid/results/limma_GSEA_xcell_tables/agingGenes_male_increasing.txt", quote="\"", comment.char="")$V1
agingGenes_male_decreasing <- read.table("~/Aging_thyroid/results/limma_GSEA_xcell_tables/agingGenes_male_decreasing.txt", quote="\"", comment.char="")$V1
agingGenes_female_increasing <- read.table("~/Aging_thyroid/results/limma_GSEA_xcell_tables/agingGenes_female_increasing.txt", quote="\"", comment.char="")$V1
agingGenes_female_decreasing <- read.table("~/Aging_thyroid/results/limma_GSEA_xcell_tables/agingGenes_female_decreasing.txt", quote="\"", comment.char="")$V1

agingGene_male = c(agingGenes_male_increasing, agingGenes_male_decreasing)
agingGene_female = c(agingGenes_female_increasing, agingGenes_male_decreasing)

# Get limma results on edges
maleEdges = data.frame(fread("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_sigGeneEdges_male.txt"))
femaleEdges = data.frame(fread("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GTEx_aging_sigGeneEdges_female.txt"))

# Significant edges connected to sex-specific aging genes
maleEdges = maleEdges[which(maleEdges$gene_name %in% agingGene_male),]
femaleEdges = femaleEdges[which(femaleEdges$gene_name %in% agingGene_female),]

# Significant edges on Chr 1-22 or X
maleSigEdges = maleEdges[which(maleEdges$P.Value <0.05 & substr(maleEdges$chr_TF, start=1, stop=3) == "chr"),]
femaleSigEdges = femaleEdges[which(femaleEdges$P.Value <0.05 & substr(femaleEdges$chr_TF, start=1, stop=3) == "chr"),]

# Chromosomal distribution of TFs connected to significant edges
ng = 100
male_tab = table(maleSigEdges$chr_TF[1:ng])/length(maleSigEdges$chr_TF[1:ng])*100
female_tab = table(femaleSigEdges$chr_TF[1:ng])/length(femaleSigEdges$chr_TF[1:ng])*100

fulldata = data.frame(percentage=c(male_tab, female_tab))
fulldata$chromosome = factor(c(names(male_tab), names(female_tab)))
fulldata$sex = factor(c(rep("Male", length(male_tab)), rep("Female", length(female_tab))), levels = c("Female", "Male"))
head(fulldata)

# barplot to show position of most significant TFs
library(ggplot2)
maintitle= paste("Position of TFs corresponding to top", ng, "edges\nconnected to significant aging genes", sep = " ")
ggplot(data=fulldata, aes(y=chromosome, x=percentage, fill=sex)) +
  geom_bar(stat="identity", position = "dodge")+theme_minimal() + ggtitle(maintitle)

# Top TFs with most significant edges
table(maleSigEdges$TF[1:100])

table(femaleSigEdges$TF[1:100])

# Check common and distinct TFs
TF_male = maleSigEdges$TF[1:100]
TF_female = femaleSigEdges$TF[1:100]

intersect(TF_male, TF_female)

setdiff(TF_male, TF_female)

setdiff(TF_female, TF_male)

################################################
# Barplot of effect sizes by chromosome
library(dplyr)
male_tab = maleSigEdges %>%
  group_by(chr_TF) %>%
  summarize(mean_effect = mean(t, na.rm = TRUE))

female_tab = femaleSigEdges %>%
  group_by(chr_TF) %>%
  summarize(mean_effect = mean(t, na.rm = TRUE))

fulldata = data.frame(percentage=c(male_tab$mean_effect, female_tab$mean_effect))
fulldata$chromosome = factor(c(male_tab$chr_TF, female_tab$chr_TF))
fulldata$sex = factor(c(rep("Male", nrow(male_tab)), rep("Female", nrow(female_tab))), levels = c("Female", "Male"))
head(fulldata)

# barplot to show position of most significant TFs
library(ggplot2)
maintitle= paste("Mean effect sizes of chromosomes corresponding to \ntop", ng, "edges connected to significant aging genes", sep = " ")
ggplot(data=fulldata, aes(y=chromosome, x=percentage, fill=sex)) +
  geom_bar(stat="identity", position = "dodge")+theme_minimal() + ggtitle(maintitle) +xlab("mean effect size")


