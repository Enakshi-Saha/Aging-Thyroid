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

bonobo = data.frame(fread("/home/esaha/Aging_thyroid/networks/BonoboPanda_GSE165724.txt"))
genes = bonobo$V1
bonobo = bonobo[,-1]
bonobo_indices =  colnames(bonobo)

# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)

# Discard Y genes
### Get GenecIDs for Y genes 
Y_genes.GTEx = fdata.GTEx$gene_name[which(fdata.GTEx$seqnames == "chrY")]

bonobo = bonobo[-which(genes %in% Y_genes.GTEx),]
genes = genes[-which(genes %in% Y_genes.GTEx)]

# Get phenotypic data
phenotypes <- read.csv("/home/esaha/Aging_thyroid/data/validation_data/GSE165724/GSE165724_phenotypes_normalThyroid.txt", sep="")

age = as.numeric(phenotypes$age)
gender = factor(phenotypes$gender, levels = c("M", "F"))
tissue = factor(phenotypes$tissue_type, levels = c("Normal", "Normal_adjacent"))

###### Fit limma model ######
indegree = bonobo
rownames(indegree) = genes
#### Design Matrix
design = model.matrix(~ age*tissue + gender + tissue)

# Male-only model
# indegree = indegree[,which(gender == "M")]
# age = age[which(gender == "M")]
# tissue = tissue[which(gender == "M")]
# design = model.matrix(~age*tissue)

fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for age effect
# tb = topTable(fit,coef="age:genderF",number=Inf)
tb = topTable(fit,coef="age",number=Inf)
colnames(design) = make.names(colnames(design))
contrast.matrix <- makeContrasts(ageF = age + age.genderF, levels=design)
colnames(fit$coefficients) = rownames(contrast.matrix)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# tb = topTable(fit2,coef="ageF",number=Inf)

# Add chromosome information
tb$chr = fdata.GTEx$seqnames[match(rownames(tb),fdata.GTEx$gene_name)]
tb$gene_name = rownames(tb)
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

# save(fgseaRes, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE165724_Thyroid_Aging_SexDiff_agingSexInteraction_GSEA.RData")
# save(fgseaRes, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE165724_Thyroid_Aging_SexDiff_agingMALE_GSEA.RData")
# save(fgseaRes, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE165724_Thyroid_Aging_SexDiff_agingFEMALE_GSEA.RData")
# save(fgseaRes, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE165724_normalThyroid_Aging_GSEA.RData")
# save(fgseaRes, file = "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE165724_normalThyroid_Aging_GSEA_maleOnlyModel.RData")

# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_aging_sex_interaction.txt")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_aging_male.txt")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_aging_female.txt")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_normalThyroid_aging.txt")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_normalThyroid_aging_maleOnlyModel.txt")

# Save table for difference between normal and normal adjacent
tb = topTable(fit,coef="tissueNormal_adjacent",number=Inf)
# Add chromosome information
tb$chr = fdata.GTEx$seqnames[match(rownames(tb),fdata.GTEx$gene_name)]
tb$gene_name = rownames(tb)
head(tb)

# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_normalVSnormalAdjacent.txt")
# write.table(tb, "/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/limma_GSE165724_normalVSnormalAdjacent_maleOnlyModel.txt")

