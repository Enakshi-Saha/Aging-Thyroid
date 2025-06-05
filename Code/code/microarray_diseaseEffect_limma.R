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


# Get gene data
GTEx_thyroid = get(load("/home/esaha/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
fdata.GTEx = as.data.frame(GTEx_thyroid@rowRanges)
gene_names = fdata.GTEx$gene_name

##############################
# Set datasetID to be one of these three: GSE29315, GSE33630, GSE27155
datasetID = "GSE27155"
setwd("/home/esaha/Aging_thyroid/networks/")
file=paste("BonoboPanda_", datasetID, ".txt", sep = "")

bonobo = data.frame(fread(file))
genes = bonobo$V1
bonobo = bonobo[,-1]
bonobo_indices =  colnames(bonobo)
rownames(bonobo) = genes

# Get phenotypic data
phenotypes <- read.csv(paste("~/Aging_thyroid/data/validation_data/",datasetID, "/", datasetID, "_phenotypes.txt", sep = ""), sep="")
phenotypes = phenotypes[match(bonobo_indices, rownames(phenotypes)),]

gender = factor(phenotypes$sex)
tissue = factor(phenotypes$tissue_type)

normal_tissue_labels = c("hyperplasia", "Normal Thyroid", "patient-matched non-tumor control")

# Set normal thyroid as reference level for tissue
reference_tissue = levels(tissue)[which(levels(tissue) %in% normal_tissue_labels)]
tissue = relevel(tissue, ref = reference_tissue)

#### Design Matrix
design = model.matrix(~ gender + tissue)
fit <- lmFit(bonobo, design)
fit <- eBayes(fit)

# Save table for disease
diseases = levels(tissue)[-1]
for (d in diseases){
  coef_name = paste("tissue",d,sep="")
  tb = topTable(fit,coef=coef_name,number=Inf)
  tb$chr = fdata.GTEx$seqnames[match(rownames(tb),gene_names)]
  tb$gene_name = fdata.GTEx$gene_name[match(rownames(tb),gene_names)]
  
  # save limma table
  write.table(tb, file = paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/", datasetID, "/limma_", datasetID, "_", "diseaseEffect_", d, ".txt",sep = ""))
  
  # Rank genes in limma table
  indegree_rank <- setNames(object=tb[,"t"], tb$gene_name)
  # Load KEGG pathways
  pathways <- gmtPathways("/home/ubuntu/c2.cp.kegg.v7.1.symbols.gmt")
  fgseaRes <- fgsea(pathways, indegree_rank, minSize=15, maxSize=500)
  
  # save GSEa table
  save(fgseaRes, file = paste("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/", datasetID, "/", datasetID, "_", "diseaseEffect_", d, "_GSEA.RData",sep = ""))
}
  



