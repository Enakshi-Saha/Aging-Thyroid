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


####### GTEx disease pathways ##########
GSE27155_gsea = get(load("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE27155/sexAgnosticGSEA/GSE27155_diseaseEffect_Papillary Thyroid Carcinoma_GSEA.RData"))

####### GSE33630 disease pathways ##########
GSE33630_gsea = get(load("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE33630/sexAgnosticGSEA/GSE33630_diseaseEffect_papillary thyroid carcinoma (PTC)_GSEA.RData"))

#### Validated pathways in female (males can't be validated) ####
GSE27155_sig_paths = GSE27155_gsea$pathway[which(GSE27155_gsea$padj<0.05)]
GSE33630_sig_paths = GSE33630_gsea$pathway[which(GSE33630_gsea$padj<0.05)]

# pathways significant in larger data
common_paths = union(GSE33630_sig_paths, GSE27155_gsea)

# check if direction of effect is the same for both datasets (i.e. the product of NES is positive)
effect_product = GSE27155_gsea$NES[match(common_paths,GSE27155_gsea$pathway)]*
  GSE33630_gsea$NES[match(common_paths,GSE33630_gsea$pathway)]

# validated paths
validated_paths = unlist(common_paths[which(effect_product>0)])

PTC_info = list()
PTC_info[["dataset"]] = "GSE27155"
PTC_info[["validated_pathways"]] = validated_paths

save(PTC_info, file = "/home/esaha/Aging_thyroid/results/disease_specific_validated_pathways_list/PTC_pathways.RData.")

