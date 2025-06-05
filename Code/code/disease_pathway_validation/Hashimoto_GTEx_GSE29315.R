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
GTEx_gsea = get(load("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_hashimoto_GSEA.RData"))

####### GSE29315 disease pathways ##########
GSE29315_gsea = get(load("~/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE29315/sexAgnosticGSEA/GSE29315_diseaseEffect_hashimoto thyroiditis_GSEA.RData"))

#### Validated pathways in female (males can't be validated) ####
GTEx_sig_paths = GTEx_gsea$pathway[which(GTEx_gsea$padj<0.05)]
GSE29315_sig_paths = GSE29315_gsea$pathway[which(GSE29315_gsea$padj<0.05)]

# pathways significant in both data
common_paths = union(GTEx_sig_paths, GSE29315_sig_paths)

# check if direction of effect is the same for both datasets (i.e. the product of NES is positive)
effect_product = GTEx_gsea$NES[match(common_paths,GTEx_gsea$pathway)]*
  GSE29315_gsea$NES[match(common_paths,GSE29315_gsea$pathway)]

# validated paths
validated_paths = common_paths[which(effect_product>0)]

hashimoto_info = list()
hashimoto_info[["dataset"]] = "GTEx"
hashimoto_info[["validated_pathways"]] = validated_paths

save(hashimoto_info, file = "/home/esaha/Aging_thyroid/results/disease_specific_validated_pathways_list/hashimoto_validated_pathways.RData.")

