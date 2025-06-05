# 
remotes::install_github("grst/immunedeconv")

library(recount)
library(recount3)
library(immunedeconv)

# Gene expression should be TPM, not log transformed, row names should be HGNC symbols
GTEx_tpm = read.delim("~/BONOBO/data/GTEx_thyroid/GTEx_thyroid_allSex.txt")
rownames(GTEx_tpm) = GTEx_tpm$gene
GTEx_tpm = GTEx_tpm[,-1]
# Transform data to non log scale
GTEx_tpm = 2^(GTEx_tpm)-1

info.GTEx = get(load("~/BONOBO/data/GTEx_thyroid/GTEx_thyroid_RSE.RData"))
gene_info = data.frame(info.GTEx@rowRanges)
gene_id = gsub("\\..*", "", gene_info$gene_id)
gene_names = gene_info$gene_name[match(rownames(GTEx_tpm), gene_id)]

# For rows with same gene name, take the mean as expression
GTEx_TPM = aggregate(GTEx_tpm, by=list(Category=gene_names), FUN=mean)
rownames(GTEx_TPM) = GTEx_TPM$Category
GTEx_TPM = GTEx_TPM[,-1]

xcell_gtex <- data.frame(immunedeconv::deconvolute(GTEx_TPM, "xcell"))

write.table(xcell_gtex, "/home/esaha/Aging_thyroid/data/immuneInfiltration/GTEx_xcell.txt", row.names = T, col.names = T)

