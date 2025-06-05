library(ggplot2)
# Load GSEA tables 
GTEx = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_hashimotoExcluded_samples/GTEx_Thyroid_aging_GSEA.RData"))
GSE165724 = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GSE165724_normalThyroid_Aging_GSEA.RData"))

tab = data.frame(GTEx$NES)
rownames(tab) = GTEx$pathway
tab$GSE165724 = GSE165724$NES[match(rownames(tab), GSE165724$pathway)]
colnames(tab) = c("GTEx", "GSE165724")
head(tab)

tab_pval = data.frame(GTEx$padj)
rownames(tab_pval) = GTEx$pathway
tab_pval$GSE165724 = GSE165724$padj[match(rownames(tab_pval), GSE165724$pathway)]
colnames(tab_pval) = c("GTEx", "GSE165724")
head(tab_pval)

# Heatmap with significant pathways
sig_pathways = GTEx$pathway[which(GTEx$padj <0.05)]

tab_subset = tab[match(sig_pathways, rownames(tab)),]
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(tab_subset)

# A matrix of characters (stars) to indicate which results are significant
lab = matrix(NA, nrow(tab_pval), ncol(tab_pval))
lab[which(tab_pval < 0.05)] = "*"
lab[which(tab_pval >= 0.05)] = ""
lab_subset = lab[match(sig_pathways, rownames(tab_pval)),]

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)


mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_thyroid_hashimotoExcluded/GTEx_GSE165724normal_thyroid_aging_validation_noStar.pdf",h=25,w=16)
heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1.5, cexCol=1.5, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", notecol="black", notecex=5)
dev.off()

