library(ggplot2)
# Load GSEA tables for MALEaging, FEMALaging and SexAgeInteraction
interaction = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_Aging_SexDiff_agingSexInteraction_GSEA.RData"))

male = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_Aging_SexDiff_agingMALE_GSEA.RData"))
female = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_Aging_SexDiff_agingFEMALE_GSEA.RData"))

tab = data.frame(interaction$NES)
rownames(tab) = interaction$pathway
tab$male = male$NES[match(rownames(tab), male$pathway)]
tab$female = female$NES[match(rownames(tab), female$pathway)]
colnames(tab) = c("Age effect: female-male", "Age effect: male","Age effect: female")
head(tab)

tab_pval = data.frame(interaction$padj)
rownames(tab_pval) = interaction$pathway
tab_pval$male = male$padj[match(rownames(tab_pval), male$pathway)]
tab_pval$female = female$padj[match(rownames(tab_pval), female$pathway)]
colnames(tab_pval) = c("Age effect: female-male", "Age effect: male","Age effect: female")
head(tab_pval)

# Heatmap with significant pathways
sigPath_interaction = interaction$pathway[which(interaction$padj <0.05)]
sigPath_male = male$pathway[which(male$padj <0.05)]
sigPath_female = female$pathway[which(female$padj <0.05)]

sig_pathways = unique(c(sigPath_female, sigPath_male))
# check if any pathway has significant interaction but not significant in either sex
sigPath_interaction[-which(sigPath_interaction %in% sig_pathways)]

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
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_thyroid_aging_SexDiff_heatmap.pdf",h=25,w=16)
heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1.5, cexCol=1.5, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", cellnote=lab_subset, notecol="black", notecex=5)
dev.off()

##########Different plots for sex and interaction ###########
mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_thyroid_aging_SexDiff_heatmap_noInteractionColumn.pdf",h=25,w=16)
heatmap.2(as.matrix(tab_subset[,1:2]),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1.5, cexCol=1.5, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", cellnote=lab_subset, notecol="black", notecex=5)
dev.off()



