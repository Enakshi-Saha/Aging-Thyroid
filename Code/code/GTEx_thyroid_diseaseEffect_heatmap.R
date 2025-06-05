######### Heatmap for disease effect #########

library(ggplot2)
# Load GSEA tables
GTEx_age = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_Aging_SexDiff_ageEffect_withoutSexInteraction_GSEA.RData"))
GTEx_hashimoto = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_hashimoto_GSEA.RData"))
GTEx_goiter = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_goiter_GSEA.RData"))
GTEx_adenoma = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_diseaseEffect_adenoma_GSEA.RData"))

tab = data.frame(GTEx_age$NES)
rownames(tab) = GTEx_age$pathway
tab$hashimoto = GTEx_hashimoto$NES[match(rownames(tab), GTEx_hashimoto$pathway)]
tab$goiter = GTEx_goiter$NES[match(rownames(tab), GTEx_goiter$pathway)]
tab$adenoma = GTEx_adenoma$NES[match(rownames(tab), GTEx_adenoma$pathway)]
colnames(tab) = c("age effect", "hashimoto effect", "goiter effect", "adenoma")
head(tab)

tab_pval = data.frame(GTEx_age$padj)
rownames(tab_pval) = GTEx_age$pathway
tab_pval$hashimoto = GTEx_hashimoto$padj[match(rownames(tab_pval), GTEx_hashimoto$pathway)]
tab_pval$goiter = GTEx_goiter$padj[match(rownames(tab_pval), GTEx_goiter$pathway)]
tab_pval$adenoma = GTEx_adenoma$padj[match(rownames(tab_pval), GTEx_adenoma$pathway)]
colnames(tab_pval) = c("age effect", "hashimoto effect", "goiter effect", "adenoma")
head(tab_pval)

# Heatmap with significant pathways
sig_pathways_h = c(GTEx_age$pathway[which(GTEx_age$padj <0.05)], GTEx_hashimoto$pathway[which(GTEx_hashimoto$padj <0.05)])
sig_pathways_g = c(GTEx_age$pathway[which(GTEx_age$padj <0.05)], GTEx_goiter$pathway[which(GTEx_goiter$padj <0.05)])
sig_pathways_a = c(GTEx_age$pathway[which(GTEx_age$padj <0.05)], GTEx_adenoma$pathway[which(GTEx_adenoma$padj <0.05)])

tab_subset_h = tab[match(sig_pathways_h, rownames(tab)),1:2]
rownames(tab_subset_h) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset_h), split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(tab_subset_h)

tab_subset_g = tab[match(sig_pathways_g, rownames(tab)),c(1,3)]
rownames(tab_subset_g) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset_g), split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(tab_subset_g)

tab_subset_a = tab[match(sig_pathways_a, rownames(tab)),c(1,4)]
rownames(tab_subset_a) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset_a), split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(tab_subset_a)

# A matrix of characters (stars) to indicate which results are significant
lab = matrix(NA, nrow(tab_pval), ncol(tab_pval))
lab[which(tab_pval < 0.05)] = "*"
lab[which(tab_pval >= 0.05)] = ""
lab_subset_h = lab[match(sig_pathways_h, rownames(tab_pval)),1:2]
lab_subset_g = lab[match(sig_pathways_g, rownames(tab_pval)),c(1,3)]
lab_subset_a = lab[match(sig_pathways_a, rownames(tab_pval)),c(1,3)]

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)


mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_thyroid_diseaseEffect_hashimoto_heatmap.pdf",h=25,w=16)
heatmap.2(as.matrix(tab_subset_h),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1, cexCol=1, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", cellnote=lab_subset_h, notecol="black", notecex=5)
dev.off()

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_thyroid_diseaseEffect_goiter_heatmap.pdf",h=25,w=16)
heatmap.2(as.matrix(tab_subset_g),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1, cexCol=1, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", cellnote=lab_subset_g, notecol="black", notecex=5)
dev.off()

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_thyroid_diseaseEffect_adenoma_heatmap.pdf",h=25,w=16)
heatmap.2(as.matrix(tab_subset_a),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1, cexCol=1, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", cellnote=lab_subset_a, notecol="black", notecex=5)
dev.off()
