library(ggplot2)
# Load GSEA tables
GTEx_age = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_Aging_SexDiff_ageEffect_withoutSexInteraction_GSEA.RData"))
GTEx_sex = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_SexDiff_FEMALE.RData"))

tab = data.frame(GTEx_age$NES)
rownames(tab) = GTEx_age$pathway
tab$GTEx_sex = GTEx_sex$NES[match(rownames(tab), GTEx_sex$pathway)]
colnames(tab) = c("age effect", "sex effect (F-M)")
head(tab)

tab_pval = data.frame(GTEx_age$padj)
rownames(tab_pval) = GTEx_age$pathway
tab_pval$GTEx_sex = GTEx_sex$padj[match(rownames(tab_pval), GTEx_sex$pathway)]
colnames(tab_pval) = c("age_effect", "sex effect (F-M)")
head(tab_pval)

# Heatmap with significant pathways
sig_pathways = GTEx_age$pathway[which(GTEx_age$padj <0.05)]

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
pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_thyroid_age_sex_heatmap.pdf",h=25,w=16)
# pdf("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_GSE165724normal_thyroid_aging_validation_maleOnlyModel.pdf",h=25,w=16)
heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1.5, cexCol=1.5, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", cellnote=lab_subset, notecol="black", notecex=5)
dev.off()

