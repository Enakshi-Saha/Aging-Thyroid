######### Heatmap for  effect #########
library(ggplot2)
# Load GSEA tables
GTEx_age = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_Aging_SexDiff_ageEffect_withoutSexInteraction_GSEA.RData"))
GTEx_below50 = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_ageEffect_female_below50_GSEA.RData"))
GTEx_above50 = get(load("/home/esaha/Aging_thyroid/results/limma_GSEA_xcell_tables/GTEx_Thyroid_ageEffect_female_above50_GSEA.RData"))

tab = data.frame(GTEx_age$NES)
rownames(tab) = GTEx_age$pathway
tab$below50 = GTEx_below50$NES[match(rownames(tab), GTEx_below50$pathway)]
tab$above50 = GTEx_above50$NES[match(rownames(tab), GTEx_above50$pathway)]
colnames(tab) = c("age effect", "age effect (age<=50)", "age effect (age>50)")
head(tab)

tab_pval = data.frame(GTEx_age$padj)
rownames(tab_pval) = GTEx_age$pathway
tab_pval$below50 = GTEx_below50$padj[match(rownames(tab_pval), GTEx_below50$pathway)]
tab_pval$above50 = GTEx_above50$padj[match(rownames(tab_pval), GTEx_above50$pathway)]
colnames(tab_pval) = c("age effect", "age effect (age<=50)", "age effect (age>50)")
head(tab_pval)

# Heatmap with significant pathways
sig_pathways = unique(c(GTEx_below50$pathway[which(GTEx_below50$padj <0.05)], GTEx_above50$pathway[which(GTEx_above50$padj <0.05)]))

tab_subset = tab[match(sig_pathways, rownames(tab)),-1]
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))
head(tab_subset)

# A matrix of characters (stars) to indicate which results are significant
lab = matrix(NA, nrow(tab_pval), ncol(tab_pval))
lab[which(tab_pval < 0.05)] = "*"
lab[which(tab_pval >= 0.05)] = ""
lab_subset = lab[match(sig_pathways, rownames(tab_pval)),-1]

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)


mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf(paste("/home/esaha/Aging_thyroid/results/plots/heatmap/GTEx_thyroid_ageEffect_female_byAgeGroup_heatmap.pdf", sep = ""),h=25,w=16)
heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
          cexRow=1, cexCol=1, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", cellnote=lab_subset, notecol="black", notecex=5)
dev.off()
